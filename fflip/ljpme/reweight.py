# -*- coding: utf-8 -*-

import warnings
# from fflip.omm.reweightprop import *
from fflip.omm.util import check_and_make_dir
from fflip.ljpme.reweightfuncs import *
from fflip.ljpme.rdfhandle import *
from fflip.ljpme.util import *
from fflip.ljpme.scheme import LipidScheme


class ReweightTarget(object):
    def __init__(self, name, temperature, property_dir, energy_dir,
                 property_file_template, result_dir,
                 exp, lipid, exp_dir, parse_groups,
                 sim_rdf_r_range=(0.0015, 0.999), sim_rdf_r_intvl=0.003):
        """
        Args:
            name:
            temperature:
            property_dir:
            energy_dir:
            property_file_template:
            result_dir:
            exp:
            lipid:
            exp_dir:
            sim_rdf_r_range:
            sim_rdf_r_intvl:
        """
        self.name = name
        self.temperature = temperature
        self.property_dir = property_dir
        self.energy_dir = energy_dir
        self.property_file_template = property_file_template
        self.result_dir = result_dir
        self.exp = exp
        self.exp_dir = exp_dir
        self.sim = None
        self.rew = None
        self.lipid = lipid
        self.lipid_scheme = LipidScheme(lipid)
        self.groups = parse_groups
        self.sim_rdf_r_range = sim_rdf_r_range
        self.sim_rdf_r_intvl = sim_rdf_r_intvl
        if not os.path.isdir(self.result_dir):
            os.system("mkdir -p {}".format(self.result_dir))

    @property
    def done_reweighting(self):
        abs_path_for_rew_result = os.path.abspath(self.result_dir)
        if os.path.isfile(
                os.path.join(abs_path_for_rew_result, self.name + '.rew')
        ) or type(self.rew) == np.ndarray:
            return True
        else:
            return False

    @property
    def ngroups(self):
        """
        Get the numebr of perturbed parameters
        """
        if self.groups == 'all':
            return len(self.lipid_scheme.parse())
        elif isinstance(self.groups, list):
            return len(
                self.lipid_scheme.parse(groups=self.groups)
            )
        else:
            raise Exception

    @property
    def property_type(self):
        if 'rdf' in str(self.property_dir):
            return 2
        else:
            return 1

    @property
    def sim_x(self):
        if self.property_type == 2:
            return np.arange(
                self.sim_rdf_r_range[0],
                self.sim_rdf_r_range[1],
                self.sim_rdf_r_intvl
            )
        else:
            return None

    @property
    def exp_x(self):
        # used for rdf
        if self.property_type == 2:
            # TODO: the file name should be provided by the user, not fixed
            return np.loadtxt(os.path.join(self.exp_dir, 'r.exp'))
        else:
            return None

    def reweight(self, starting_traj, ending_traj,
                 traj_interval_energy, traj_interval_prop=None,
                 use_cluster=False, partition='ivy,sbr',
                 **kwargs):
        param_ids = self.lipid_scheme.param_ids(**kwargs)
        if not use_cluster:
            original, perturbed = reweight_many_params(
                param_ids, self.temperature,
                os.path.join(
                    self.property_dir, 'block_data', self.property_file_template
                ),
                os.path.join(self.energy_dir, 'block_data/original_{}.dat'),
                os.path.join(self.energy_dir, 'block_data/perturbed_{}_{}.dat'),
                starting_traj, ending_traj,
                traj_interval_energy, traj_interval_prop
            )
            self.sim = np.array(original)
            self.rew = np.array(perturbed)
        else:
            from dask.distributed import Client
            from dask_jobqueue import SLURMCluster
            cluster = SLURMCluster(
                queue=partition,
                cores=1,
                walltime="00:30:00",
                shebang='#!/usr/bin/bash',
                memory="4 GB"
            )
            cluster.scale(jobs=1)
            client = Client(cluster)
            future_result = client.submit(
                reweight_many_params,
                param_ids, self.temperature,
                os.path.join(
                    self.property_dir, 'block_data', self.property_file_template
                ),
                os.path.join(self.energy_dir, 'block_data/original_{}.dat'),
                os.path.join(self.energy_dir, 'block_data/perturbed_{}_{}.dat'),
                starting_traj, ending_traj,
                traj_interval_energy, traj_interval_prop
            )
            self.sim = future_result.result()[0]
            self.rew = future_result.result()[1]

    def save_reweighted(self):
        check_and_make_dir(os.path.abspath(self.result_dir))
        to_save = np.atleast_1d(self.sim)
        np.savetxt(os.path.join(self.result_dir, self.name + '.sim'), to_save)
        to_save = np.atleast_1d(self.rew)
        np.savetxt(os.path.join(self.result_dir, self.name + '.rew'), to_save)

    def robustness_analysis(self, first_trj, last_trj,
                            trj_interval_energy, trj_interval_prop=None,
                            block_size=50, use_cluster=True,
                            partition='ivy,sbr',
                            **kwargs):
        org = []
        ptb = []
        diff = []
        if not use_cluster:
            warnings.warn("This part of code is not up to date, use cluster!")
            param_ids = self.lipid_scheme.param_ids(**kwargs)
            for starting_trj in range(first_trj, last_trj, block_size):
                original, perturbed = reweight_many_params(
                    param_ids, self.temperature,
                    os.path.join(
                        self.property_dir, 'block_data',
                        self.property_file_template
                    ),
                    os.path.join(self.energy_dir, 'block_data/original_{}.dat'),
                    os.path.join(
                        self.energy_dir, 'block_data/perturbed_{}_{}.dat'
                    ),
                    starting_trj, starting_trj + block_size - 1,
                    trj_interval_energy, trj_interval_prop
                )
                org.append(original)
                ptb.append(perturbed)
                diff.append(perturbed - original)
            # self.robustness = (np.mean(np.array(diff), axis = 0) / np.std(
            # np.array(diff), axis = 0))**2
            return (np.mean(np.array(diff), axis=0) /
                    np.std(np.array(diff), axis=0)) ** 2
        else:
            from dask.distributed import Client
            from dask_jobqueue import SLURMCluster
            dask_futures = []
            param_ids = self.lipid_scheme.param_ids(**kwargs)
            for starting_trj in range(first_trj, last_trj, block_size):
                cluster = SLURMCluster(
                    queue=partition,
                    cores=1,
                    walltime="00:30:00",
                    shebang='#!/usr/bin/bash',
                    memory="4 GB"
                )
                cluster.scale(jobs=1)
                client = Client(cluster)
                dask_futures.append(
                    client.submit(
                        reweight_many_params,
                        param_ids, self.temperature,
                        os.path.join(
                            self.property_dir, 'block_data',
                            self.property_file_template
                        ),
                        os.path.join(
                            self.energy_dir, 'block_data/original_{}.dat'
                        ),
                        os.path.join(
                            self.energy_dir, 'block_data/perturbed_{}_{}.dat'
                        ),
                        starting_trj, starting_trj + block_size - 1,
                        trj_interval_energy, trj_interval_prop
                    )
                )
                # this is ugly, change later ...
            for future_result in dask_futures:
                original = future_result.result()[0]
                perturbed = future_result.result()[1]
                if not ('area' in self.name or 'scd' in self.name):
                    evaluator = SensitivityEvaluator(
                        self.ngroups, self.exp_x, self.exp, self.sim_x,
                        original, perturbed, sens_type=self.property_type,
                        n_peaks=2
                    )
                    one_diff = evaluator.diff_rew_sim
                    diff.append(one_diff)
                else:
                    # we can also use the sensitivity_evaluator here!
                    diff.append(perturbed - original)
            return diff

    def retrieve_sim_rew(self):
        if not self.sim:
            try:
                self.sim = np.loadtxt(
                    os.path.join(self.result_dir, self.name + '.sim')
                )
            except NoSimDataError:
                return 0
        if not self.rew:
            try:
                self.rew = np.loadtxt(
                    os.path.join(self.result_dir, self.name + '.rew')
                )
            except NoRewDataError:
                return 0
        return 1

    def add_sensitivity_evaluator(self):
        """This method doesn't check if the reweighting is done or not"""
        if not (type(self.sim) == np.ndarray and type(self.rew) == np.ndarray):
            self.retrieve_sim_rew()
        else:
            self.sim = np.array(self.sim)
            self.rew = np.array(self.rew)
        if self.name == 'O2-OW' or self.name == 'Ob-OW':
            evaluator = SensitivityEvaluator(
                self.ngroups, self.exp_x, self.exp, self.sim_x, self.sim,
                self.rew, sens_type=self.property_type, n_peaks=1
            )
        else:
            evaluator = SensitivityEvaluator(
                self.ngroups, self.exp_x, self.exp, self.sim_x, self.sim,
                self.rew, sens_type=self.property_type
            )
        self.sensitivity_evaluator = evaluator

    @property
    def diff_sim_exp(self):
        # this might be wrong?
        if self.sensitivity_evaluator:
            return self.sensitivity_evaluator.diff_sim_exp

    @property
    def rel_diff_sim_exp(self):
        if self.sensitivity_evaluator:
            return self.sensitivity_evaluator.rel_diff_sim_exp

    @property
    def diff_rew_sim(self):
        if self.sensitivity_evaluator:
            return self.sensitivity_evaluator.diff_rew_sim

    @property
    def rel_diff_rew_sim(self):
        if self.sensitivity_evaluator:
            return self.sensitivity_evaluator.rel_diff_rew_sim


class SensitivityEvaluator(object):
    def __init__(self, ngroups, exp_x, exp, sim_x, sim, rew,
                 sens_type=1, n_peaks=2, n_foots=1):
        """
        Args:
            -- exp: the experimental value(s)
            -- sim: the simulated value(s)
            -- rew: the reweighted value(s)
            -- sens_type: the catagory of the property (1: area/scd/db, 2: rdf)
        """
        self.ngroups = ngroups
        self.sim = sim
        self.exp = exp
        if sens_type == 2:
            self.exp_x = exp_x
            self.sim_x = sim_x
        self.rew = rew
        self.sens_type = sens_type
        self.n_peaks = n_peaks
        self.n_foots = n_foots

    @property
    def diff_sim_exp(self):
        if self.sens_type == 1:
            """
            Area / Order parameter / DB (thickness)
            """
            return self.sim - self.exp
        elif self.sens_type == 2:
            """
            Things like rdf which contain both positions and magnitudes
            """
            x_exp, y_exp = find_rdf_peaks_and_foots(
                self.exp, self.exp_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=3
            )
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim, self.sim_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=1
            )
            return x_sim - x_exp, y_sim - y_exp

    @property
    def rel_diff_sim_exp(self):
        if self.sens_type == 1:
            """
            Area / Order parameter / DB
            """
            return (self.sim - self.exp) / self.exp
        elif self.sens_type == 2:
            """
            Things like rdf which contain both positions and magnitudes
            """
            x_exp, y_exp = find_rdf_peaks_and_foots(
                self.exp, self.exp_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=3
            )
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim, self.sim_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=1
            )
            return (x_sim - x_exp)/x_exp, (y_sim - y_exp)/y_exp

    @property
    def diff_rew_sim(self):
        if self.sens_type == 1:
            rew = self.rew
            sim_tiled = np.tile(self.sim, self.ngroups)
            return rew - sim_tiled
        elif self.sens_type == 2:
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim, self.sim_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=1
            )
            r_list, peak_foot_value_list = find_rdf_peaks_and_foots(
                self.rew, self.sim_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=1
            )
            diff_list = []
            for i, (r, pfv) in enumerate(zip(r_list, peak_foot_value_list)):
                diff_list.append(np.array([r - x_sim, pfv - y_sim]))
            return np.array(diff_list)

    @property
    def rel_diff_rew_sim(self):
        if self.sens_type == 1:
            rew = self.rew
            sim_tiled = np.tile(self.sim, self.ngroups)
            return (rew - sim_tiled) / self.exp
        elif self.sens_type == 2:
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim, self.sim_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=1
            )
            r_list, peak_foot_value_list = find_rdf_peaks_and_foots(
                self.rew, self.sim_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=1
            )
            x_exp, y_exp = find_rdf_peaks_and_foots(
                self.exp, self.exp_x, first_n_peaks=self.n_peaks,
                first_n_foots=self.n_foots, smooth_window_size=3
            )
            rel_diff_list = []
            for i, (r, pfv) in enumerate(zip(r_list, peak_foot_value_list)):
                rel_diff_list.append(
                    np.array([(r - x_sim) / x_exp, (pfv - y_sim) / y_exp])
                )
            return np.array(rel_diff_list)

            







