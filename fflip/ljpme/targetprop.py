# -*- coding: utf-8 -*-

# Contains master class of training target (area, rdf, scd, db ...)

from fflip.ljpme.reweight import *
from fflip.omm.genclac import OmmJobGenerator
from fflip.ljpme.scheme import *
from fflip.ljpme.param import gen_param_offset


class KaGenerator(object):
    def __init__(self, target_properties, **kwargs):
        """
        Args:
            properties: a list contains at least two TargetProperty object
        Returns:
            the deviation and sensitivity
        """
        if 'delta' in kwargs:
            self.deltagamma = kwargs['delta']
        else:
            self.deltagamma = 10
        self.target_properties = target_properties

    def gen_sim(self):
        delta_a = self.target_properties[1].reweight_target.sim - \
                  self.target_properties[2].reweight_target.sim
        a_0 = self.target_properties[0].reweight_target.sim
        return 2 * a_0 * self.deltagamma / delta_a

    def gen_rew(self):
        delta_a = self.target_properties[1].reweight_target.rew - \
                 self.target_properties[2].reweight_target.rew
        a_0 = self.target_properties[0].reweight_target.rew
        return 2 * a_0 * self.deltagamma / delta_a

    def gen_sensitivity(self):
        a_0 = self.target_properties[0].reweight_target.sim
        delta_a = self.target_properties[1].reweight_target.sim - \
            self.target_properties[2].reweight_target.sim
        sens_a0 = self.target_properties[0].sensitivity
        sens_delta_a = self.target_properties[1].sensitivity - \
            self.target_properties[2].sensitivity
        return 2 * self.deltagamma * (
                sens_a0 / delta_a - a_0 * sens_delta_a / delta_a ** 2
        )


class DeltaGenerator(object):
    def __init__(self, target_properties, **kwargs):
        """
        Currently used for temperature dependence
        Args:
            target_properties:  a list contains two TargetProperty object
            **kwargs:
        Returns
            the deviation and sensitivity
        """
        self.target_properties = target_properties

    def gen_sim(self):
        delta = self.target_properties[1].reweight_target.sim - \
            self.target_properties[0].reweight_target.sim
        return delta

    def gen_rew(self):
        delta = self.target_properties[1].reweight_target.rew - \
            self.target_properties[0].reweight_target.rew
        return delta

    def gen_sensitivity(self):
        sens_delta = self.target_properties[1].sensitivity - \
            self.target_properties[0].sensitivity
        return sens_delta


class TargetSystem(object):
    def __init__(self,
                 system_type,
                 lipid_name,
                 num_lipids,
                 temperature,
                 surface_tension,
                 psf_file,
                 crd_file,
                 pot_template,
                 sim_template=None,
                 option_scheme=SimOptScheme,
                 naming_scheme=FolderNamingScheme,
                 ff='c36'):
        self.system_type = system_type
        self.lipid_name = lipid_name
        self.num_lipids = num_lipids
        self.temperature = temperature
        self.surface_tension = surface_tension
        self.psf_file = psf_file
        self.crd_file = crd_file
        self.pot_template = pot_template
        self.sim_template = sim_template
        self.option_scheme = option_scheme(self)
        self.folder_naming = naming_scheme(self)
        self.ff = ff

    def simulate(self, iteration, trj_folder, last_seqno=None,
                 boxx=None, boxz=None, zmode=None,
                 barostat=None, integrator=None,
                 change_para=False, solution_file=None, torfix_file=None,
                 overwrite=False, start=False, verbose=0):
        if self.sim_template is None:
            return 0  # exit

        if change_para:
            assert solution_file is not None  # and torfix_file is not None

        trj_loc = self.folder_naming.trajectory_folder(iteration, trj_folder)
        if verbose >= 1:
            print("Runnning in {}".format(trj_loc))
        calc = OmmJobGenerator(
            self.crd_file, self.psf_file,
            template=self.sim_template, work_dir=trj_loc
        )
        # constructing options dictionary
        options = dict()
        # special treatment to the option file name
        options["mdo"] = "input.mdo"
        options["surface_tension"] = self.surface_tension
        if boxx is not None:
            options["boxx"] = float(boxx)
        else:
            options["boxx"] = self.option_scheme.boxx
        if boxz is not None:
            options["boxz"] = boxz
        elif self.option_scheme.boxz is not None:
            # for systems, this can be empty
            options["boxz"] = self.option_scheme.boxz
        options["temperature"] = self.temperature
        if zmode is not None:
            options["zmode"] = int(zmode)
        elif self.option_scheme.zmode is not None:
            # for systems, this can be empty
            options["zmode"] = self.option_scheme.zmode
        if barostat is not None:
            options["barostat"] = barostat
        else:
            options["barostat"] = self.option_scheme.barostat
        if integrator is not None:
            options["intgrt"] = integrator
        else:
            options["intgrt"] = self.option_scheme.intgrt
        options["change_para"] = 'yes' if change_para else 'no'
        if change_para:
            options["sfile"] = solution_file
            options["tfile"] = torfix_file
        options["psf"] = self.psf_file
        options["crd"] = self.crd_file
        options["lipid"] = self.lipid_name
        if last_seqno is not None:
            options["last_seqno"] = int(last_seqno)
        else:
            options["last_seqno"] = self.option_scheme.critical_seqno[1]
        # get the job run
        job = calc(0, options, overwrite=overwrite)
        if start:
            job("rflow submit sdyn.sh")


class SpecialProperty(TargetSystem):
    def __init__(self, name, temperature, weight_factor, root_dir,
                 parent_properties, generator, exp_rel_dir='exp', **kwargs):
        """
        The (sensitivity) evaluator for compressibility of membrane, 
        based on differential of Ka = 2A0 * (d_gamma / dA)
        Args:
            name: name of the property
            temperature: temperature
            weight_factor: weight factor in solving best parameter set
            root_dir: root directory to do the optimization
            parent_properties: the A0, the A under minus surface tension,
                               and the A under positive surface tension.
            exp_rel_dir: experimental dat dir
            **kwargs: 
        """
        self.name = name
        self.temperature = temperature
        self.root_dir = root_dir
        self.exp_dir = os.path.join(root_dir, exp_rel_dir)
        self.parent_properties = parent_properties
        self.generator = generator
        self.options = kwargs
        self._scaling = make_guess_of_scaling(self.name)
        self._app_weight = weight_factor

    @property
    def exp(self):
        return np.loadtxt(os.path.join(self.exp_dir, self.name + '.exp'))

    def update_scaling_using_exp(self):
        self._scaling = 1 / self.exp
        pass

    @property
    def app_weight_factor(self):
        return self._app_weight

    @property
    def scaling(self):
        return self._scaling

    @property
    def weight_factor(self):
        return self.app_weight_factor * self.scaling

    def get_sensitivity(self, iteration):
        for p in self.parent_properties:
            p.get_sensitivity(iteration)
        generator = self.generator(self.parent_properties, **self.options)
        self.rew = generator.gen_rew()
        self.sim = generator.gen_sim()
        self.deviation = self.exp - self.sim
        self.sensitivity = generator.gen_sensitivity()


class TargetProperty(TargetSystem):
    def __init__(self,
                 name,
                 prop_type,
                 system_type,
                 lipid_name,
                 lipid,
                 num_lipids,
                 weight_factor,
                 temperature,
                 surface_tension,
                 psf_file,
                 crd_file,
                 root_dir,
                 pot_template,
                 obs_template,
                 obs_file_format,
                 sim_template=None,
                 sim_scheme=SimOptScheme,
                 naming_scheme=FolderNamingScheme,
                 ff='c36',
                 parse_groups='all',
                 **misc_kwargs):
        """
        misc_kwargs may contain:
        exp_rel_dir, ... (more to come)
        """
        self.name = name
        self.prop_type = prop_type
        self.system_type = system_type
        self.lipid_name = lipid_name
        self.num_lipids = num_lipids
        self.temperature = temperature
        self.surface_tension = surface_tension
        self.psf_file = psf_file
        self.crd_file = crd_file
        self.root_dir = root_dir
        self.obs_template = obs_template
        self.pot_template = pot_template
        self.sim_template = sim_template
        self.property_file_format = obs_file_format
        self.option_scheme = sim_scheme(self)
        self.folder_naming = naming_scheme(self)
        if not 'exp_rel_dir' in misc_kwargs:
            self.exp_dir = self.folder_naming.exp_folder()
        else:
            self.exp_dir = self.folder_naming.exp_folder(
                rel_dir=misc_kwargs['exp_rel_dir']
            )
        self.reweight_dir = self.folder_naming.reweight_dir()
        self.robdir = self.folder_naming.robustness_dir()
        self.first_trj, self.last_trj = \
            make_guess_of_trajectory_range(self.name)
        # self.trj_intvl_e, self.trj_intvl_p = make_guess_of_intervals(self.name)
        self.lipid = lipid
        self.lipid_scheme = LipidScheme(self.lipid)
        self.ff = ff
        # CHARMM integer charge groups
        self.groups = parse_groups
        self._prop_block_size = make_guess_of_block_size(0, self.prop_type)
        self._pot_block_size = make_guess_of_block_size(1, self.prop_type)
        self._scaling = make_guess_of_scaling(self.name)
        self._app_weight = weight_factor
        self._perturbation = 1.0

    def update_block_sizes(self, **kwargs):
        if "observable" in kwargs:
            self._prop_block_size = kwargs["observable"]
        if "potential" in kwargs:
            self._pot_block_size = kwargs["potential"]

    def update_amount_of_perturbation(self, amount):
        self._perturbation = round(float(amount), 4)

    def update_first_last_trj(self, first_trj, last_trj):
        self.first_trj = first_trj
        self.last_trj = last_trj

    def update_app_weight(self, w):
        """ Update the apparent weight factor"""
        self._app_weight = w

    def update_scaling_using_exp(self):
        self._scaling = 1 / self.exp
        pass

    @property
    def prop_block_size(self):
        return self._prop_block_size

    @property
    def pot_block_size(self):
        return self._pot_block_size

    @property
    def perturbation(self):
        return self._perturbation

    @property
    def app_weight_factor(self):
        return self._app_weight

    @property
    def scaling(self):
        return self._scaling

    @property
    def weight_factor(self):
        return self.app_weight_factor * self.scaling

    def charge_perturbation(self):
        pass

    def parameter_offsets(self):
        parameter_sets = self.lipid.parse_groups(groups=self.groups)
        offsets = [
            gen_param_offset(
                ps, amount=self.perturbation
            ) for ps in parameter_sets
        ]
        return offsets

    @property
    def _exp(self):
        if 'peak' in self.name or 'foot' in self.name:
            name_exp = self.name.split('_')[0]
        else:
            name_exp = self.name
        return np.loadtxt(os.path.join(self.exp_dir, name_exp + '.exp'))
    
    @property
    def exp(self):
        if not hasattr(self, 'exchanged_exp'):
            return extract_exp(self.name, self._exp)
        else:
            return extract_exp(self.name, self.exchanged_exp)

    @property
    def parameters(self):
        if not isinstance(self.lipid, DrudeLipid):
            return self.lipid.parse_groups(groups=self.groups)
        else:
            return self.lipid.parse_groups(id_allowed=self.groups)

    @property
    def num_parameters(self):
        return len(self.parameters)

    def reverse_para(self):
        pass

    def simulation(self, iteration, trj_folder, last_seqno=None):
        super().simulate(
            self, iteration, trj_folder, last_seqno=last_seqno,
            change_para=False, solution_file=None, torfix_file=None
        )

    def recalc_energy(self, iteration, traj_root, toppar_path,
                      overwrite=False, wait=True,
                      last_solution=None, torfix=None):

        print("parameter perturbation: ~ {}".format(self.perturbation))

        trj_loc = self.folder_naming.trajectory_folder(iteration, traj_root)

        work_dir = self.folder_naming.potential_data_folder(iteration)

        calc = OmmJobGenerator(
            self.crd_file, self.psf_file,
            template=self.pot_template, work_dir=work_dir
        )

        options = dict()
        options["option_file"] = "potcalc.inp"
        options["name"] = self.lipid_name
        options["amount"] = self.perturbation
        options["trj_location"] = trj_loc
        options["first_trj"] = self.first_trj
        options["last_trj"] = self.last_trj
        options["block_size"] = self.pot_block_size
        options["toppar_path"] = toppar_path
        if last_solution is not None:
            assert os.path.isfile(last_solution)
            # os.system("cp {} {}".format(last_solution, work_dir))
        if torfix is not None:
            assert os.path.isfile(torfix)
            # os.system("cp {} {}".format(torfix, work_dir))
        options["solution"] = last_solution
        options["torfix"] = torfix
        job = calc(1, options, overwrite=overwrite)
        job('python submit.py')
        if wait:
            pass

    def calc_observable(self, iteration, traj_root, wait=True,
                        overwrite=False):

        trj_loc = self.folder_naming.trajectory_folder(iteration, traj_root)

        work_dir = self.folder_naming.property_data_folder(iteration)

        calc = OmmJobGenerator(
            self.crd_file, self.psf_file,
            template=self.obs_template, work_dir=work_dir
        )
        # constructing options dictionary
        options = dict()
        # special treatment to the option file name
        options["option_file"] = "obscalc.inp"
        options["nlip"] = self.num_lipids
        options["trj_location"] = trj_loc
        options["first_trj"] = self.first_trj
        options["last_trj"] = self.last_trj
        options["block_size"] = self.prop_block_size
        options["name"] = self.lipid_name
        # get the job run
        job = calc(2, options, overwrite=overwrite)
        job('python submit.py')
        if wait:
            pass
            #   TODO: finish this part
            # fr = FutureResult()

    def gen_reweight_target(self, iteration):

        if 'peak' in self.name or 'foot' in self.name:
            name_to_reweight = self.name.split('_')[0]
        elif 'scd' in self.name:
            name_to_reweight = self.name
        else:
            name_to_reweight = self.folder_naming.reweighting_file_name()

        self.reweight_target = ReweightTarget(
            name=name_to_reweight,
            temperature=self.temperature,
            property_dir=self.folder_naming.property_data_folder(iteration),
            energy_dir=self.folder_naming.potential_data_folder(iteration),
            property_file_template=self.property_file_format,
            result_dir=self.folder_naming.reweighting_folder(iteration),
            exp=self._exp,
            exp_dir=self.exp_dir,
            lipid=self.lipid,
            parse_groups=self.groups
        )

    def reweight(
            self, use_cluster=True, partition=None,
            force_to=False, save_result=True,
            quit=False, **kwargs
    ):
        """
        Args:
            use_cluster:
            force_to:
            save_result:
            quit: if we want to quit if there is no available result

        Returns: None
        """
        if use_cluster:
            assert partition is not None
        if (not self.reweight_target.done_reweighting and not quit) or force_to:
            self.reweight_target.reweight(
                self.first_trj, self.last_trj,
                self.pot_block_size, self.prop_block_size,
                use_cluster=use_cluster, partition=partition,
                **kwargs
            )
            if save_result:
                self.reweight_target.save_reweighted()
        elif not self.reweight_target.done_reweighting and quit:
            raise Exception(
                "Quitting ... run sensitivity calculation for property {} "
                "first!".format(self.name)
            )
        else:
            pass

    def get_robustness(self, iteration, fromfile=True, **kwargs):
        """
        The robustness here is the reciprocal of the standard DEVIATION of
        the sensitivity, which is not a rigorous definition ... but an useful
        one.
        Note: block_size can be passed through kwargs.
        """
        if not hasattr(self, 'reweight_target'):
            self.gen_reweight_target(iteration)
        if 'first_trj' not in kwargs:
            first = self.first_trj
        else:
            first = kwargs['first_trj']
        if 'last_trj' not in kwargs:
            last = self.last_trj
        else:
            last = kwargs['last_trj']
        if not fromfile:
            if not os.path.isdir(
                    self.folder_naming.robustness_folder(iteration)
            ):
                os.system(
                    "mkdir -p {}".format(
                        self.folder_naming.robustness_folder(iteration)
                    )
                )
            diff = self.reweight_target.robustness_analysis(
                first, last, self.pot_block_size, self.prop_block_size, **kwargs
            )
            if 'area' in self.name or 'scd' in self.name or 'db' in self.name:
                self.robustness = np.abs(
                    np.mean(np.array(diff), axis=0) /
                    np.std(np.array(diff), axis=0)
                )
                self.diff = np.mean(np.array(diff), axis=0)
                self.average_diff = np.mean(np.abs(np.array(diff)))
                self.uncertainty = np.std(np.array(diff), axis=0)
                np.savetxt(
                    self.folder_naming.robustness_diff_file(iteration),
                    np.mean(np.array(diff), axis=0)
                )
                np.savetxt(
                    self.folder_naming.robustness_std_file(iteration),
                    np.std(np.array(diff), axis=0)
                )
            else:
                new_diff = []
                order = order_peak_foot(self.name)
                for one_of_diff in diff:
                    new_one_diff = one_of_diff[:, 1, order]
                    new_diff.append(new_one_diff)
                self.robustness = np.abs(
                    np.mean(np.array(new_diff), axis=0) /
                    np.std(np.array(new_diff), axis=0)
                )
                self.diff = np.mean(np.array(new_diff), axis=0)
                self.average_diff = np.mean(np.abs(np.array(new_diff)))
                self.uncertainty = np.std(np.array(new_diff), axis=0)
                np.savetxt(
                    self.folder_naming.robustness_diff_file(iteration),
                    np.mean(np.array(new_diff), axis=0)
                )
                np.savetxt(
                    self.folder_naming.robustness_std_file(iteration),
                    np.std(np.array(new_diff), axis=0)
                )
        else:
            diff = np.loadtxt(
                self.folder_naming.robustness_diff_file(iteration)
            )
            std = np.loadtxt(
                self.folder_naming.robustness_std_file(iteration)
            )
            self.robustness = np.abs(diff / std)  # This is our definition!
            self.average_diff = np.mean(np.abs(diff))
            self.diff = diff
            self.uncertainty = std

    def get_sensitivity(
            self, iteration, force_redo=False,
            use_cluster=True, partition=None, quit=False
    ):
        if not hasattr(self, 'reweight_target'):
            # which is not the usual case anyway
            self.gen_reweight_target(iteration)
        self.reweight(
            force_to=force_redo, use_cluster=use_cluster, partition=partition, quit=quit
        )
        self.reweight_target.add_sensitivity_evaluator()
        # get raw
        if not ('peak' in self.name or 'foot' in self.name):
            self.deviation = - self.reweight_target.diff_sim_exp
            self.rel_deviation = - self.reweight_target.rel_diff_sim_exp
            self.sensitivity = self.reweight_target.diff_rew_sim / \
                self.perturbation
            self.rel_sensitivity = self.reweight_target.rel_diff_rew_sim / \
                self.perturbation
        else:
            # testing, the name of rdf peak should be
            # "atom-pair_peak/foot_#ofpeak/foot
            order = order_peak_foot(self.name)
            self.deviation = - np.array(
                self.reweight_target.diff_sim_exp
            )[1][order]
            self.rel_deviation = - np.array(
                self.reweight_target.rel_diff_sim_exp
            )[1][order]
            self.sensitivity = \
                self.reweight_target.diff_rew_sim[:, 1, order] / \
                self.perturbation
            self.rel_sensitivity = \
                self.reweight_target.rel_diff_rew_sim[:, 1, order] / \
                self.perturbation


