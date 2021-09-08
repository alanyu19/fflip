# -*- coding: utf-8 -*-

import scipy.optimize as sopt
from fflip.omm.torsionfuncs import *
from fflip.omm.util import *


def get_dihedral_parameters(atoms, psf_file, parameter_files):
    target = DihedralTarget(
        atoms, psf_file, parameter_files
    )
    target.get_cosine_series()
    return target.k, target.multp, target.phase


def compute_dihedral_energy(k, m, p, dihedrals):
    """
    Args:
        k: force constants
        m: multiplicities
        p: phases
        dihedrals: dihedral series (array like)

    Returns: energies associated with the dihedral series
    """
    dihedrals = np.array(dihedrals)
    dihedral_function = DihedralFunction(k, m, p)
    energies = dihedral_function(dihedrals)
    return energies


class ObjfuncScd(object):
    def __init__(
            self, m_dict, p_dict, ref_scds, scd_data, dihedral_data,
            original_energy, mcount_dict, dihedral_names_sorted,
            temperature, starting_param, scale
    ):
        self.m_dict = m_dict
        self.p_dict = p_dict
        self.ref_scds = ref_scds
        self.scd_data = scd_data
        self.dihedral_data = dihedral_data
        self.original_energy = original_energy
        self.mcount_dict = mcount_dict
        self.dihedral_names_sorted = dihedral_names_sorted
        self.temperature = temperature
        self.starting_param = starting_param
        self.scale = scale

    def __call__(self, x):
        x = list(x)
        k_dict = separate_k(x, self.mcount_dict, self.dihedral_names_sorted)
        total_energy = None
        for dihedral_name in self.dihedral_names_sorted:
            if isinstance(self.dihedral_data[dihedral_name][0], list):
                for dih_data in self.dihedral_data[dihedral_name]:
                    energies = compute_dihedral_energy(
                        k_dict[dihedral_name], self.m_dict[dihedral_name], self.p_dict[dihedral_name], dih_data
                    )
                    if total_energy is None:
                        total_energy = energies
                    else:
                        total_energy += energies
            else:
                energies = compute_dihedral_energy(
                    k_dict[dihedral_name], self.m_dict[dihedral_name], self.p_dict[dihedral_name],
                    self.dihedral_data[dihedral_name]
                )
                if total_energy is None:
                    total_energy = energies
                else:
                    total_energy += energies
        beta = beta_kjmol(self.temperature)
        o_energy = np.array(self.original_energy)
        p_energy = np.array(total_energy)
        tune = beta * np.mean(p_energy) - beta * np.mean(o_energy)
        ssr = 0
        for scd in self.scd_data:
            scd_data_ = self.scd_data[scd]
            obs = np.sum(scd_data_ * np.exp(
                -beta * p_energy + beta * o_energy + tune
            ))
            ptf = np.sum(np.exp(-beta * p_energy + beta * o_energy + tune))
            scd_data_rew = obs / ptf
            ssr += (scd_data_rew - self.ref_scds[scd])**2
        ssr_scd = copy.deepcopy(ssr)
        ssr += np.sum((x - np.array(self.starting_param))**2) * self.scale
        print(ssr, ssr_scd)
        return ssr


def separate_k(force_constants, mcount_dict, dihedral_names_sorted):
    k_dict = dict()
    start = 0
    end = 0
    for i, dihedral_name in enumerate(dihedral_names_sorted):
        if i == 0:
            end += mcount_dict[dihedral_name]
        else:
            start += previous_count
            end += mcount_dict[dihedral_name]
        k_dict[dihedral_name] = force_constants[start:end]
        previous_count = mcount_dict[dihedral_name]
    return k_dict


class ScdOptimizer:
    """
    An CHARMM dihedral optimizer that allows you to change the member
    (depending on the reweighter provided) of multiplicities in the original FF
    """
    def __init__(
            self, ref_scds, dihedral_dict, sim_scd_path, sim_dih_path,
            psf_file, parameter_files, scd_template, temperature, scale,
            method='BFGS', options={'eps': 1e-5, 'gtol': 1e-04},
    ):
        """
        Args:
            ref_scds: dict
            dihedral_dict:
            sim_scd_path:
            sim_dih_path:
            method:
            options:
        """
        # self.obj_func = ObjfuncDihedral(
        #     sim_dihedrals, ref_dihedrals, reweighter
        # )
        self.ref_scds = ref_scds
        self.sim_scd_path = sim_scd_path
        self.sim_dih_path = sim_dih_path
        self.dihedral_dict = dihedral_dict
        # self.reweighter = reweighter
        # self.sim_dihedrals = sim_dihedrals
        # self.ref_dihedrals = ref_dihedrals
        self.psf_file = psf_file
        self.parameter_files = parameter_files
        self.scd_template = scd_template
        self.temperature = temperature
        self.method = method
        self.dihedral_names_sorted = list(dihedral_dict.keys())
        self.dihedral_names_sorted.sort()
        self.scale = scale
        self.options = options
        self.initialized = False
        self.data_loaded = False

    def initialize_parameters(self):
        self.k_dict = dict()
        self.m_dict = dict()
        self.p_dict = dict()
        self.mcount_dict = dict()
        for dihedral_name in self.dihedral_names_sorted:
            if isinstance(self.dihedral_dict[dihedral_name][0], str):
                k_, m_, p_ = get_dihedral_parameters(
                    self.dihedral_dict[dihedral_name], self.psf_file, self.parameter_files
                )
            elif isinstance(self.dihedral_dict[dihedral_name][0], list):
                k_, m_, p_ = get_dihedral_parameters(
                    self.dihedral_dict[dihedral_name][0], self.psf_file, self.parameter_files
                )
            else:
                raise Exception('Improper dihedral dictionary!')
            self.k_dict[dihedral_name] = k_
            self.m_dict[dihedral_name] = m_
            self.p_dict[dihedral_name] = p_
            self.mcount_dict[dihedral_name] = len(k_)
        self.initial_k = []
        for dihedral_name in self.dihedral_names_sorted:
            for k in self.k_dict[dihedral_name]:
                self.initial_k.append(k)
        self.initialized = True

    def load_data(self):
        # part1, dihedral data
        self.dihedral_data = dict()
        for dihedral_name in self.dihedral_names_sorted:
            if isinstance(self.dihedral_dict[dihedral_name][0], str):
                self.dihedral_data[dihedral_name] = np.loadtxt(
                    os.path.join(self.sim_dih_path, "{}.dat".format(dihedral_name))
                )
            elif isinstance(self.dihedral_dict[dihedral_name][0], list):
                self.dihedral_data[dihedral_name] = []
                for i in range(len(self.dihedral_dict[dihedral_name])):
                    self.dihedral_data[dihedral_name].append(
                        np.loadtxt(
                            os.path.join(self.sim_dih_path, "{}_{}.dat".format(dihedral_name, i+1))
                        )
                    )
        # part2, original energies
        total_original_energy = None
        for dihedral_name in self.dihedral_names_sorted:
            if isinstance(self.dihedral_data[dihedral_name][0], list):
                for dih_data in self.dihedral_data[dihedral_name]:
                    energies = compute_dihedral_energy(
                        self.k_dict[dihedral_name], self.m_dict[dihedral_name], self.p_dict[dihedral_name],
                        dih_data
                    )
                    if total_original_energy is None:
                        total_original_energy = energies
                    else:
                        total_original_energy += energies
            else:
                energies = compute_dihedral_energy(
                    self.k_dict[dihedral_name], self.m_dict[dihedral_name], self.p_dict[dihedral_name],
                    self.dihedral_data[dihedral_name]
                )
                if total_original_energy is None:
                    total_original_energy = energies
                else:
                    total_original_energy += energies
        self.total_original_energy = total_original_energy
        # part3, scd data
        self.scd_data = dict()
        for scd in self.ref_scds:
            self.scd_data[scd] = np.loadtxt(
                os.path.join(self.sim_scd_path, self.scd_template.format(scd))
            )
        self.data_loaded = True

    def __call__(self):
        # the dimension of the optimization is bound to the reweighter,
        # which is provided to the objective function
        # start_k = np.array(self.reweighter.dihfunc2.k)
        if not self.initialized:
            self.initialize_parameters()
        if not self.data_loaded:
            self.load_data()
        self.obj_func = ObjfuncScd(
            self.m_dict, self.p_dict, self.ref_scds, self.scd_data, self.dihedral_data,
            self.total_original_energy, self.mcount_dict, self.dihedral_names_sorted,
            self.temperature, self.initial_k, self.scale
        )
        optimum = sopt.minimize(
            self.obj_func, self.initial_k, method=self.method, options=self.options
        )
        return optimum
