# -*- coding: utf-8 -*-

import scipy.optimize as sopt
from fflip.omm.torsionfuncs import *
# from fflip.ljpme.util import * ## change to specific class/function(s)


def generate_weights(
        energy_series, energy_series_cross=None, criterion="boltzmann",
        cutoff=None, temperature=303.15, extra_weights=None):
    energy_series = np.array(energy_series)
    energy_series = energy_series - energy_series.min()
    energy_series_cross = np.array(energy_series_cross)
    energy_series_cross = energy_series_cross - energy_series_cross.min()
    if criterion == 'cross':
        assert energy_series_cross is not None
        assert energy_series_cross.shape == energy_series.shape
        if cutoff is not None:
            assert cutoff > 0
            new_series = np.minimum(energy_series, cutoff * np.ones(energy_series.shape[0]))
            new_series_cross = np.minimum(energy_series_cross, cutoff * np.ones(energy_series_cross.shape[0]))
            w = np.exp(-1.0 * new_series / (0.001987 * temperature)) + \
                np.exp(-1.0 * new_series_cross / (0.001987 * temperature))
        else:
            w = np.exp(-1.0 * energy_series / (0.001987 * temperature)) + \
                np.exp(-1.0 * energy_series_cross / (0.001987 * temperature))
    elif criterion == 'boltzmann':
        if cutoff is not None:
            assert cutoff > 0
            new_series = np.minimum(energy_series, cutoff * np.ones(energy_series.shape[0]))
            w = np.exp(-1.0 * new_series / (0.001987 * temperature))
        else:
            w = np.exp(-1.0 * energy_series / (0.001987 * temperature))
    elif criterion == 'uniform':
        w = np.ones(energy_series.shape[0])
    else:
        raise Exception("Criterion '{}' not accepted!".format(criterion))
    if extra_weights is not None:
        extra_weights = np.array(extra_weights)
        assert w.shape == extra_weights.shape
        return w * extra_weights
    else:
        return w


def read_dihedral_series(file):
    with open(file, 'r') as f:
        header = f.readline()
        atoms = header.split()
    dihedrals = np.loadtxt(file, skiprows=1)
    return atoms, dihedrals


def mm_energy(dihedrals, ks, ms):
    energy_series = np.zeros(dihedrals.shape[0])
    for k, m in zip(ks, ms):
        energy_series += k * (1 + np.cos(m * np.pi * dihedrals / 180))
    return energy_series


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


def rmsd_qm_mm(energy_series_1, energy_series_2, weights, offset_method):
    energy_series_1 = np.array(energy_series_1)
    energy_series_1 = energy_series_1 - energy_series_1.min()
    energy_series_2 = np.array(energy_series_2)
    energy_series_2 = energy_series_2 - energy_series_2.min()
    if offset_method == 'global_min':
        offset = 0
    elif offset_method == 'weight_guided':
        offset = (np.sum(weights * energy_series_2) - np.sum(weights * energy_series_1)) / np.sum(weights)
    sd = np.sum(weights * (energy_series_1 - energy_series_2 + offset)**2)
    msd = sd / weights.sum()
    rmsd = np.sqrt(msd)
    return rmsd


class ObjfuncDihFit(object):
    def __init__(self, dihedral_dict, m_dict, dihedral_names_sorted, qme, mme, weights, offset_method):
        """
        Args:
            dihedral_list: dictionary of dihedral series
            m_list: dictionary of multiplicities allowed by each fitted dihedral
        """
        self.dihedral_dict = dihedral_dict
        self.m_dict = m_dict
        self.dihedral_names_sorted = dihedral_names_sorted
        self.mcount_dict = dict()
        self.qme = qme
        self.mme = mme
        self.weights = weights
        self.offset_method = offset_method
        for dn in self.dihedral_names_sorted:
            self.mcount_dict[dn] = len(self.m_dict[dn])
        dn0 = self.dihedral_names_sorted[0]
        for dn in self.dihedral_names_sorted:
            assert self.dihedral_dict[dn].shape == self.dihedral_dict[dn0].shape
        assert self.qme.shape == self.dihedral_dict[dn0].shape
        assert self.mme.shape == self.dihedral_dict[dn0].shape

    def __call__(self, x):
        dn0 = self.dihedral_names_sorted[0]
        mme = copy.deepcopy(self.mme)
        k_dict = separate_k(x, self.mcount_dict, self.dihedral_names_sorted)
        for dn in self.dihedral_names_sorted:
            energy = mm_energy(self.dihedral_dict[dn], k_dict[dn], self.m_dict[dn])
            mme += energy
        if self.weights is None:  # indicating cross is used
            # TODO: wrapper for this
            weights = generate_weights(
                self.qme, mme, criterion="cross", cutoff=None, temperature=303.15, extra_weights=None
            )
        else:
            weights = self.weights
        return rmsd_qm_mm(self.qme, mme, weights, self.offset_method)


# generate_weights(energy_series, criterion="boltzmann", cutoff=None, temperature=303.15, extra_weights=None):

class DihedralFitter(object):
    def __init__(self, dihedral_files, allowed_m, qme, mme, temperature, weight_criterion, offset_method,
                 energy_cutoff=None, extra_weights=None):
        """
        Args:
            dihedral_files: list of file names
            allowed_m: list of list of multiplicities
            qme: qm energy series
            mme: mm energy series
            temperature: temperature in K
            weight_criterion: 'boltzmann' or 'uniform'
            energy_cutoff: upper bound used for boltzmann weighting
            extra_weights: extra weights to be multiplied to the normal weights
        """
        self.dihedral_dict = dict()
        self.m_dict = dict()
        self.dimension = 0
        dihedral_names = []
        for file, ms in zip(dihedral_files, allowed_m):
            atoms, dihedral_series = read_dihedral_series(file)
            dihedral_name = "{}-{}-{}-{}".format(atoms[0], atoms[1], atoms[2], atoms[3])
            self.dihedral_dict[dihedral_name] = dihedral_series
            self.m_dict[dihedral_name] = ms
            dihedral_names.append(dihedral_name)
            self.dimension += len(ms)
        dihedral_names.sort()
        self.dihedral_names_sorted = dihedral_names
        self.qme = qme
        self.mme = mme
        self.offset_method = offset_method
        self.weight_criterion = weight_criterion
        if weight_criterion != 'cross':
            self.weights = generate_weights(
                self.qme, self.mme, criterion=weight_criterion, cutoff=energy_cutoff,
                temperature=temperature, extra_weights=extra_weights
            )
        else:
            self.weights = None

    def __call__(self):
        self.obj_func = ObjfuncDihFit(
            self.dihedral_dict, self.m_dict, self.dihedral_names_sorted, self.qme, self.mme,
            self.weights, self.offset_method
        )
        # print(self.dimension)
        optimum = sopt.basinhopping(
            self.obj_func, np.zeros(self.dimension), minimizer_kwargs={"method": "BFGS"}
        )
        return optimum
