# -*- coding: utf-8 -*-

import scipy.optimize as sopt
import nlopt
from fflip.ljpme.torsionfuncs import *
# from fflip.ljpme.util import * ## change to specific class/function(s)


def generate_weights(
        energy_series, energy_series_cross=None, criterion="boltzmann",
        cutoff=None, temperature=303.15, extra_weights=None
):
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


def read_dihedral_series(filename):
    with open(filename, 'r') as f:
        header = f.readline()
        atoms = header.split()
    dihedrals = np.loadtxt(filename, skiprows=1)
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
        offset = (
            np.sum(weights * energy_series_2) -
            np.sum(weights * energy_series_1)
        ) / np.sum(weights)
    sd = np.sum(weights * (energy_series_1 - energy_series_2 + offset)**2)
    msd = sd / weights.sum()
    rmsd = np.sqrt(msd)
    return rmsd


def phase_penalty(phase_allowed, k, pforce):
    """
    currently only support 0 (-1) and 180 (1), or no restriction (0)
    :param phase_allowed: list of allowed phase (-1, 0, or 1)
    :param k: list of force constants
    :return: penalty
    """
    penalty = 0
    for p, k, f in zip(phase_allowed, k, pforce):
        if p * k <= 0:
            pass
        else:
            penalty += f * p**2 * k**2
    return penalty


def correct_phase(phase_allowed, k):
    """
    currently only support 0 (-1) and 180 (1), or no restriction (0)
    :param phase_allowed: list of allowed phase (-1, 0, or 1)
    :param k: list of force constants
    :return: correct force constant
    """
    penalty = 0
    new_k = []
    for p, k in zip(phase_allowed, k):
        if p * k <= 0:
            new_k.append(k)
        else:
            new_k.append(-k)
    return new_k


class ObjfuncDihFit(object):
    def __init__(
        self, dihedral_dict, m_dict, p_dict, pforce_dict,
        dihedral_names_sorted, qme, mme, weights, extra_weights,
        offset_method
    ):
        """
        Args:
            dihedral_list: dictionary of dihedral series
            m_list: dictionary of multiplicities allowed by each fitted dihedral
        """
        self.dihedral_dict = dihedral_dict
        self.m_dict = m_dict
        self.p_dict = p_dict
        self.pforce_dict = pforce_dict
        self.dihedral_names_sorted = dihedral_names_sorted
        self.mcount_dict = dict()
        self.qme = qme
        self.mme = mme
        self.weights = weights
        self.extra_weights = extra_weights  # only used for "cross"
        self.offset_method = offset_method
        for dn in self.dihedral_names_sorted:
            self.mcount_dict[dn] = len(self.m_dict[dn])
        dn0 = self.dihedral_names_sorted[0]
        for dn in self.dihedral_names_sorted:
            assert self.dihedral_dict[dn].shape == self.dihedral_dict[dn0].shape
        assert self.qme.shape == self.dihedral_dict[dn0].shape
        assert self.mme.shape == self.dihedral_dict[dn0].shape

    def rmsd(self, x):
        dn0 = self.dihedral_names_sorted[0]
        mme = copy.deepcopy(self.mme)
        k_dict = separate_k(x, self.mcount_dict, self.dihedral_names_sorted)
        for dn in self.dihedral_names_sorted:
            energy = mm_energy(self.dihedral_dict[dn], k_dict[dn], self.m_dict[dn])
            # new_k = correct_phase(self.p_dict[dn], k_dict[dn])
            # energy = mm_energy(self.dihedral_dict[dn], new_k, self.m_dict[dn])
            mme += energy
        if self.weights is None:  # indicating cross is used
            # TODO: wrapper for this
            weights = generate_weights(
                self.qme, mme, criterion="cross", cutoff=None,
                temperature=303.15, extra_weights=self.extra_weights
            )
        else:
            weights = self.weights
        rmsd = rmsd_qm_mm(self.qme, mme, weights, self.offset_method)
        return rmsd
    
    def __call__(self, x, grad):
        dn0 = self.dihedral_names_sorted[0]
        mme = copy.deepcopy(self.mme)
        penalty = 0
        k_dict = separate_k(x, self.mcount_dict, self.dihedral_names_sorted)
        for dn in self.dihedral_names_sorted:
            penalty += phase_penalty(self.p_dict[dn], k_dict[dn], self.pforce_dict[dn])
            energy = mm_energy(self.dihedral_dict[dn], k_dict[dn], self.m_dict[dn])
            # new_k = correct_phase(self.p_dict[dn], k_dict[dn])
            # energy = mm_energy(self.dihedral_dict[dn], new_k, self.m_dict[dn])
            mme += energy
        if self.weights is None:  # indicating cross is used
            # TODO: wrapper for this
            weights = generate_weights(
                self.qme, mme, criterion="cross", cutoff=None,
                temperature=303.15, extra_weights=self.extra_weights
            )
        else:
            weights = self.weights
        rmsd = rmsd_qm_mm(self.qme, mme, weights, self.offset_method)
        total = rmsd + penalty
        # print(rmsd, penalty)
        # return rmsd
        return total


# generate_weights(energy_series, criterion="boltzmann", cutoff=None,
# temperature=303.15, extra_weights=None):

class DihedralFitter(object):
    def __init__(self, dihedral_files, allowed_m, phase, pforce, qme, mme,
                 temperature, weight_criterion, offset_method,
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
        self.p_dict = dict()
        self.pforce_dict = dict()
        self.mcount_dict = dict()
        self.dimension = 0
        dihedral_names = []
        for file, ms, ps, pf in zip(dihedral_files, allowed_m, phase, pforce,):
            atoms, dihedral_series = read_dihedral_series(file)
            dihedral_name = "{}-{}-{}-{}".format(atoms[0], atoms[1], atoms[2], atoms[3])
            self.dihedral_dict[dihedral_name] = dihedral_series
            self.m_dict[dihedral_name] = ms
            self.p_dict[dihedral_name] = ps
            self.pforce_dict[dihedral_name] = pf
            dihedral_names.append(dihedral_name)
            self.dimension += len(ms)
        dihedral_names.sort()
        self.dihedral_names_sorted = dihedral_names
        for dn in self.dihedral_names_sorted:
            self.mcount_dict[dn] = len(self.m_dict[dn])
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
            self.extra_weights = extra_weights
            self.weights = None
        self.optimum = None

    def nlopt_fit(self, method, lower_bounds=None, upper_bounds=None, maxiter=None):
        self.obj_func = ObjfuncDihFit(
            self.dihedral_dict, self.m_dict, self.p_dict, self.pforce_dict,
            self.dihedral_names_sorted, self.qme, self.mme,
            self.weights, self.extra_weights, self.offset_method
        )
        # print(self.dimension)
        # optimum = sopt.basinhopping(
        #    self.obj_func, np.zeros(self.dimension), minimizer_kwargs={"method": "BFGS"}
        # )
        opt = nlopt.opt(method, self.dimension)
        if lower_bounds is None:
            lower_bounds = - np.ones(self.dimension) * 3
        opt.set_lower_bounds(lower_bounds)
        if upper_bounds is None:
            upper_bounds = np.ones(self.dimension) * 3
        opt.set_upper_bounds(upper_bounds)
        if maxiter is not None:
            opt.set_maxeval(maxiter)
        opt.set_min_objective(self.obj_func)
        opt.set_xtol_rel(0.005)  # hard-coded
        if method is nlopt.G_MLSL_LDS or method is nlopt.G_MLSL:
            opt.set_local_optimizer(nlopt.opt(nlopt.LN_SBPLX, self.dimension))
        # x = opt.optimize(np.zeros(self.dimension))
        # print(max(lower_bounds), min(upper_bounds))
        x = opt.optimize(
            np.random.uniform(
                low=max(lower_bounds), high=min(upper_bounds), size=(self.dimension)
            )
        )
        minf = opt.last_optimum_value()
        # print(x, minf)
        self.optimum = x
        self.opt = opt

    def simulated_annealing(self, lower_bounds=None, upper_bounds=None, nsteps=10000, tempr0=1000):
        """
        Args:
            lower_bounds (array like): lower bounds for k
            upper_bounds (array like): upper bounds for k
            nsteps (int): total steps to run
            tempr0 (float): initial temperature for the simulated annealing algorithm
        """
        import logging
        import random
        logger = logging.getLogger(".")
        logger.setLevel(logging.DEBUG)
        if os.path.isfile("rmsd.log"):
            os.system("rm rmsd.log")
        file_handler = logging.FileHandler('rmsd.log')
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)
        if lower_bounds == None:
            lower_bounds = -3 * np.ones(self.dimension)
        if upper_bounds == None:
            upper_bounds = 3 * np.ones(self.dimension)
        x = np.random.uniform(
            low=max(lower_bounds), high=min(upper_bounds), size=(self.dimension)
        )
        assert np.array(lower_bounds).shape == np.array(upper_bounds).shape == x.shape
        self.obj_func = ObjfuncDihFit(
            self.dihedral_dict, self.m_dict, self.p_dict, self.pforce_dict,
            self.dihedral_names_sorted, self.qme, self.mme,
            self.weights, self.extra_weights, self.offset_method
        )
        # initialize with large rmsd
        rmsd_old = 1e2
        rmsd_best = 1e2
        step = 0
        while step < nsteps:
            step += 1
            tempr = tempr0 * np.exp(-(step/(nsteps/4.0)))
            rmsd = self.obj_func(x, 0)
            if step == 0:
                p =0.0
            else:
                drmsd = rmsd - rmsd_old
                boltz = -1.0 * drmsd / (0.001987*tempr)
                p  = np.exp(boltz)
            p0 = random.uniform(0.0, 1.0)
            if p0 < p:
                accepted = 1
            else:
                accepted = 0
            if accepted:
                rmsd_old = rmsd
                if rmsd_old < rmsd_best:
                    rmsd_best = rmsd_old
                    x_best = x
            logger.info("{}".format(round(rmsd_best, 4)))
            for i in range(len(x)):
                # the max random move size is hard-coded
                x[i] = x[i] + random.uniform(-0.5, 0.5)
                x[i] = max(lower_bounds[i], x[i])
                x[i] = min(upper_bounds[i], x[i])
        self.optimum = x_best
        self.rmsd_best = rmsd_best
    
    def quick_rmsd(self, x):
        return self.obj_func.rmsd(x)

    def show_optimum(self, correct_phase=False):
        assert self.optimum is not None
        # k_dict = separate_k(self.optimum.x, self.mcount_dict, self.dihedral_names_sorted)
        k_dict = separate_k(self.optimum, self.mcount_dict, self.dihedral_names_sorted)
        for dn in self.dihedral_names_sorted:
            if correct_phase:
                k = correct_phase(self.p_dict[dn], k_dict[dn])
            else:
                k = k_dict[dn]
            print(dn, k)
