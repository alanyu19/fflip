# -*- coding: utf-8 -*-

from coffe.omm.util import filter_solution
from fflip.ljpme.moreutil import *
import numpy as np
import pandas as pd


class Optimizer(object):
    def __init__(self, target_properties, special_properties,
                 perturbation_baseline = 1):
        pass

    def exchange_exp_value(self, pname1, pname2):
        for p in self.all_properties:
            if p.name == pname1:
                exp1 = p._exp
            if p.name == pname2:
                exp2 = p._exp
        for p in self.all_properties:
            if p.name == pname1:
                p.exchanged_exp = exp2
            if p.name == pname2:
                p.exchanged_exp = exp1

    def gen_target_properties(self):
        self.target_properties = 1

    def load_last_solution(self, file):
        """
        Get the solution of parameter from last iteration
        :param file:
        :return: a list contains the simplified solution
        """
        self.last_solution = filter_solution(file)


class PropertyLinearEstimator(Optimizer):
    def __init__(self, target_properties, special_properties,
                 perturbation_baseline = 1, uncertainty_scaling = 300):
        """
        Test linear estimator for properties based on sensitivity info
        Args:
            targets_to_reweight: a list of reweight_target objects
            startpars: the starting vector
            algorithm: optimization algorithm
        """
        self.perturbation_baseline = perturbation_baseline
        self.target_properties = target_properties
        self.special_properties = special_properties
        self.targets = []  # might not be useful in the end
        self.uncertainty_scaling = uncertainty_scaling

    @property
    def all_properties(self):
        return self.target_properties + self.special_properties

    @property
    def target_prop_perturbations(self):
        """
        Returns:
            a list contains the percentages of perturbation (for rescaling
            sensitivity)
        """
        tpp = []
        for prop in self.target_properties:
            tpp.append(prop.perturbation)
        return tpp

    def save_sensitivity_table(self):
        for target in self.targets:
            pass

    def plot_sensitivity(self):
        pass

    @property
    def num_properties(self):
        return len(self.target_properties)

    @property
    def num_special_properties(self):
        return len(self.special_properties)

    @property
    def num_all_properties(self):
        return self.num_properties + self.num_special_properties

    @property
    def parameter_info(self):
        # In the current version of code, the parameter info is uniform for
        # all properties studied
        prop = self.target_properties[0]
        return prop.parameters

    @property
    def num_parameters(self):
        return len(self.parameter_info)

    @property
    def targets_sensitivity(self):
        """
        The difference between the reweighted and the simulated,
        it depends on the change of parameters.

        Returns: a dictionary with the property names as keys

        """
        sensitivity_dict = {}
        for target in self.targets:
            sensitivity_dict[target.name] = target.diff_rew_sim
        return sensitivity_dict

    @property
    def targets_sensitivity_rel(self):
        """
        The relative (divided experiments) difference between the reweighted
        and the simulated,
        it depends on the change of parameters.

        Return: a dictionary with the property names as keys
        """
        sensitivity_dict = {}
        for target in self.targets:
            sensitivity_dict[target.name] = target.rel_diff_rew_sim
        return sensitivity_dict

    @property
    def targets_deviation(self):
        """
        The difference between last run (sim) to exp
        :return: a dictionary that uses the property name as the keys
        """
        dev_dict = {}
        for target in self.targets:
            dev_dict[target.name] = target.diff_sim_exp
        return dev_dict

    @property
    def targets_deviation_rel(self):
        """
        The relative difference between last run (sim) to exp
        :return: a dictionary that uses the property name as the keys
        """
        dev_dict = {}
        for target in self.targets:
            dev_dict[target.name] = target.rel_diff_sim_exp
        return dev_dict

    @property
    def robustness(self):
        roblist = []
        weights = []
        for prop in self.target_properties:
            if hasattr(prop, 'robustness'):
                weights.append(prop.weight_factor)
                roblist.append(
                    prop.weight_factor * prop.robustness
                )
        mean_weight = np.mean(np.array(weights))
        return np.mean(np.array(roblist), axis=0) / mean_weight

    @property
    def uncertainty(self):
        uncerlist = []
        weights = []
        for prop in self.target_properties:
            if hasattr(prop, 'robustness'):
                # The weights here is used to scale the robustness!
                weights.append(prop.weight_factor)
                uncerlist.append(
                    prop.weight_factor * prop.uncertainty
                )
        mean_weight = np.mean(np.array(weights))
        return self.uncertainty_scaling * np.mean(np.array(uncerlist), axis=0) \
            / mean_weight

    def get_qm_charge_from_file(self, filename):
        """

        Args:
            filename: the file name stores the QM charges, should be .csv file
        Returns:
        """
        charge_info = pd.read_csv(filename)
        self.qm_info = {}
        for atom, charge in zip(
                charge_info['atom'].tolist(), charge_info['charge'].tolist()
        ):
            self.qm_info[atom] = charge
        num_qm_charge = 0
        self.qm_charges = {}
        for i, p in enumerate(self.parameter_info):
            if p.par_type == 'charge':
                self.qm_charges[p.center_names[0]] = self.qm_info[
                    p.center_names[0]]
                num_qm_charge += 1
        self.num_qm = num_qm_charge

    def get_weight_matrix(
            self, qm_weight=0.01,
            hard_bounds={'sigma': 0.05, 'epsilon': 0.05, 'charge': 0.02},
            drop_bounds={'sigma': 999, 'epsilon': 999, 'charge': 999},
            forbid={'epsilon': ['C12', 'H11A', 'P', 'C2', 'N', 'HS', 'C13'],
                    'sigma': ['C12', 'H11A', 'P', 'C2', 'N', 'HS', 'C13']}
    ):
        """
        Generate the weight matrix
        Args:
            qm_weight: float, the weight for QM charges.
            hard_bounds: dictionary, the lower bound of the robustness,
            key should be the parameter type ('charge',
            ...).
            drop_bounds: dictionary, the upper bound for the uncertainty,
            if exceeded, that parameter won't change.
            forbid: dictionary, parameters we don't want to touch at all.
        """
        rows = self.num_all_properties + self.num_qm + self.num_parameters
        cols = rows
        matrix = np.zeros((rows, cols))
        # Part1. properties
        for i, p in enumerate(self.target_properties + self.special_properties):
            matrix[i][i] = p.weight_factor
        # Part2. QM calculations for partial charges,
        # Todo: write a more generic class for handling different parts of
        #  the weight maxtrix.
        i = self.num_all_properties - 1
        for p in self.parameter_info:
            if p.par_type == 'charge':
                i += 1
                matrix[i][i] = qm_weight
        # Part3. deviation from original parameter set
        for i0 in range(self.num_parameters):
            ptype = self.parameter_info[i0].par_type
            i = i0 + self.num_all_properties + self.num_qm
            if ptype not in forbid:
                if not self.robustness[i0] < drop_bounds[ptype]:
                    print(
                        "{0:>5s} {1:>8s} {2:>10f}".format(
                            self.parameter_info[i0].center_names[0], ptype,
                            round(
                                max(self.uncertainty[i0], hard_bounds[ptype]), 6
                            )
                        )
                    )
                    matrix[i][i] = max(self.uncertainty[i0], hard_bounds[ptype])
                else:
                    print(
                        "{0:>5s} {1:>8s} {2:<30s}".format(
                            self.parameter_info[i0].center_names[0], ptype,
                            "  ** exceeds drop bound **"
                        )
                    )
                    matrix[i][i] = 999999
            elif ptype in forbid:
                if self.robustness[i0] < drop_bounds[ptype]:
                    allow_change = False
                    exceed_bound = True
                else:
                    exceed_bound = False
                    allow_change = True
                    for atom_name in self.parameter_info[i0].center_names:
                        if atom_name in forbid[ptype]:
                            allow_change = False
                        else:
                            continue
                if allow_change:
                    print(
                        "{0:>5s} {1:>8s} {2:>10f}".format(
                            self.parameter_info[i0].center_names[0], ptype,
                            round(
                                max(self.uncertainty[i0], hard_bounds[ptype]), 6
                            )
                        )
                    )
                    matrix[i][i] = max(self.uncertainty[i0], hard_bounds[ptype])
                elif exceed_bound:
                    print(
                        "{0:>5s} {1:>8s} {2:<30s}".format(
                            self.parameter_info[i0].center_names[0], ptype,
                            "  ** exceeds drop bound **"
                        )
                    )
                    matrix[i][i] = 999999
                else:
                    print(
                        "{0:>5s} {1:>8s} {2:<30s}".format(
                            self.parameter_info[i0].center_names[0], ptype,
                            "  ** not allowed for changing **"
                        )
                    )
                    matrix[i][i] = 999999
            else:
                print("???")
                matrix[i][i] = 999999
        self.W = matrix

    def get_sensitivity_matrix(self):
        rows = self.num_all_properties + self.num_qm + self.num_parameters
        cols = self.num_parameters
        matrix = np.zeros((rows, cols))
        # Part1. properties
        for i, prop in enumerate(
                self.target_properties + self.special_properties
        ):
            matrix[i] = prop.sensitivity
        # Part2. QM (not diagonal)
        # initialize
        qindex = {}; count = 0; qparam = []  # used to record the order
        for p in self.parameter_info:
            if p.par_type == 'charge':
                qparam.append(p)
                qindex[p.center_names[0]] = count
            # increase count no matter what par type!
            count += 1
        i = self.num_all_properties - 1
        # loop over all parameters to add up the sensitivity
        for j, p in enumerate(self.parameter_info):
            if p.par_type == 'charge':
                i += 1
                # 1.self
                matrix[i][j] += 0.01/len(p.center_names) * \
                                get_sign(p.original_p)
                # 2. neighbors
                for neib in p.neighbors:
                    # only consider those selected as centers?
                    if neib in qindex:
                        matrix[i][qindex[neib]] += -0.01/len(p.neighbors) * \
                                                   get_sign(p.original_p)
        # Part3. C36
        for i, j in zip(
                range(self.num_all_properties + self.num_qm, rows),
                range(self.num_parameters)
        ):
            if self.parameter_info[j].par_type == 'charge':
                matrix[i][j] = 1
            else:
                if hasattr(self, 'last_solution'):
                    matrix[i][j] = 1 + 0.01 * self.last_solution[j]
                else:
                    matrix[i][j] = 1
        self.S = matrix

    def get_deviation_vector(self):
        vector = []
        # Part1. properties
        for p in (self.target_properties + self.special_properties):
            vector.append(p.deviation)
        # Part2. QM
        qindex = {}
        count = 0
        qparam = []  # used to record the order
        for p in self.parameter_info:
            if p.par_type == 'charge':
                qparam.append(p)
                qindex[p.center_names[0]] = count
                # increase count only when par type is 'charge'
                count += 1
        # prepare an empty list
        qmc_vector = list(np.zeros(len(qparam)))
        # loop over all parameters to add up the deviation
        # index i is used to count the charge parameter
        # index j is used to retrieve the last solution
        i = -1
        for j, p in enumerate(self.parameter_info):
            if p.par_type == 'charge':
                i += 1
                # 1.self
                qmc_vector[i] += self.qm_charges[p.center_names[0]] - \
                                 p.original_p
                if hasattr(self, 'last_solution'):
                    # self continued
                    qmc_vector[i] += -0.01/len(p.center_names) * \
                                     get_sign(p.original_p) * \
                                     self.last_solution[j]
                    # 2. neighbors
                    for neib in p.neighbors:
                        # only consider those selected as centers
                        if neib in qindex:
                            qmc_vector[qindex[neib]] += \
                                0.01/len(p.neighbors) * \
                                get_sign(p.original_p) * \
                                self.last_solution[j]
        for qv in qmc_vector:
            vector.append(qv)
        # Part3. C36
        for i in range(self.num_parameters):
            if hasattr(self, 'last_solution'):
                vector.append(0)
                # vector.appned(-self.last_solution[i]/2)
                # vector.append(-self.last_solution[i])
            else:
                vector.append(0)
        self.F = np.array(vector)

    def update_weight(
        self,
        hard_bounds={'sigma': 0.05, 'epsilon': 0.05, 'charge': 0.02},
        lower_bound=0.3, use_last_solution=False, soft_upper_bound=5, factor=2
    ):
        # Only Part3. deviation from original parameter set
        count = 0
        for i0 in range(self.num_parameters):
            i = i0 + self.num_all_properties + self.num_qm
            if hasattr(self, 'solution'):
                if np.abs(self.solution[i0]) < lower_bound:
                    self.W[i][i] = 100000
                else:
                    count = count + 1
                    ptype = self.parameter_info[i0].par_type
                    self.W[i][i] = max(self.uncertainty[i0], hard_bounds[ptype])
                    if use_last_solution and hasattr(self, 'last_solution'):
                        sub = soft_upper_bound
                        if self.last_solution[i0] + self.solution[i0] > sub:
                            # currently very crude (only use plus, even for lj)
                            print(
                                "Putting more restriction on {}-{}".format(
                                    self.parameter_info[i0].center_names[0],
                                    self.parameter_info[i0].par_type,
                                )
                            )
                            self.W[i][i] = self.W[i][i] * factor
            else:
                raise Exception('Error in update_weight!')
        print('{} parameters changed in total'.format(count))

    def __call__(
            self,
            save_result=False,
            result_file='linear_solver_result.txt',
            ssr_file='ssr.png'):
        a = np.matmul(self.W, self.S)
        b = np.matmul(self.W, self.F)
        solution = np.linalg.lstsq(a, b, rcond=None)
        if save_result:
            residual = np.matmul(self.S, solution[0])
            residual_dict = dict()
            residual_dict['before'] = self.F[:self.num_all_properties]
            residual_dict['changed'] = residual[:self.num_all_properties]
            residual_dict['remained'] = self.F[:self.num_all_properties] - \
                residual[:self.num_all_properties]
            residual_dict['%remained'] = \
                1 - residual[:self.num_all_properties] \
                / self.F[:self.num_all_properties]
            to_write = pd.DataFrame(residual_dict)
            row_renaming = rename_row_col(
                [p.name for p in self.target_properties +
                 self.special_properties]
            )
            to_write.rename(index=row_renaming, inplace=True)
            to_write.to_csv(result_file)
            # ...................................
            # Recoding SSR from different sources:
            weighted_residual = np.matmul(a, solution[0]) - b
            from_prop = weighted_residual[:self.num_all_properties]
            from_qm = weighted_residual[self.num_all_properties:
                                        self.num_all_properties + self.num_qm]
            from_c36 = weighted_residual[self.num_all_properties+self.num_qm:]
            print(
                "Summary of weighted residues:",
                np.sum(from_prop**2) + np.sum(from_qm**2) + np.sum(from_c36**2)
            )
            print("from properties:", np.sum(from_prop**2))
            print("from QM:", np.sum(from_qm**2))
            print("from C36:", np.sum(from_c36**2))
            # ................................................................
            # Detailed for each target property (including special properties):
            indexes = []
            names = []
            sr = []
            err_before = []
            err_after = []
            for i, prop in enumerate(self.all_properties):
                indexes.append(i*1.2)
                names.append(prop.name)
                sr.append((from_prop**2)[i])
                err_before.append(
                    np.sqrt(((residual_dict['before']*prop.scaling)**2)[i])
                )
                err_after.append(
                    np.sqrt(((residual_dict['remained']*prop.scaling)**2)[i])
                )
            from matplotlib import pyplot as plt
            fig, ax = plt.subplots(figsize=(22, 5))
            ax.bar(
                np.array(indexes) - 0.3, sr, width=0.3,
                label='contribution to optimization residue'
            )
            ax.bar(
                np.array(indexes), err_before, width=0.3,
                label='scaled error before'
            )
            ax.bar(
                np.array(indexes) + 0.3, err_after, width=0.3,
                label='scaled error remaining'
            )
            plt.xticks(indexes, names, rotation='vertical')
            plt.legend(fontsize=14)
            plt.savefig(ssr_file, bbox_inches='tight')
        self.solution = solution[0]
        return solution[0]

