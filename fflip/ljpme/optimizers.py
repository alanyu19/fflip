# -*- coding: utf-8 -*-

import warnings
from fflip.omm.util import filter_solution
import numpy as np
import pandas as pd
from fflip.ljpme.util import get_sign, rename_row_col


class Optimizer(object):
    def __init__(self, target_properties, special_properties,
                 model_compounds, uncertainty_scaling, parameters=None, qmc=None):
        """
        The general optimizer that takes information from targets and handles the experiments
        Args:
            target_properties (type: list): TargetProperty objects
            special_properties (type: list): SpecialProperty objects
            model_compounds（type: dict): model compounds, names as keys
            qmc (type: dict): QM partial charges for available atoms
            uncertainty_scaling (float): controls the importance of error 
        """
        self.uncertainty_scaling = uncertainty_scaling
        self.target_properties = target_properties
        self.special_properties = special_properties
        self.model_compounds = model_compounds
        if parameters is not None:
            self.parameters = parameters
        else:
            try:
                self.parameters = self.target_properties[0].parameters
            except Exception("Please provide the parameters changed!"):
                pass
        self.qmc = dict()
        if qmc is not None:
            for i, p in enumerate(self.parameters):
                if p.par_type == 'charge':
                    self.qmc[p.center_names[0]] = qmc[p.center_names[0]]
        self.W = None
        self.S = None
        self.T = None

    @property
    def model_compound_names(self):
        if self.model_compounds:
            return sorted(list(self.model_compounds.keys()))
        else:
            return list()

    @property
    def parameter_info(self):
        return self.parameters

    @property
    def all_properties(self):
        return self.target_properties + self.special_properties

    @property
    def num_qmc(self):
        return len(self.qmc)
    
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
    def num_model_compounds(self):
        return len(self.model_compounds)

    @property
    def num_parameters(self):
        return len(self.parameters)

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

        Return: a dictionary that uses the property name as the keys
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
        self.num_qmc = num_qm_charge

    def zero_qm_charge(self):
        """
        generate 0 "QM" charges, needed when weights of QM is 0
        and QM charge file is not provied
        """
        num_qm_charge = 0
        self.qm_charges = {}
        for i, p in enumerate(self.parameter_info):
            if p.par_type == 'charge':
                self.qm_charges[p.center_names[0]] = 0
                num_qm_charge += 1
        self.num_qmc = num_qm_charge

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

    def gen_weight_matrix(self, hard_bounds, drop_bounds, forbid={},
                          qmc_weight=None, qmscan_weights=None, verbose=0):
        """
        Generate the weight matrix
        Args:
            hard_bounds (type: dictionary): the lower bound of the restraint's
            force constant. If the parameter's uncertainty is less than the
            hard bound provided for its parameter type (sigma, charge, etc.),
            then the hard bound will be used as the force constant. Otherwise
            the uncertainty will be used as the force constant.
            An example could be {'sigma': 1, 'epsilon': 1, ...}.

            drop_bounds (type: dictionary): the upper bound for the
            weighted robustness. If exceeded, that parameter won't change.

            forbid （type: dictionary): parameters we don't want to touch at all.
            
            qmc_weight (type: float): the weight for QM charges.
            
            qmscan_weights (type: dictionary): the weights for model
            compounds' conformational energies
        """
        self.drop_bounds = drop_bounds
        self.forbid = forbid
        self.hard_bounds = hard_bounds
        rows = self.num_all_properties + self.num_qmc + self.num_parameters +\
               self.num_model_compounds
        cols = rows
        matrix = np.zeros((rows, cols))
        # Part1. properties
        for i, p in enumerate(self.target_properties + self.special_properties):
            matrix[i][i] = p.weight_factor
        # Part2. QM partial charges,
        i = self.num_all_properties - 1
        for p in self.parameters:
            if p.par_type == 'charge':
                i += 1
                matrix[i][i] = qmc_weight
        # Part3. deviation from original parameter set
        for i0 in range(self.num_parameters):
            ptype = self.parameters[i0].par_type
            i = i0 + self.num_all_properties + self.num_qmc
            if ptype not in forbid:
                if self.robustness[i0] >= drop_bounds[ptype]:
                    if verbose:
                        print(
                            "{0:>5s} {1:>8s} {2:>10f}".format(
                                self.parameters[i0].center_names[0], ptype,
                                round(
                                    max(self.uncertainty[i0], hard_bounds[ptype]), 6
                                )
                            )
                        )
                    matrix[i][i] = max(self.uncertainty[i0], hard_bounds[ptype])
                else:
                    if verbose:
                        print(
                            "{0:>5s} {1:>8s} {2:<30s}".format(
                                self.parameters[i0].center_names[0], ptype,
                                "  ** exceeds drop bound **"
                            )
                        )
                    matrix[i][i] = 1e5
            elif ptype in forbid:
                if self.robustness[i0] < drop_bounds[ptype]:
                    allow_change = False
                    exceed_bound = True
                else:
                    exceed_bound = False
                    for atom_name in self.parameters[i0].center_names:
                        if atom_name in forbid[ptype]:
                            allow_change = False
                        else:
                            allow_change = True
                if allow_change:
                    if verbose:
                        print(
                            "{0:>5s} {1:>8s} {2:>10f}".format(
                                self.parameters[i0].center_names[0], ptype,
                                round(
                                    max(self.uncertainty[i0], hard_bounds[ptype]), 6
                                )
                            )
                        )
                    matrix[i][i] = max(self.uncertainty[i0], hard_bounds[ptype])
                elif exceed_bound:
                    if verbose:
                        print(
                            "{0:>5s} {1:>8s} {2:<30s}".format(
                                self.parameters[i0].center_names[0], ptype,
                                "  ** exceeds drop bound **"
                            )
                        )
                    matrix[i][i] = 1e5
                else:
                    if verbose:
                        print(
                            "{0:>5s} {1:>8s} {2:<30s}".format(
                                self.parameters[i0].center_names[0], ptype,
                                "  ** not allowed for changing **"
                            )
                        )
                    matrix[i][i] = 1e5
            else:
                warnings.warn("Unexpected Situation!")
                matrix[i][i] = 1e5
        # Part4. QM scans of model compounds.
        for i0 in range(self.num_model_compounds):
            i = i0 + self.num_all_properties + self.num_qmc + \
                self.num_parameters
            matrix[i][i] = qmscan_weights[self.model_compound_names[i0]]
        self.W = matrix

    def gen_sensitivity_matrix(self, scale=1):
        """
        Args:
             scale: the scaling factor for the sensitivity. If using the
             potential energy template provided by fflip, the sensitivity
             is calculated based on (0.01 * perturbation) perturbation size
             in the "fflip obspot" step.
        """
        rows = self.num_all_properties + self.num_qmc + self.num_parameters +\
               self.num_model_compounds
        cols = self.num_parameters + self.num_model_compounds
        matrix = np.zeros((rows, cols))
        # Part1. properties
        for i, prop in enumerate(
                self.target_properties + self.special_properties
        ):
            matrix[i] = list(scale * prop.sensitivity) + \
                [0 for _ in range(self.num_model_compounds)]
        # Part2. QM (not diagonal)
        if self.num_qmc > 0:
            # initialize
            qindex = {};
            count = 0;
            qparam = []  # used to record the order
            for p in self.parameters:
                if p.par_type == 'charge':
                    qparam.append(p)
                    qindex[p.center_names[0]] = count
                # increase count no matter what par type!
                count += 1
            i = self.num_all_properties - 1
            # loop over all parameters to add up the sensitivity
            for j, p in enumerate(self.parameters):
                if p.par_type == 'charge':
                    i += 1
                    # 1.self
                    matrix[i][j] += 1 / len(p.center_names)
                    # 2. neighbors
                    for neib in p.neighbors:
                        if neib in qindex:
                            matrix[i][qindex[neib]] += 1 / len(p.neighbors)
        # Part3. C36
        for i, j in zip(
            range(
                self.num_all_properties + self.num_qmc,
                self.num_all_properties + self.num_qmc + self.num_parameters
            ),
            range(self.num_parameters)
        ):
            if self.parameters[j].par_type == 'charge':
                matrix[i][j] = 1
            else:
                if hasattr(self, 'last_solution'):
                    matrix[i][j] = 1 + 0.01 * self.last_solution[j]
                else:
                    matrix[i][j] = 1
        # Part4. model compounds' QM
        for i, j in zip(
            range(
                self.num_all_properties + self.num_qmc + self.num_parameters,
                rows
            ),
            range(self.num_parameters, self.num_parameters + self.num_model_compounds)
        ):
            matrix[i][j] = 1
        self.S = matrix

    def gen_target_vector(self):
        vector = []
        # Part1. properties
        for p in (self.target_properties + self.special_properties):
            vector.append(p.deviation)
        # Part2. QM
        if len(self.qmc) > 0:
            qindex = {}
            count = 0
            qparam = []  # used to record the order
            for p in self.parameters:
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
            for j, p in enumerate(self.parameters):
                if p.par_type == 'charge':
                    i += 1
                    # 1.self
                    qmc_vector[i] += self.qmc[p.center_names[0]] - \
                                     p.original_p
                    if hasattr(self, 'last_solution'):
                        # self continued
                        qmc_vector[i] += -1 / len(p.center_names) * \
                                         self.last_solution[j]
                        # 2. neighbors
                        for neib in p.neighbors:
                            # only consider those selected as centers
                            if neib in qindex:
                                qmc_vector[qindex[neib]] += \
                                    1 / len(p.neighbors) * \
                                    self.last_solution[j]
            for qv in qmc_vector:
                vector.append(qv)
        # Part3. Original parameters
        for i in range(self.num_parameters):
            if hasattr(self, 'last_solution'):
                vector.append(0)
                # vector.appned(-self.last_solution[i]/2)
                # vector.append(-self.last_solution[i])
            else:
                vector.append(0)
        self.T = np.array(vector)

    def load_last_solution(self, file):
        self.last_solution = filter_solution(file)


class PropertyLinearEstimator(Optimizer):
    def __init__(self, target_properties, special_properties,
                 uncertainty_scaling=300):
        """
        Test linear estimator for properties based on sensitivity info
        Args:
            targets_to_reweight: a list of reweight_target objects
            startpars: the starting vector (array like)
            algorithm: optimization algorithm
            uncertainty_scaling: the relative weight of parameter restraint
        """
        super().__init__(
            target_properties, special_properties, [], uncertainty_scaling
        )
        self.targets = []  # might not be useful in the end

    @property
    def target_prop_perturbations(self):
        """
        Returns:
            a list contains the amounts of perturbation (for rescaling
            sensitivity)
        """
        tpp = []
        for prop in self.target_properties:
            tpp.append(prop.perturbation)
        return tpp

    def update_weight(
        self,
        hard_bounds={'sigma': 0.05, 'epsilon': 0.05, 'thole': 0.05,
                     'alpha': 0.05, 'charge': 0.02},
        min_change=0.001, use_last_solution=False, soft_upper_bound=5, factor=2
    ):
        """
        min_change is the minimum change of parameters required to make the change
        really happen. This value will be compared with the parameter changes in
        the trial optimization. If the trial optimization gives change less than
        this value for a parameter, that parameter won't be changed.
        """
        self.hard_bounds = hard_bounds
        self.min_change = min_change
        # Only Part3. deviation from original parameter set
        count = 0
        for i0 in range(self.num_parameters):
            i = i0 + self.num_all_properties + self.num_qmc
            if hasattr(self, 'solution'):
                if np.abs(self.solution[i0]) < min_change:
                    # don't have too much sense to change this parameter
                    self.W[i][i] = 1e5
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

    def __call__(self, save_result=False, result_file='result.csv', 
                 ssr_file='ssr.png'):
        a = np.matmul(self.W, self.S)
        b = np.matmul(self.W, self.T)
        solution = np.linalg.lstsq(a, b, rcond=None)
        if save_result:
            residual = np.matmul(self.S, solution[0])
            residual_dict = dict()
            residual_dict['Relative error before optimization'] = \
                - self.T[:self.num_all_properties]
            residual_dict['Changed'] = residual[:self.num_all_properties]
            residual_dict['After'] = - self.T[:self.num_all_properties] + \
                residual[:self.num_all_properties]
            residual_dict['Remained'] = \
                1 - residual[:self.num_all_properties] \
                / self.T[:self.num_all_properties]
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
                                        self.num_all_properties + self.num_qmc]
            from_c36 = weighted_residual[self.num_all_properties+self.num_qmc:]
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
                    np.sqrt(((
                        residual_dict[
                            'Relative error before optimization'
                        ]*prop.scaling
                    )**2)[i])
                )
                err_after.append(
                    np.sqrt(((residual_dict['After']*prop.scaling)**2)[i])
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

