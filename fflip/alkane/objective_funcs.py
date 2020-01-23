# -*- coding: utf-8 -*-

from __future__ import division, print_function


from fflip.alkane.misfunc import *
from fflip.alkane.TrainingTarget import *


class ObjectiveFunction(object):

    def __init__(
            self, targets, driver_path, objfunc_counter, x, prm_names,
            initial_var_for_table, lb=None, ub=None, fit_dihedral=True,
            use_soft_constraint=False, **kwargs
    ):
        for target in targets:
            assert isinstance(target, TrainingTarget)
        self.targets = targets
        self.objfunc_counter = objfunc_counter
        self.driver_path = driver_path
        assert len(list(x)) == len(list(prm_names))
        self.x = x
        self.prm_names = prm_names
        self.initial_var_for_table = initial_var_for_table
        self.fit_dihedral = fit_dihedral
        self.options = kwargs
        self.use_soft_constraint = use_soft_constraint
        self.lb = lb
        self.ub = ub

    @property
    def dimension(self):
        return len(list(self.x))

    def __call__(self):

        import copy

        x_copy = copy.deepcopy(self.x)

        if self.objfunc_counter == 1:
            for t in self.targets:
                t.create_table_with_exp(self.initial_var_for_table)
                t.create_folder()
    
        for target in self.targets:
            target.clean_up_previous_iteration()

        for target in self.targets:
            for name, value in zip(self.prm_names, x_copy):
                target.dic[name] = [value]

        if self.use_soft_constraint:
            out_of_boundary = apply_constraint(x_copy, self.lb, self.ub)
            if out_of_boundary:
                for target in self.targets:
                    target.dic["ssr_sum"] = [10]
                    target.write_info_to_table(
                        self.driver_path + "table/", self.objfunc_counter
                    )
                return 10
                # return a big enough figure as lost function to
                # 'reject' the current trial
        
        previous_counter = self.objfunc_counter - 1
        # a customizable fitting function, currently put in util
        if self.fit_dihedral:
            # THE fit_dihedral function should be rewritten, the first
            # argument should be the specific function for the system.
            fit_dihedrals(
               self.driver_path, self.dimension, self.options['replacelines'],
               self.prm_names, x_copy, previous_counter
            )
        else:
            change_parm_file(
                self.driver_path, self.options['replacelines'],
                self.prm_names, x_copy
            )
        for target in self.targets:
            target.simulate(self.objfunc_counter)

        for target in self.targets:
            target.get_density_kappa(self.objfunc_counter)

        for target in self.targets:
            target.unfold_and_gas(self.objfunc_counter)

        for target in self.targets:
            assert isinstance(target, TrainingTarget)
            target.get_heat_of_vaporization(self.objfunc_counter)

        # Diffusion (continued)
        for target in self.targets:
            target.get_diffusion_constant(self.objfunc_counter)

        for target in self.targets:
            target.assert_diffusion_is_done(self.objfunc_counter)

        # Heat of Vaporization (continued)
        for target in self.targets:
            target.assert_hov_is_done(self.objfunc_counter)

        # Start to Gather the Simulated/Calculated Results
        ssr_sum = 0
        for target in self.targets:
            # Is there any counter here?
            ssr = target.get_ssr(self.objfunc_counter)
            ssr_sum += ssr

        for target in self.targets:
            target.dic["ssr_sum"] = [ssr_sum]
            target.write_info_to_table(
                self.driver_path + "table/", self.objfunc_counter
            )

        return ssr_sum


class ObjectiveFunctionContinue(object):
    def __init__(
            self, targets, driver_path, objfunc_counter, x, prm_names,
            initial_var_for_table, **kwargs
    ):
        self.targets = targets
        self.objfunc_counter = objfunc_counter
        self.driver_path = driver_path
        assert len(list(x)) == len(list(prm_names))
        self.x = x
        self.prm_names = prm_names
        self.initial_var_for_table = initial_var_for_table
        self.options = kwargs

    @property
    def dimension(self):
        return len(list(self.x))

    def __call__(self):

        import copy
        x_copy = copy.deepcopy(self.x)

        os.chdir(self.driver_path)
        ssr_table = np.loadtxt("previous.txt")

        if self.objfunc_counter == 1:
            for t in self.targets:
                t.create_table_with_exp(self.initial_var_for_table)
                t.create_folder()
        
        for target in self.targets:
            target.clean_up_previous_iteration()
            target.only_empty_dic(self.initial_var_for_table)
            for name, value in zip(self.prm_names, x_copy):
                target.dic[name] = [value]

            target.dic["ssr_sum"] = [ssr_table[self.objfunc_counter - 1]]
            target.write_info_to_table(
                self.driver_path + "table/", self.objfunc_counter
            )

        return ssr_table[self.objfunc_counter - 1]
