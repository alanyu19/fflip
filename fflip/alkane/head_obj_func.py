# -*- coding: utf-8 -*-


class HeadObjFunc(object):
    def __init__(self, targets, driver_path, prm_names, initial_var_for_table,
                 lb=None, ub=None, fit_dihedral=True,
                 use_soft_constraint=False, broken_point=0, **kwargs):
        self.objfunc_counter = 0
        self.targets = targets
        self.driver_path = driver_path
        self.broken_point = broken_point
        self.initial_var_for_table = initial_var_for_table
        self.lb = lb
        self.up = ub
        self.fit_dihedral = fit_dihedral
        self.prm_names = prm_names
        self.use_soft_constraint = use_soft_constraint
        self.options = kwargs

    def __call__(self, x, grad):
        from fflip.alkane.objective_funcs import ObjectiveFunction,\
            ObjectiveFunctionContinue
        self.objfunc_counter += 1
        # print("objfunc_counter: ", self.objfunc_counter)
        if self.objfunc_counter <= int(self.broken_point):
            obj_func = ObjectiveFunctionContinue(
                self.targets, self.driver_path, self.objfunc_counter, x,
                self.prm_names, self.initial_var_for_table, **self.options
            )
        else:
            obj_func = ObjectiveFunction(
                self.targets, self.driver_path, self.objfunc_counter, x,
                self.prm_names, self.initial_var_for_table, self.lb, self.up,
                self.fit_dihedral, self.use_soft_constraint,
                **self.options
            )
        return obj_func()


class HeadObjFuncScipy(object):
    def __init__(self, targets, driver_path, prm_names, initial_var_for_table,
                 lb=None, ub=None, fit_dihedral=True,
                 use_soft_constraint=False, broken_point=0, **kwargs):
        self.objfunc_counter = 0
        self.targets = targets
        self.driver_path = driver_path
        self.broken_point = broken_point
        self.initial_var_for_table = initial_var_for_table
        self.fit_dihedral = fit_dihedral
        self.lb = lb
        self.up = ub
        self.prm_names = prm_names
        self.use_soft_constraint = use_soft_constraint
        self.options = kwargs

    def __call__(self, x):
        from fflip.alkane.objective_funcs import ObjectiveFunction, \
            ObjectiveFunctionContinue
        self.objfunc_counter += 1
        # print("objfunc_counter: ", self.objfunc_counter)
        if self.objfunc_counter <= int(self.broken_point):
            obj_func = ObjectiveFunctionContinue(
                self.targets, self.driver_path, self.objfunc_counter,
                x, self.prm_names, self.initial_var_for_table, **self.options
            )
        else:
            obj_func = ObjectiveFunction(
                self.targets, self.driver_path, self.objfunc_counter,
                x, self.prm_names, self.initial_var_for_table,
                self.lb, self.up, self.fit_dihedral, self.use_soft_constraint,
                **self.options
            )
        return obj_func()
