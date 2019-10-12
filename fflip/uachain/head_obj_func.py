# -*- coding: utf-8 -*-

from fflip.uachain.misfunc import *
from fflip.uachain.TrainingTarget import *

class head_obj_func(object):
    def __init__(self, targets, driver_path, prm_names, lb=None, ub=None, usesoftconstraint=False, broken_point = 0, **kwargs):
        self.objfunc_counter = 0
        self.targets = targets
        self.driver_path = driver_path
        self.broken_point = broken_point
        self.lb = lb
        self.up = ub
        self.prm_names = prm_names
        self.usesoftconstraint = usesoftconstraint
        self.options = kwargs

    def __call__(self, x, grad):
        from fflip.uachain.objective_funcs import objective_function, objective_function_continue
        self.objfunc_counter += 1
        #print("objfunc_counter: ", self.objfunc_counter)
        if self.objfunc_counter <= int(self.broken_point):
            obj_func = objective_function_continue(
                self.targets, self.driver_path, self.objfunc_counter, x, self.prm_names, **self.options
            )
        else:
            obj_func = objective_function(
            self.targets, self.driver_path, self.objfunc_counter, x, self.prm_names,
            self.lb, self.up, self.usesoftconstraint, **self.options
            )
        return obj_func()


class head_obj_func_scipy(object):
    def __init__(self, targets, driver_path, prm_names, lb=None, ub=None, usesoftconstraint=False, broken_point = 0, **kwargs):
        self.objfunc_counter = 0
        self.targets = targets
        self.driver_path = driver_path
        self.broken_point = broken_point
        self.lb = lb
        self.up = ub
        self.prm_names = prm_names
        self.usesoftconstraint = usesoftconstraint
        self.options = kwargs

    def __call__(self, x):
        from fflip.uachain.objective_funcs import objective_function, objective_function_continue
        self.objfunc_counter += 1
        #print("objfunc_counter: ", self.objfunc_counter)
        if self.objfunc_counter <= int(self.broken_point):
            obj_func = objective_function_continue(
                self.targets, self.driver_path, self.objfunc_counter, x, self.prm_names, **self.options
            )
        else:
            obj_func = objective_function(
            self.targets, self.driver_path, self.objfunc_counter, x, self.prm_names,
            self.lb, self.up, self.usesoftconstraint, **self.options
            )
        return obj_func()

# startpars = [-0.118, 2.192, -0.175, 2.192, -0.115, 2.133]

#####################################################################################
# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene #
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane    #
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane    #
#####################################################################################

# lower_bounds = [-0.142, 1.754, -0.210, 1.754, -0.138, 1.706]
# upper_bounds = [-0.094, 2.630, -0.140, 2.630, -0.092, 2.560]
