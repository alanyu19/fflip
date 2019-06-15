# -*- coding: utf-8 -*-

from fflip.tail.misfunc import *
from fflip.tail.TrainingTarget import *

class head_obj_func(object):
    def __init__(self, targets, driver_path, broken_point):
        self.objfunc_counter = 0
        self.targets = targets
        self.driver_path = driver_path
        self.broken_point = broken_point

    def __call__(self, x, grad):
        from fflip.tail.objective_funcs import objective_function_2, objective_function_1
        self.objfunc_counter += 1
        print("objfunc_counter: ", self.objfunc_counter)
        if self.objfunc_counter <= int(self.broken_point):
            obj_func = objective_function_2(self.targets, self.driver_path, self.objfunc_counter, x)
        else:
            obj_func = objective_function_1(self.targets, self.driver_path, self.objfunc_counter, x)
        # call the sub_obj_function
        return obj_func()

# startpars = [-0.118, 2.192, -0.175, 2.192, -0.115, 2.133]

#####################################################################################
# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene #
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane    #
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane    #
#####################################################################################

# lower_bounds = [-0.142, 1.754, -0.210, 1.754, -0.138, 1.706]
# upper_bounds = [-0.094, 2.630, -0.140, 2.630, -0.092, 2.560]
