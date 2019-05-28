# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import nlopt
import copy

from fflip.tail.misfunc import *
from fflip.tail.TrainingTarget import *

# the targets should be user defined, like:
# from targets import *

class optimize(object):
    """
    NLOPT optimizer to find the best parameter set
    """
    def __init__(self, obj_func, startpars, lower_bounds, upper_bounds, algorithm = nlopt.LN_SBPLX):
        # moved to objfunc
        self.obj_func = obj_func
        self.startpars = startpars
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.algorithm = algorithm
        assert len(list(self.lower_bounds)) == len(list(self.upper_bounds)) == len(list(startpars)),\
            "Optimizer error: lower bounds dimension doesn't match upper bounds!"
        self.dimension = len(list(self.lower_bounds))

    def __call__(self, **kwargs):
        opt = nlopt.opt(self.algorithm, self.dimension)
        opt.set_lower_bounds(self.lower_bounds)
        opt.set_upper_bounds(self.upper_bounds)
        opt.set_min_objective(self.obj_func)
        if "xtol_rel" in kwargs:
            opt.set_xtol_rel(kwargs["xtol_rel"])
        else:
            opt.set_xtol_rel(0.0001)
        x = opt.optimize(self.startpars)
        minf = opt.last_optimum_value()
        # print(x, minf)
        return x, minf


# startpars = [-0.118, 2.192, -0.175, 2.192, -0.115, 2.133]

#####################################################################################
# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene #
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane    #
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane    #
#####################################################################################

# lower_bounds = [-0.142, 1.754, -0.210, 1.754, -0.138, 1.706]
# upper_bounds = [-0.094, 2.630, -0.140, 2.630, -0.092, 2.560]
