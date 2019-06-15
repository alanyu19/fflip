# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import nlopt
import copy

from fflip.tail.misfunc import *
from fflip.tail.TrainingTarget import *
from fflip.tail.main_obj_func import *

# the targets should be user defined, like:
# from targets import *

class optimize(object):
    """
    NLOPT optimizer to find the best parameter set
    """
    def __init__(self, obj_func, startpars, algorithm = nlopt.LN_SBPLX):
        # moved to objfunc
        self.obj_func = obj_func
        self.startpars = startpars
        self.algorithm = algorithm
        self.dimension = len(list(self.startpars))

    def __call__(self, **kwargs):
        opt = nlopt.opt(self.algorithm, self.dimension)
        if "lower_bounds" in kwargs:
            opt.set_lower_bounds(kwargs["lower_bounds"])
            assert len(list(kwargs["lower_bounds"])) == len(list(self.startpars)), \
                "Optimizer error: lower bounds dimension doesn't match parameter set size!"
        if "upper_bounds" in kwargs:
            opt.set_upper_bounds(kwargs["upper_bounds"])
            assert len(list(kwargs["upper_bounds"])) == len(list(self.startpars)), \
                "Optimizer error: upper bounds dimension doesn't match parameter set size!"
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
