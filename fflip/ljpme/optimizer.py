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
        #self.lower_bounds = lower_bounds
        #self.upper_bounds = upper_bounds
        self.algorithm = algorithm
        self.dimension = len(list(self.startpars))

    def __call__(self, **kwargs):
        opt = nlopt.opt(self.algorithm, self.dimension)
        if "lower_bounds" in kwargs:
            opt.set_lower_bounds(kwargs["lower_bounds"])
            assert len(list(self.lower_bounds)) == len(list(startpars)),\
            "Optimizer error: lower bounds dimension doesn't match startpars!"
        if "upper_bounds" in kwargs:
            opt.set_upper_bounds(kwargs["upper_bounds"])
            assert len(list(self.upper_bounds)) == len(list(startpars)),\
            "Optimizer error: upper bounds dimension doesn't match startpars!"
        opt.set_min_objective(self.obj_func)
        if "xtol_rel" in kwargs:
            opt.set_xtol_rel(kwargs["xtol_rel"])
        else:
            opt.set_xtol_rel(0.02)
        if "xtol_abs" in kwargs:
            print("干他妈的")
            opt.set_xtol_abs(kwargs["xtol_abs"])
        else:
            opt.set_xtol_abs(0.001)
            print("Warning: there isn't a effective stopping criteria specified!")
        print("Starting Optimization ...")
        print("abs xtol is {}".format(opt.get_xtol_abs()))
        print("rel xtol is {}".format(opt.get_xtol_rel()))
        x = opt.optimize(self.startpars)
        minf = opt.last_optimum_value()
        # print(x, minf)
        return x, minf


