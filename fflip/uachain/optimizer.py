# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import copy
import nlopt
import scipy.optimize as sopt

from fflip.uachain.misfunc import *
from fflip.uachain.TrainingTarget import *
from fflip.uachain.head_obj_func import *

class Optimize(object):
    """The parent optimizer of all others
    """
    def __init__(self, startpars, objfunc):
        self.startpars = startpars
        self.objfunc = objfunc
    @property
    def dimension(self):
        if type(self.startpars)==int or type(self.startpars)==float or type(self.startpars)==np.float:
            return 1
        else:
            return len(list(self.startpars))
    def other_functions(self):
        pass


class NloptOptimize(Optimize):
    """
    NLOPT optimizer to search for the best parameter set
    """
    def __init__(self, objfunc, startpars, algorithm = nlopt.LN_SBPLX):
        self.algorithm = algorithm
        super().__init__(startpars, objfunc)

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
        opt.set_min_objective(self.objfunc)
        if "xtol_rel" in kwargs:
            opt.set_xtol_rel(kwargs["xtol_rel"])
        else:
            opt.set_xtol_rel(0.0001)
        x = opt.optimize(self.startpars)
        minf = opt.last_optimum_value()
        # print(x, minf)
        return x, minf

class ScipyOptimize(Optimize):
    """
    Scipy optimizer to search for the best parameter set
    """
    def __init__(self, objfunc, starpars, method, **options):
        """
        :param objfunc: the main objective funciton to use
        :param starpars: the starting parameter
        :param method: a method supported by scipy.optimize.minimize
        :param options: a dictionary that can include xtol, bounds, jac...
        """
        self.method = method
        self.options = options
        super().__init__(starpars, objfunc)

    def __call__(self):
        optimum = sopt.minimize(self.objfunc, self.startpars, method=self.method, **self.options)
        return optimum

class ScipyBrute(Optimize):
    def __init__(self):
        pass



# startpars = [-0.118, 2.192, -0.175, 2.192, -0.115, 2.133]

#####################################################################################
# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene #
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane    #
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane    #
#####################################################################################

# lower_bounds = [-0.142, 1.754, -0.210, 1.754, -0.138, 1.706]
# upper_bounds = [-0.094, 2.630, -0.140, 2.630, -0.092, 2.560]
