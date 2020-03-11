# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import nlopt
import scipy.optimize as sopt


class Optimize(object):
    """The parent optimizer of all others
    """
    def __init__(self, startpars, objfunc):
        self.startpars = startpars
        self.objfunc = objfunc

    @property
    def dimension(self):
        if type(self.startpars) == int or type(self.startpars) == float \
                or type(self.startpars) == np.float:
            return 1
        else:
            return len(list(self.startpars))

    def other_functions(self):
        pass


class NloptOptimize(Optimize):
    """
    NLOPT optimizer to search for the best parameter set
    """
    def __init__(self, objfunc, startpars, algorithm=nlopt.LN_SBPLX):
        self.algorithm = algorithm
        super().__init__(startpars, objfunc)

    def __call__(self, **kwargs):
        opt = nlopt.opt(self.algorithm, self.dimension)
        if "lower_bounds" in kwargs:
            opt.set_lower_bounds(kwargs["lower_bounds"])
            assert \
                len(list(kwargs["lower_bounds"])) == len(list(self.startpars)),\
                "Optimizer error: lower bounds dimension " \
                "doesn't match parameter set size!"
        if "upper_bounds" in kwargs:
            opt.set_upper_bounds(kwargs["upper_bounds"])
            assert \
                len(list(kwargs["upper_bounds"])) == len(list(self.startpars)),\
                "Optimizer error: upper bounds dimension " \
                "doesn't match parameter set size!"
        opt.set_min_objective(self.objfunc)
        if "xtol_rel" in kwargs:
            opt.set_xtol_rel(kwargs["xtol_rel"])
        else:
            opt.set_xtol_rel(0.0001)
        print(self.startpars)
        x = opt.optimize(self.startpars)
        min_f = opt.last_optimum_value()
        # print(x, min_f)
        return x, min_f


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
        optimum = sopt.minimize(
            self.objfunc, self.startpars, method=self.method, **self.options
        )
        return optimum


class ScipyBrute:
    def __init__(self, objfunc, lowers, uppers, intervals):
        self.objfunc = objfunc
        assert len(lowers) == len(uppers) and len(lowers) == len(intervals)
        rrange = ()
        for l, u, i in zip(lowers, uppers, intervals):
            rrange = rrange + (slice(l, u, i),)
        print(rrange)
        self.rrange = rrange

    def __call__(self):
        optimum = sopt.brute(self.objfunc, self.rrange, full_output=True)
        return optimum
