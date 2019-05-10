#!/usr/bin/env python

from __future__ import division, print_function

import nlopt
import copy

from fflip.tail.misfunc import *
from fflip.tail.TrainingTarget import *

# the targets should be user defined, like:
# from targets import *

class objective_function(object):

    def __init__(self, targets):

        self.targets = targets
        self.objfunc_counter = 0
        self.driver_path = os.path.dirname(os.path.abspath(__file__))

    def __call__(self, x, grad):

        self.objfunc_counter += 1
        print("objfunc_counter: ", self.objfunc_counter)


        ############ this whole thing should be written into a user-defined function to increase flexibility ###########
        ############ possible inputs: the sigma/epsilons, the toppar dir, the line number (found elsewhere) ... ########
        e2, s2, e3, s3, e1, s1 = copy.deepcopy(x)

        substringch_1 = "CH1E\t0.0\t%.4f\t%.3f\n" % (e1, s1)
        substringch_2 = "CH2E\t0.0\t%.3f\t%.3f\n" % (e2, s2)
        substringch_3 = "CH3E\t0.0\t%.3f\t%.3f\n" % (e3, s3)

        os.chdir(self.driver_path + "toppar")
        # The line number here should be found by the program
        replace(self.driver_path + "toppar/c36ua.str", 1038, substringch_1)
        replace(self.driver_path + "toppar/c36ua.str", 1039, substringch_2)
        replace(self.driver_path + "toppar/c36ua.str", 1040, substringch_3)

        # Fit the Dihedral Parameters Using the Updated LJ Parameters

        previous_counter = self.objfunc_counter - 1

        fit_dihedral(self.driver_path, substringch_2, substringch_3, previous_counter)
        fit_dihedral_2d(self.driver_path, substringch_1, previous_counter)

        for target in self.targets:
            target.dic["CH1E_sigma"] = [s1]
            target.dic["CH1E_epsilon"] = [e1]
            target.dic["CH2E_sigma"] = [s2]
            target.dic["CH2E_epsilon"] = [e2]
            target.dic["CH3E_sigma"] = [s3]
            target.dic["CH3E_epsilon"] = [e3]
        ################################################################################################################

        for target in self.targets:
            target.CleanUpPreviousIteration()
            target.Simulate()

        for target in self.targets:
            target.GetDensitykappa(self.objfunc_counter)

        for target in self.targets:
            target.Unfold_and_Gas(self.objfunc_counter)

        for target in self.targets:
            target.GetHeatofVaporization(self.objfunc_counter)

        # Diffusion (continued)
        for target in self.targets:
            target.GetDiffusionConstant(self.objfunc_counter)

        for target in self.targets:
            target.AssertDiffusionIsDone(self.objfunc_counter)

        # Heat of Vaporization (continued)
        for target in self.targets:
            target.AsserthovIsDone(self.objfunc_counter)

        # Start to Gather the Simulated/Calculated Results
        ssr_sum = 0
        for target in self.targets:
            # Is there any counter here?
            ssr = target.GetSSR(self.objfunc_counter)
            ssr_sum += ssr

        for target in self.targets:
            target.dic["ssr_sum"] = [ssr_sum]
            target.WriteInfoToTable(self.driver_path + "table/", self.objfunc_counter)

        return ssr_sum


# global variables
class optimize(object):
    """
    NLOPT optimizer to find the best parameter set
    """
    def __init__(self, targets, obj_func, startpars, lower_bounds, upper_bounds, algorithm = nlopt.LN_SBPLX):
        self.targets = targets
        self.obj_func = obj_func
        self.startpars = startpars
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.algorithm = algorithm
        assert len(list(self.lower_bounds)) == len(list(self.upper_bounds)) == len(list(startpars)),\
            "Optimizer error: lower bounds dimension doesn't match upper bounds!"
        self.dimension = len(list(self.lower_bounds))

    def __call__(self, **kwargs):
        for t in self.targets:
            t.CreateTableWithExp()
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
