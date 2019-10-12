# -*- coding: utf-8 -*-

from __future__ import division, print_function


from fflip.uachain.misfunc import *
from fflip.uachain.TrainingTarget import *


class objective_function(object):

    def __init__(self, targets, driver_path, objfunc_counter, x, prm_names, lb=None, ub=None, usesoftconstraint=False, **kwargs):
        self.targets = targets
        self.objfunc_counter = objfunc_counter
        self.driver_path = driver_path
        assert len(list(x)) == len(list(prm_names))
        self.x = x
        self.prm_names = prm_names
        self.options = kwargs
        self.usesoftconstraint = usesoftconstraint
        self.lb = lb
        self.ub = ub

    @property
    def dimension(self):
        return len(list(self.x))

    def __call__(self):

        import copy

        xcopy = copy.deepcopy(self.x)

        if self.objfunc_counter==1:
            for t in self.targets:
                assert self.dimension==4 or self.dimension==6
                t.CreateTableWithExp(self.dimension)
                t.create_folder()
    
        for target in self.targets:
            target.CleanUpPreviousIteration()

        for target in self.targets:
            for name, value in zip(self.prm_names, xcopy):
                target.dic[name] = [value]

        if self.usesoftconstraint:
            out_of_boundary = applyconstraint(xcopy, self.lb, self.ub)
            if out_of_boundary:
                for target in self.targets:
                    target.dic["ssr_sum"] = [10]
                    target.WriteInfoToTable(self.driver_path + "table/", self.objfunc_counter)
                return 10   #return a big enough figure as lost function to 'reject' the current trial
        
        previous_counter = self.objfunc_counter - 1
        # a customizable fitting function, currently put in util
        fit_dihedrals(self.driver_path, self.dimension, self.options['replacelines'], self.prm_names, xcopy, previous_counter)

        for target in self.targets:
            target.Simulate(self.objfunc_counter)

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


class objective_function_continue(object):
    def __init__(self,  targets, driver_path, objfunc_counter, x, prm_names, **kwargs):
        self.targets = targets
        self.objfunc_counter = objfunc_counter
        self.driver_path = driver_path
        assert len(list(x)) == len(list(prm_names))
        self.x = x
        self.prm_names = prm_names
        self.options = kwargs

    @property
    def dimension(self):
        return len(list(self.x))

    def __call__(self):

        import copy
        xcopy = copy.deepcopy(self.x)

        os.chdir(self.driver_path)
        ssrtable = np.loadtxt("previous.txt")

        if self.objfunc_counter==1:
            for t in self.targets:
                assert self.dimension==4 or self.dimension==6
                t.CreateTableWithExp(self.dimension)
                t.create_folder()
        
        for target in self.targets:
            target.CleanUpPreviousIteration()
            target.OnlyEmptyDic(self.dimension)
            for name, value in zip(self.prm_names, xcopy):
                target.dic[name] = [value]

            target.dic["ssr_sum"] = [ssrtable[self.objfunc_counter - 1]]
            target.WriteInfoToTable(self.driver_path + "table/", self.objfunc_counter)

        return ssrtable[self.objfunc_counter - 1]


# startpars = [-0.118, 2.192, -0.175, 2.192, -0.115, 2.133]

#####################################################################################
# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene #
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane    #
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane    #
#####################################################################################

# lower_bounds = [-0.142, 1.754, -0.210, 1.754, -0.138, 1.706]
# upper_bounds = [-0.094, 2.630, -0.140, 2.630, -0.092, 2.560]
