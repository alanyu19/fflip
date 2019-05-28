# -*- coding: utf-8 -*-

from __future__ import division, print_function


from fflip.tail.misfunc import *
from fflip.tail.TrainingTarget import *


class objective_function_1(object):

    def __init__(self, targets, driver_path, objfunc_counter, x):

        self.targets = targets
        self.objfunc_counter = objfunc_counter
        self.driver_path = driver_path
        self.x = x

    def __call__(self):

        import copy
        e2, s2, e3, s3, e1, s1 = copy.deepcopy(self.x)

        if self.objfunc_counter==1:
            for t in self.targets:
                t.CreateTableWithExp()
                t.create_folder()
    
        for target in self.targets:
            target.CleanUpPreviousIteration()

        ############ this whole thing should be written into a user-defined function to increase flexibility ###########
        ############ possible inputs: the sigma/epsilons, the toppar dir, the line number (found elsewhere) ... ########
        if self.objfunc_counter==1: 
            os.chdir(self.driver_path + "c5_fitting")
            os.system("rm -r olds")
            os.chdir(self.driver_path + "c6_fitting")
            os.system("rm -r olds")
            os.chdir(self.driver_path + "2d_fitting")
            os.system("rm -r olds")

        for target in self.targets:
            target.dic["CH1E_sigma"] = [s1]
            target.dic["CH1E_epsilon"] = [e1]
            target.dic["CH2E_sigma"] = [s2]
            target.dic["CH2E_epsilon"] = [e2]
            target.dic["CH3E_sigma"] = [s3]
            target.dic["CH3E_epsilon"] = [e3]
        
        substringch_1 = "CH1E\t0.0\t%.5f\t%.4f\n" % (e1, s1)
        substringch_2 = "CH2E\t0.0\t%.5f\t%.4f\n" % (e2, s2)
        substringch_3 = "CH3E\t0.0\t%.5f\t%.4f\n" % (e3, s3)

        os.chdir(self.driver_path + "toppar")
        # The line number here should be found by the program
        replace(self.driver_path + "toppar/c36ua.str", 1038, substringch_1)
        replace(self.driver_path + "toppar/c36ua.str", 1039, substringch_2)
        replace(self.driver_path + "toppar/c36ua.str", 1040, substringch_3)

        # Fit the Dihedral Parameters Using the Updated LJ Parameters

        previous_counter = self.objfunc_counter - 1

        fit_dihedral(self.driver_path, substringch_2, substringch_3, previous_counter)
        fit_dihedral_2d(self.driver_path, substringch_1, previous_counter)
        
        ################################################################################################################
        for target in self.targets:
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


class objective_function_2(object):
    # To be finished
    def __init__(self, targets, driver_path, objfunc_counter, x):
        self.targets = targets
        self.objfunc_counter = objfunc_counter
        self.driver_path = driver_path
        self.x = x

    def __call__(self):

        import copy
        e2, s2, e3, s3, e1, s1 = copy.deepcopy(self.x)
        ssrtable = np.loadtxt("previous.txt")

        for target in self.targets:
            target.CleanUpPreviousIteration()
            target.OnlyEmptyDic()
            target.dic["CH1E_sigma"] = [s1]
            target.dic["CH1E_epsilon"] = [e1]
            target.dic["CH2E_sigma"] = [s2]
            target.dic["CH2E_epsilon"] = [e2]
            target.dic["CH3E_sigma"] = [s3]
            target.dic["CH3E_epsilon"] = [e3]

            target.dic["ssr_sum"] = [ssrtable[self.objfunc_counter - 1]]
            target.WriteInfoToTable(self.driver_path, self.objfunc_counter)

        return ssrtable[self.objfunc_counter - 1]


# startpars = [-0.118, 2.192, -0.175, 2.192, -0.115, 2.133]

#####################################################################################
# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene #
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane    #
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane    #
#####################################################################################

# lower_bounds = [-0.142, 1.754, -0.210, 1.754, -0.138, 1.706]
# upper_bounds = [-0.094, 2.630, -0.140, 2.630, -0.092, 2.560]
