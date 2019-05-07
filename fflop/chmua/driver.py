#!/usr/bin/env python

from __future__ import division, print_function

import nlopt
import numpy as np
import MDAnalysis as mda
import pandas as pd

import os
import time
import tempfile
import shutil
import sys
import copy
from TrainingTarget import *
from misfunc import *
from targets import *

def objfunc(x, grad):
    global objfunc_counter
    global lower_bounds
    global upper_bounds
    global y
    global driver_path
    global kT
    global nframes

    objfunc_counter += 1
    print("objfunc_counter: ", objfunc_counter)

    a, b = copy.deepcopy(x)
    #d = c/1.05
    
    # Make Changes to Stream Files
    substringch1 = "CH1E\t0.0\t%.4f\t%.3f\n" % (a, b)
    os.chdir(driver_path + "toppar")        
    replace(driver_path+"toppar/c36ua.str", 1038, substringch1)  # check line number

    previous_counter = objfunc_counter - 1
    fit_dihedral_2d(driver_path, substringch1, previous_counter)
    
    for system in syslist:
        system.CleanUpPreviousIteration()
        system.Simulate()
        system.dic["CH1E_sigma"] = [a]
        system.dic["CH1E_epsilon"] = [b]
        
    for system in syslist:
        system.GetDensitykappa(objfunc_counter)
                    
    for system in syslist:
        system.Unfold_and_Gas(objfunc_counter)
        
    for system in syslist:
        system.GetHeatofVaporization(objfunc_counter)
        
    # Diffusion (continued)    
    for system in syslist:
        system.GetDiffusionConstant(objfunc_counter)

    for system in syslist:
        system.AssertDiffusionIsDone(objfunc_counter)
        
    # Heat of Vaporization (continued)
    for system in syslist:
        system.AsserthovIsDone(objfunc_counter)
              
    # Start to Gather the Simulated/Calculated Results
    ssr_sum = 0
    for system in syslist:
        # Is there any counter here?
        ssr = system.GetSSR(objfunc_counter)
        ssr_sum += ssr
    
    for system in syslist:
        system.dic["ssr_sum"] = [ssr_sum]
        system.WriteInfoToTable(driver_path, objfunc_counter)
    
    return ssr_sum


# global variables
startpars = [-0.115, 2.133]
#driver_path = "/u/alanyu/optim/trial_1/"
objfunc_counter = 0
lower_bounds = [-0.15, 1.800]
upper_bounds = [-0.08, 2.450]
#initial_steps = [0.1, 0.1, 0.007, 0.007]

# CH1E      0.0       -0.115    2.133     0.0   0.0   2.133   ! CH  (sp2) but-2-ene
# CH2E      0.0       -0.118    2.192     0.0   0.0   2.192   ! CH2 (sp3) butane
# CH3E      0.0       -0.175    2.192     0.0   0.0   2.192   ! CH3 (sp3) butane
opt = nlopt.opt(nlopt.LN_SBPLX, 2)
opt.set_lower_bounds(lower_bounds)
opt.set_upper_bounds(upper_bounds)
#opt.set_initial_step(initial_steps)
opt.set_min_objective(objfunc)
opt.set_xtol_rel(0.0002)
#opt.set_maxeval(1)


# create empty (but containing exp data) csv
for system in syslist:
    system.CreateTableWithExp()

x = opt.optimize(startpars)
minf = opt.last_optimum_value()
print(x, minf)

# end of script
