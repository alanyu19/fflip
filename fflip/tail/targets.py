#!/usr/bin/env python

from __future__ import division, print_function


from TrainingTarget import TrainingTarget
# Taken from ForceBalance

driver_path = "/lustre/alanyu17/optim/unsat/v2/"

syslist = []
syslist.append(TrainingTarget(name = 'hexene', temp = 298.15, psf = 'hexene.psf', 
                              number_molecules = 128, molmass = 84.162, dcd = 3, first_dcd = 3,
                              last_dcd = 3, replica =1, addcalc = ['hov',],
                              prop = ['rho', 'hov'], weights = [1,1],
                              root_dir = driver_path, hovnb = 15, fftx = 36 ,
                              resname = 'C6C2', crd = 'hexene.crd', guess_box = 33,                       
                              exp_dic = {'rho':682.5 , 'hov':7.52 , 'visc': 999}))

syslist.append(TrainingTarget(name = 'hexene', temp = 293.15, psf = 'hexene.psf', 
                              number_molecules = 128, molmass = 84.162, dcd = 3, first_dcd = 3,
                              last_dcd = 3, replica =1, addcalc = ['hov',],
                              prop = ['rho', 'hov'], weights = [1,1],
                              root_dir = driver_path, hovnb = 15, fftx = 36 ,
                              resname = 'C6C2', crd = 'hexene.crd', guess_box = 33,                       
                              exp_dic = {'rho':687.2 , 'hov':7.70 , 'visc': 999}))

syslist.append(TrainingTarget(name = 'undecene', temp = 293.15, psf = 'undecene.psf', 
                              number_molecules = 256, molmass = 154.292, dcd = 3, first_dcd = 3,
                              last_dcd = 3, replica =1, addcalc = ['hov',],
                              prop = ['rho'], weights = [1],
                              root_dir = driver_path, hovnb = 22, fftx = 48,
                              resname = 'T1C5', crd = 'undecene.crd', guess_box = 42,                       
                              exp_dic = {'rho':753.7 , 'hov':999 , 'visc':999}))

syslist.append(TrainingTarget(name = 'undecene', temp = 348, psf = 'undecene.psf', 
                              number_molecules = 256, molmass = 154.292, dcd = 3, first_dcd = 3,
                              last_dcd = 3, replica =1, addcalc = ['hov',],
                              prop = ['hov'], weights = [1],
                              root_dir = driver_path, hovnb = 22, fftx = 48,
                              resname = 'T1C5', crd = 'undecene.crd', guess_box = 44,                       
                              exp_dic = {'rho':999, 'hov':12.28, 'visc':999}))
