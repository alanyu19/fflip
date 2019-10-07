# -*- coding: utf-8 -*-

from __future__ import division, print_function


import os
import time
import tempfile
import shutil
import numpy as np

def check_dir(fname):
    if not os.path.isdir("./{}".format(fname)):
        os.system("mkdir {}".format(fname))


def bzavg(obs,boltz):
    # Get the Boltzmann average of an observable.
    if obs.ndim == 2:          # ndim is the dimension of data set, a method of numpy
        if obs.shape[0] == len(boltz) and obs.shape[1] == len(boltz):   # number of rows and columns """
            raise Exception('Error - both dimensions have length equal to number of snapshots, now confused!')
        elif obs.shape[0] == len(boltz):
            return np.sum(obs*boltz.reshape(-1,1),axis=0)/np.sum(boltz)  # now the average is got
        elif obs.shape[1] == len(boltz):
            return np.sum(obs*boltz,axis=1)/np.sum(boltz)
        else:
            raise Exception('The dimensions are wrong!')  # when none of them has the correct length
    elif obs.ndim == 1:
        return np.dot(obs,boltz)/sum(boltz)
    else:
        raise Exception('The number of dimensions can only be 1 or 2!')


def calc_kappa(b=None, **kwargs):
    if b is None: b = np.ones(kwargs["nframes"], dtype=float)
    if 'v_' in kwargs:
        v_ = kwargs['v_']
        # return bar_unit / kT * (bzavg(v_**2,b)-bzavg(v_,b)**2)/bzavg(v_,b)
        return (bzavg(v_**2,b)-bzavg(v_,b)**2)/bzavg(v_,b)*1.0e-30/ kwargs["kT"]


def replace(source_file_path, linen, substring):  # replace the target line with the given substring
    fh, target_file_path = tempfile.mkstemp()
    with open(target_file_path, "w") as target_file:
         with open(source_file_path, "r") as source_file:
            for i, line in enumerate(source_file):
                if i == linen - 1:
                    target_file.write(line.replace(line, substring))
                else:
                    target_file.write(line.replace(line, line))
    os.rename(source_file_path, source_file_path + "-last")
    shutil.move(target_file_path, source_file_path) # move the temprary to the source file path


def replace2(filename, text_to_search, replacement_text):
    import fileinput
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(text_to_search, replacement_text), end='')


def fit_dihedral(dp, counter):
    """
    open the C5 AND C6 dihedral fitting folders and fix the dihedral force constants
    """
    ## write this into a funtion and import it as we might not need it for other parameterizing process ...
    os.chdir(dp + "c5_fitting")
    check_dir("olds") 
    os.system("cp dihe_2223.str ./olds/iter{}.str".format(counter))
    os.system("rm -f done.fit *.out *.mme *.ene *.dat")
    os.system("sbatch fit.csh")
    while not os.path.isfile("done.fit"):  # check W's minimize .sh to see where this file is generated
        time.sleep(6)
    os.chdir(dp + "c7_fitting")
    check_dir("olds") 
    os.system("cp dihe_2222.str ./olds/iter{}.str".format(counter))
    os.system("rm -f done.fit *.out *.mme *.ene *.dat")
    os.system("sbatch fit.csh")
    while not os.path.isfile("done.fit"):  # check W's minimize .sh to see where this file is generated
        time.sleep(6)

def fit_dihedral_2d(dp, counter):
    # if we stream the torsion parameters independently, we might want to delete info in the above stream file.
    # TODO: the path should be one of the input?
    os.chdir(dp + "2d_fitting")
    check_dir("olds")
    os.system("cp dihe_11222.str ./olds/iter{}.str".format(counter))
    """
    open the 2d fitting folders
    """
    ## write this into a funtion and import it as we might not need it for other parameterizing process ...
    os.system("rm -f done.fit *.out *.mme *.ene *.dat")
    os.system("sbatch fit.csh")
    while not os.path.isfile("done.fit"):  # check W's minimize .sh to see where this file is generated
        time.sleep(30)


def fit_dihedrals(driver_path, dimension, line_numbers, prm_names, x, counter):
    parameter = {}
    for pn, p in zip(prm_names, x):
        parameter[pn] = p
    if dimension == 6:
        os.chdir(driver_path + "c5_fitting")
        os.system("rm -r olds")
        os.chdir(driver_path + "c7_fitting")
        os.system("rm -r olds")
        os.chdir(driver_path + "2d_fitting")
        os.system("rm -r olds")
        substringch_1 = "CH1E\t0.0\t%.5f\t%.4f\n" % (parameter['CH1E_epsilon'], parameter['CH1E_sigma'])
        substringch_2 = "CH2E\t0.0\t%.5f\t%.4f\n" % (parameter['CH2E_epsilon'], parameter['CH2E_sigma'])
        substringch_3 = "CH3E\t0.0\t%.5f\t%.4f\n" % (parameter['CH3E_epsilon'], parameter['CH3E_sigma'])
        os.chdir(driver_path + "toppar")
        replace(driver_path + "toppar/c36ua.str", line_numbers['CH1E'], substringch_1)
        replace(driver_path + "toppar/c36ua.str", line_numbers['CH2E'], substringch_2)
        replace(driver_path + "toppar/c36ua.str", line_numbers['CH3E'], substringch_3)
        fit_dihedral(driver_path, substringch_2, substringch_3, counter)
        fit_dihedral_2d(driver_path, counter)
    elif dimension == 4:
        parameter = {}
        for pn, p in zip(prm_names, x):
            parameter[pn] = p
        os.chdir(driver_path + "c5_fitting")
        os.system("rm -r olds")
        os.chdir(driver_path + "c7_fitting")
        os.system("rm -r olds")
        substringch_2 = "CH2E\t0.0\t%.5f\t%.4f\n" % (parameter['CH2E_epsilon'], parameter['CH2E_sigma'])
        substringch_3 = "CH3E\t0.0\t%.5f\t%.4f\n" % (parameter['CH3E_epsilon'], parameter['CH3E_sigma'])
        os.chdir(driver_path + "toppar")
        replace(driver_path + "toppar/c36ua.str", line_numbers['CH2E'], substringch_2)
        replace(driver_path + "toppar/c36ua.str", line_numbers['CH3E'], substringch_3)
        fit_dihedral(driver_path, counter)


def applyconstraint(names, values, lb, up):
    for n, v, l, u in zip(names, values, lb, up):
        if not l < v < u:
            return True
    return False
