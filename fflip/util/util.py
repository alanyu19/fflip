# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from pymbar import timeseries

def find_block_size(ineff_in_step, step_size):
    """
    Find the closest length of time (in 5 ns) to the statistical ieffeciency 'g'
    'step_size' should be in nanoseconds   
    Return the block size in ns and the statistical inefficiency
    """
    ineff_in_ns = ineff_in_step * step_size
    # print("Statistical inefficiency is {0:.2f} ns".format(ineff_in_ns))
    i = 5
    while ineff_in_ns > i:
        i+=5
    return i, ineff_in_ns

def get_equil_data(data, step_size, nskip, block_size = None):
    """
    :param data: the original data, of area/lipid, or other properties
    :param step_size: the step size of simulation in nanoseconds
    :param nskip: the step to skip between, recommended value is 500
    (make sure it's equivalent to 0.2 ~ 1 ns)
    :param block_size: the block_size for calculating average/stderr,
    if left as None, the function will find for you (recommended)
    :return: 1) the largest possible equilibrium data set that can be divided by
    the block size; 2) all equilibrium data; 3) equilibrium starting time (in ns);
    4) (counting backward from end of simulation) steps can be used for analysis;
    5) the block size it finds or you provide (you want to use uncorrelated
    blocks for averaging ...)
    """
    [t0, g, Neff_max] = timeseries.detectEquilibration(data, nskip=nskip)
    eq_start_at = t0 * step_size
    data_equil = data[t0:]
    if block_size == None:
        block_size, ineff = find_block_size(g, step_size)
    else:
        test_block_size, ineff = find_block_size(g, step_size)
        if test_block_size > block_size:
            print("Warning: the block size provided is too small!"
                  " Changing it to {} ns".format(test_block_size))
            block_size = test_block_size
    # print("Using block_size of {} ns".format(block_size))
    blocks = int(np.shape(data_equil)[0] * step_size / block_size)
    use_last_steps = int(blocks * block_size / step_size)
    data_clean = data_equil[-use_last_steps:]
    return data_clean, data_equil, eq_start_at, use_last_steps, block_size
        
def calc_block_avg(data, timestep, blocksize, time_to_skip = 0, verbose = True):
    """
    all inputs(time-related) are in nanoseconds
    """
    if time_to_skip != 0:
        print("\nSkipping {} ns of useful data ...".format(time_to_skip))
    nskip = int(time_to_skip/timestep)
    data = data[nskip:]
    block_means=[]
    steps_per_block = int(blocksize/timestep)
    block_number = int(len(data)/steps_per_block)
    steps = block_number * steps_per_block
    data = data[-steps:]
    if verbose:
        print("Using {} blocks (last {} ns) for SA calculation ...".format(
              block_number, int(block_number*blocksize)))
    for b in range(block_number):
        avg = np.mean(data[b*steps_per_block:(b+1)*steps_per_block])
        block_means.append(avg)
    blocks = np.array(block_means)
    mean = np.mean(blocks)
    error = np.std(blocks,)/np.sqrt(block_number)
    return mean, error, block_number, nskip

def execute(nlipid, temperature, block_size, step_size, force_calc = False, force_fraction = 0.6,
            file_name = 'area.dat', data_dimension = 1, data_col = -1, multiply_factor = 100,
            time_to_skip = 10):
    data0 = np.loadtxt(file_name)
    if data_dimension == 1:
        pass
    elif data_dimension == 2:
        data0 = data[:,data_col]
    data0 = data0 * multiply_factor
    kb = 1.3806
    sad, data_equil, eq_start_at, use_last_steps, block_size = \
    get_equil_data(data0, step_size = step_size, nskip = 500, block_size = block_size)
    print("Equilibrium reached after {} ns".format(eq_start_at))
    print("Block size is {} ns".format(block_size))
    if use_last_steps < 3 * block_size * step_size or force_calc == True: 
        print("Ignoring the staring point of euiqlibrium above, using last {}% data ...".format(int(force_fraction*100)))
        print("\n=.=")
        total_length = data0.shape[0]
        laststeps = int(total_length * force_fraction)
        data_to_use = data0[-laststeps:]
    else:
        data_to_use = data_equil
    mean, err, bn, skip = calc_block_avg(data_to_use, step_size, block_size, time_to_skip = time_to_skip)
    # Calculate Compressibility Modulus using largest possible data set
    N_steps= np.shape(data_to_use)[0]     # number of data points
    print("Using last {} ns data for compressibility calculation ...".format(round(N_steps*step_size,2)))
    mean2, err2, bn2, skip2 = calc_block_avg(data_to_use, step_size, block_size, time_to_skip = 0, verbose = False)
    VarS=0      # VarS is variance after the system reaches equilibrium
    for i in range(0, N_steps, 1):
        VarS=VarS+(data_to_use[i] - mean2)**2
    Variance=VarS/N_steps
    Ka = kb * temperature * mean / (nlipid * Variance * 1000) # order is corrected/recovered using 1000
    print('Final SA_avg/lip and err: %3.2f +- %3.2f' % (mean, err))
    print('Compressibility modulus is: {0:.2f}'.format(Ka))
