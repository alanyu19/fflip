# -*- coding: utf-8 -*-

import numpy as np
from pymbar import timeseries

def find_block_size(ineff_in_step, step_size, block_length_in_ns = 5):
    """
    Find the closest length of time (in 5 ns) to the statistical ieffeciency 'g'
    'step_size' should be in nanoseconds   
    Return the block size in ns and the statistical inefficiency
    """
    ineff_in_ns = ineff_in_step * step_size
    # print("Statistical inefficiency is {0:.2f} ns".format(ineff_in_ns))
    i = block_length_in_ns
    while ineff_in_ns > i:
        i += block_length_in_ns
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

