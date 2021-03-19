# -*- coding: utf-8 -*-

"""Utility and helper functions for OpenMM"""

from __future__ import absolute_import, division, print_function

import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

import mdtraj as md

from pymbar import timeseries

from fflip.omm.exceptions import *


# ---------------------------- File Handling ------------------------------
def joinTraj(Traj, BlockSize):
    """
    Use mdtraj to join short trajectories by blocksize and return long trajectory(ies);
    Return a list of long trajectories and the frame number of (each) long trajectory.
    (Not used in the new version, might be useful in the future)
    Args:
        Traj: list of short trajectories
        BlockSize: in nanoseconds
    Returns:
        long trajectories (list), the length of each trajectory
    """
    BlockSize = int(BlockSize)
    LongTrajs = []
    NumShort = len(Traj)
    ''' Test if have sufficient trajectory files (compared to block size) '''
    if (NumShort < BlockSize):
        raise TrajNumberInsufficientError(NumShort, BlockSize)
    ''' Test if number of short trajectories is exactly multiple of block size'''
    if (NumShort%BlockSize) != 0:
        print("Warning: Number of trajectory is not exactly multiple of block size")
    NumBlocks = int(NumShort/BlockSize)
    for index in range(NumBlocks):
        LongTrajs.append(md.join(iter(Traj[index*BlockSize:(index+1)*BlockSize]),
                         check_topology = True, discard_overlapping_frames = True))
    ''' Length (frames) of each long trajectory '''
    TrajLength = LongTrajs[0].n_frames
    ''' Return '''
    return LongTrajs, TrajLength


# ---------------------------- Statistical & Thermo ------------------------
def Matropolis(x):
    target = np.array(x)
    ones = np.ones(len(x))
    value = np.minimum(np.exp(-target), ones)
    return value


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
    'step_size' should be in nanoseconds
    Return the data to use, the whole equil data and the equil start time
    """
    # print("Length of simulation is {} ns".format(full_time))
    [t0, g, Neff_max] = timeseries.detectEquilibration(data, nskip=nskip)
    # print("Equil reached after {0:.2f} ns".format(t0*step_size))
    eq_start_at = t0 * step_size
    data_equil = data[t0:]
    if block_size == None:
        block_size, ineff = find_block_size(g, step_size)
    # print("Using block_size of {} ns".format(block_size))
    blocks = int(np.shape(data_equil)[0] * step_size / block_size)
    use_last_steps = int(blocks * block_size / step_size)
    data_clean = data_equil[-use_last_steps:]
    return data_clean, data_equil, eq_start_at, use_last_steps, block_size


def get_prop_blockTime(option, start, end, timestep, data, file_to_write,
                       file_to_plot, block_size, cut=0):
    """
    Get Block Averaged Surface Area.
    'start', 'end' should be in nanoseconds, same for 'timestep' and 'block_size'.
    'cut' is used in getting the block average.
    If don't know the start and end, just put -1 there
    """
    option = int(option)

    if end == -1:
        print("End time unknow, calculating it based on start time ...")
        if start == -1:
            print("Start time also unknow, supposing it is 0")
            start = 0
        end = start + FindSimLength(data, timestep)
    else:
        if start == -1:
            print("Start time unknow, calculationg it based on end time ...")
            start = end - FindSimLength(data, timestep)

    ''' Check if start and end make sense '''
    if start >= end or start < 0:
        raise StartEndTimeError(start, end)

    if option == 1 or option == 3:
        n_steps = int((end - start) / timestep)
        ''' Check if data dimensions match '''
        if n_steps != np.shape(data)[0]:
            raise xyDimensionNotEqualError(n_steps, np.shape(data)[0])

        ''' The time axis '''
        times = np.linspace(start, end, int((end - start) / timestep))
        ''' Plot '''
        plot_with_accumulative(times, data, file_to_plot, numbin=50)

    if option == 2 or option == 3:
        if os.path.isfile(file_to_write):
            print(
                "Warning from 'GetPropBlockTime()': file '{}' existed, deleting ...".format(
                    file_to_write))
            os.system("rm {}".format(file_to_write))

        ''' Block AVG '''
        block_steps = int(block_size / timestep)
        get_avg_properties_using_block_size(data, file_to_write, block_size,
                                            block_steps, cut)

    if option > 3 or option <= 0:
        raise GetPropOptionError(option)


def plot_with_accumulative(x, y, file, numbin=50):
    acc = np.zeros(np.shape(y)[0])
    for i in range(np.shape(y)[0]):
        acc[i] = np.mean(y[:i + 1])
    plt.plot(x, y)
    plt.plot(x, acc, marker='o', markersize=5)
    plt.title(label=file.split(".")[0])
    plt.savefig(file)
    plt.close()


def get_avg_properties_using_block_size(raw, file, block_size, block_steps=1000,
                                        cut=0):
    data = raw[cut:]
    num_blocks = int((np.shape(data)[0] + 1) / block_steps)
    total_steps = np.shape(data)[0]
    if total_steps % block_steps != 0:
        print("Warning: block_steps {} (block_size {}) is imperfect, \
              can't divide total steps {} by block_steps!".format(block_steps,
                                                                  block_size,
                                                                  total_steps))
    means = []
    for i in range(num_blocks):
        part = data[int(i * block_steps):int((i + 1) * block_steps)]
        mean = np.mean(part)
        std = np.std(part)
        prop = file.split(".")[0]
        with open(file, "a") as f:
            if f.tell() == 0:
                print("A new file is created ...")
                f.write("(avg){0:>14s}{1:>13s}".format(
                    prop, "STD/STE")
                )
                f.write("\n")
            else:
                pass  # print("File existed, appending")
            f.write("{0:>19.3f}{1:>13.3f}".format(
                mean, std))
            f.write("\n")
            means.append(mean)
    with open(file, "a") as f:
        f.write("{0:>8s}".format("Overall"))
        f.write("{0:>11.3f}{1:>13.3f}".format(
            np.mean(np.array(means)),
            np.std(np.array(means)) / np.sqrt(num_blocks)))
        f.write("\n")
    return np.mean(np.array(means)), np.std(np.array(means)) / np.sqrt(
        num_blocks)


def FindSimLength(data, step_size):
    """
    Find the time length of the data, 'step_size' should be in nanosecond.
    """
    n_step = np.shape(data)[0]
    duration = n_step * step_size
    return int(duration)


def get_histogram_with_reference(a, b, axis_labels, file_to_save, numbin=50):
    ''' reference is in blue circle, major is in red plus '''
    for i in range(len(a)):
        p = np.histogram(a[i], bins=numbin, density=True)
        plt.plot(p[1][:-1], p[0], 'bo')
    for i in range(len(b)):
        p = np.histogram(b[i], bins=numbin, density=True)
        plt.plot(p[1][:-1], p[0], 'r+')
    title = file_to_save.split(".")[0]
    plt.title(label=title)
    plt.xlabel(xlabel=axis_labels[0])
    plt.ylabel(ylabel=axis_labels[1])
    plt.savefig(file_to_save)
    plt.show()
    plt.close()


def get_2dlog_histogram(a, b, axis_labels, file_to_save, numbin=50):
    fig, ax = plt.subplots(tight_layout=True)
    ax.hist2d(a, b, bins=numbin, norm=colors.LogNorm())
    name = file_to_save.split(".")[0]
    plt.title(name)
    plt.xlabel(xlabel=axis_labels[0])
    plt.ylabel(ylabel=axis_labels[1])
    plt.savefig(file_to_save)
    plt.close()


# --------------------------- OpenMM Iteractions ---------------------------

def write_energies_into_files(energies, names, par_types, old_params, new_params, file_template = "./area_energies/{}_{}_{:.4f}_{:.4f}.dat"):
    """
    :param energies: array/list of arrays
    :param names: list
    :param par_types: list
    :param old_params: list/array
    :param new_params: list/array
    :param file_template:
    :return: None
    """
    for energy_data, name, par_type, old, new in zip(energies, names, par_types, old_params, new_params):
        if not os.path.isdir(file_template.split("{}")[0]):
            os.system("mkdir {}".format(file_template.split("{}")[0]))
        if os.path.isfile(file_template.format(name, par_type, old, new)) == False:
            np.savetxt(file_template.format(name, par_type, old, new), energy_data)
        else:
            data_last_step = np.loadtxt(file_template.format(name, par_type, old, new))
            to_save = np.append(data_last_step, energy_data)
            np.savetxt(file_template.format(name, par_type, old, new), to_save)


def write_area_into_file(areas, file_template = "./area_energies/areas.dat"):
    """
    :param areas:
    :param file_template:
    :return: None
    """
    if not os.path.isdir(file_template.split("areas.")[0]):
        os.system("mkdir {}".format(file_template.split("areas.")[0]))
    if os.path.isfile(file_template) == False:
        np.savetxt(file_template, areas)
    else:
        data_last_step = np.loadtxt(file_template)
        to_save = np.append(data_last_step, areas)
        np.savetxt(file_template, to_save)

