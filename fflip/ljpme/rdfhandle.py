# -*- coding: utf-8 -*-

import numpy as np
from fflip.omm.exceptions import *


# --------------------------- First things first ---------------------------

def order_peak_foot(name):
    splist = name.split('_')
    if splist[1] == 'peak':
        b = 0
    elif splist[1] == 'foot':
        b = 1
    else:
        raise OtherReweightError(
            'Name of property () is not suitable'.format(name)
        )
    a = int(splist[2])
    return 2*(a-1) + b


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def combine_peak_foot_indexes(peak_indexes, foot_indexes):
    assert 0 <= len(peak_indexes) - len(foot_indexes) <= 1
    i = 0
    comb = []
    while i < len(peak_indexes) and i < len(foot_indexes):
        comb.append(peak_indexes[i])
        comb.append(foot_indexes[i])
        i += 1
    if i < len(peak_indexes):
        comb.append(peak_indexes[i])
    return comb


def find_peak_position(x, first_n_peak=2):
    local_max_bool = np.r_[True, x[1:] > x[:-1]] & np.r_[x[:-1] > x[1:], True]
    _where_max = np.where(local_max_bool==True)
    _where_big = np.where(x > 0.2)
    where_max = list(_where_max[0])
    where_big = list(_where_big[0])
    return_list = []
    count = 1
    previous_index = -999
    for index in where_max:
        if count > first_n_peak:
            break
        if index in where_big and (index - previous_index) > 40:
            return_list.append(index)
            previous_index = index
            count += 1
    return return_list


def find_foot_position(x, first_n_foots=1):
    local_min_bool = np.r_[True, x[1:] < x[:-1]] & np.r_[x[:-1] < x[1:], True]
    _where_min = np.where(local_min_bool==True)
    _where_big = np.where(x > 0.05)
    where_min = list(_where_min[0])
    where_big = list(_where_big[0])
    return_list = []
    count = 1
    previous_index = -999
    for index in where_min:
        if count > first_n_foots:
            break
        if index in where_big and (index - previous_index) > 40:
            return_list.append(index)
            previous_index = index
            count += 1
    return return_list


def find_rdf_peaks_and_foots(rdf, r=None, first_n_peaks=2, first_n_foots=1,
                             smooth_window_size=1, no_r=False):
    """

    Args:
        r: numpy.array, (n,), radius of the rdf
        rdf: numpy.array, (m,n), the rdf
        first_n_peaks: integer, number of the nearest peaks to find ()
        first_n_foots: integer, number of the nearest foots to find ()
        smooth_window_size: integer, as named, number of data point
        no_r: bool, if return the position, if False, only return the
    Returns:

    """
    if not no_r:
        assert r is not None
        r = np.array(r)
    rdf = np.array(rdf)
    shape_original = np.shape(rdf)
    if len(shape_original) == 2:
        value_list = []
        r_list = []
        for i in range(int(shape_original[0])):
            copy_of_rdf = rdf[i, :]
            peak_indexes_smd = find_peak_position(
                smooth(copy_of_rdf, smooth_window_size), first_n_peaks
            )
            foot_indexes_smd = find_foot_position(
                smooth(copy_of_rdf, smooth_window_size), first_n_foots
            )

            # peak_indexes = find_peak_position(copy_of_rdf, first_n_peaks)
            # foot_indexes = find_foot_position(copy_of_rdf, first_n_foots)
            indexes = combine_peak_foot_indexes(
                peak_indexes_smd, foot_indexes_smd
            )
            if not no_r:
                r_list.append(
                    r[np.array(indexes)]
                )
            value_list.append(
                copy_of_rdf[np.array(indexes)]
            )
        if not no_r:
            return r_list, value_list
        else:
            return value_list
    else:
        peak_indexes_smd = find_peak_position(
            smooth(rdf, smooth_window_size), first_n_peaks
        )
        foot_indexes_smd = find_foot_position(
            smooth(rdf, smooth_window_size), first_n_foots
        )
        # peak_indexes = find_peak_position(rdf, first_n_peaks)
        # foot_indexes = find_foot_position(rdf, first_n_foots)
        indexes = combine_peak_foot_indexes(peak_indexes_smd, foot_indexes_smd)
        if not no_r:
            return r[np.array(indexes)], rdf[np.array(indexes)]
        else:
            return rdf[np.array(indexes)]


def extract_exp(name, exp):
    if 'peak' in name.lower() or 'foot' in name.lower():
        pfs = find_rdf_peaks_and_foots(
            exp, r=None, first_n_peaks=2, first_n_foots=1,
            smooth_window_size=1, no_r=True
        )
        order = order_peak_foot(name)
        return pfs[order]
    else:
        return exp
