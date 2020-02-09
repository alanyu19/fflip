# -*- coding: utf-8 -*-

import numpy as np
import os
import glob

from coffe.omm.exceptions import *
from fflip.chm import *


# --------------------------- First things first ---------------------------

def get_sign(x):
    return math.copysign(1, x)


def copy_folder(sample, destination):
    des_parent = '/'.join(destination.split('/')[:-2])
    if not os.path.isdir(des_parent):
        os.system("mkdir -p " + des_parent)
    os.system("cp -r {} {}".format(sample, destination))


def combine_solutions(solutions, param):
    addone = []
    for s in solutions:
        addone.append(0.01 * s + 1)
    current = addone[0]
    for nexts in addone[1:]:
        temp = []
        for i, j, p in zip(current, nexts, param):
            if p.par_type == 'charge':
                temp.append(i + j -1)
            else:
                temp.append(i * j)
        current = np.array(temp)
    return (current - 1) * 100


def short_layer_type(layer_type):
    if layer_type == 'bilayer':
        return 'bi'
    elif layer_type == 'monolayer':
        return 'mono'
    elif layer_type == 'bulk':
        return 'bulk'


def get_avail_exp_prop_names(
        file_template='/u/alanyu/c36ljpme/fflow/exp/*.exp'
):
    exp_props = glob.glob(file_template)
    useful_props = []
    for prop in exp_props:
        name = prop.split('/')[-1].strip().split('.')[0]
        useful_props.append(name)
    return useful_props


def get_sim_scd_names(
        file_template='/u/alanyu/c36ljpme/fflow/runner/scd_dppc/block_data/*'
):
    scd_names = []
    file_names = glob.glob(file_template)
    for name in file_names:
        scd_name = name.strip().split('/')[-1].split('_')[0]
        if scd_name not in scd_names:
            scd_names.append(scd_name)
    return scd_names


def get_rdf_names_as_properties(
        file_template='/u/alanyu/c36ljpme/fflow/runner/rdf/block_data/sparse*'
):
    def not_in(a, bb):
        for b in bb:
            if a in b:
                return False
        return True
    rdf_names = []
    file_names = glob.glob(file_template)
    for name in file_names:
        name_parts = name.strip().split('/')[-1].split('-')
        rdf_name = name_parts[1] + '-' + name_parts[2]
        if not_in(rdf_name, rdf_names) and 'Os2' not in rdf_name:
            rdf_names.append(rdf_name + '_peak_1')
            rdf_names.append(rdf_name + '_foot_1')
            if not (rdf_name == 'O2-OW' or rdf_name == 'Ob-OW'):
                rdf_names.append(rdf_name + '_peak_2')
    return rdf_names


def rename_row_col(names):
    dictt = {}
    for i, name in enumerate(names):
        dictt[i] = name
    return dictt


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
    box = np.ones(box_pts)/box_pts
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


class FolderNamingScheme(object):
    def __init__(self, target_property):
        self.tp = target_property

    def exp_folder(self):
        return os.path.join(self.tp.root_dir, 'exp')

    def robustness_dir(self):
        return os.path.join(self.tp.root_dir, "robustness")

    def reweight_dir(self):
        return os.path.join(self.tp.root_dir, "reweightings")

    def trajectory_folder(self, iteration, traj_root):
        subdir = "iter{}/{}_{}_{}_{}".format(
            iteration, self.tp.lipname, self.tp.system_type,
            self.tp.surface_tension, self.tp.temperature
        )
        return os.path.join(traj_root, subdir)

    def property_data_folder(self, iteration):
        subdir = "observables/iter{}/{}_{}_{}_{}_{}".format(
            iteration, self.tp.prop_type, self.tp.lipname, self.tp.system_type,
            self.tp.surface_tension, self.tp.temperature
        )
        return os.path.join(self.tp.root_dir, subdir)

    def potential_data_folder(self, iteration):
        subdir = "potentials/iter{}/potential_{}_{}_{}_{}_{}".format(
            iteration, self.tp.lipname, self.tp.system_type,
            self.tp.surface_tension, self.tp.temperature, self.tp.perturbation
        )
        return os.path.join(self.tp.root_dir, subdir)

    def reweighting_folder(self, iteration):
        return os.path.join(
            self.tp.reweight_dir, 'iter{}/{}'.format(
                iteration, self.tp.perturbation
            )
        )

    def reweighting_file_name(self):
        return "{}_{}_{}_{}_{}".format(
            self.tp.prop_type, self.tp.lipname, self.tp.system_type,
            self.tp.surface_tension, self.tp.temperature
        )

    def robustness_folder(self, iteration):
        return os.path.join(
            self.tp.robdir, 'iter{}/{}'.format(iteration, self.tp.perturbation)
        )

    def robustness_diff_file(self, iteration):
        return os.path.join(
            self.tp.robdir, 'iter{}/{}/{}_diff.txt'.format(
                iteration, self.tp.perturbation, self.tp.name
            )
        )

    def robustness_std_file(self, iteration):
        return os.path.join(
            self.tp.robdir, 'iter{}/{}/{}_std.txt'.format(
                iteration, self.tp.perturbation, self.tp.name
            )
        )


class SimOptScheme(object):
    def __init__(self, target_system):
        self.ts = target_system

    @property
    def last_seqno(self):
        finder = {
            "dppc_bilayer_0_323.15": 200,
            "dppc_bilayer_-5_323.15": 300,
            "dppc_bilayer_5_323.15": 300,
            "dppc_bilayer_0_333.15": 200,
            "dppc_monolayer_18_321.15": 200,
            "dppc_monolayer_40_321.15": 200,
            "dppc_monolayer_55_321.15": 200,
            "dlpc_bilayer_0_303.15": 200,
            "dmpc_bilayer_0_303.15": 200,
            "popc_bilayer_0_303.15": 200,
            "prpc_bulk_0_298.15": 100
        }
        return finder["{}_{}_{}_{}".format(
            self.ts.lipname, self.ts.system_type,
            self.ts.surface_tension, self.ts.temperature
        )]

    @property
    def boxx(self):
        finder = {
            "dppc_bilayer_0": 48,
            "dppc_bilayer_-5": 48,
            "dppc_bilayer_5": 48,
            "dppc_monolayer_18": 44,
            "dppc_monolayer_40": 48,
            "dppc_monolayer_55": 52,
            "dlpc_bilayer_0": 48,
            "dmpc_bilayer_0": 48,
            "popc_bilayer_0": 49,
            "prpc_bulk_0": 43
        }
        return finder["{}_{}_{}".format(
            self.ts.lipname, self.ts.system_type, self.ts.surface_tension
        )]

    @property
    def boxz(self):
        finder = {
            "dppc_bilayer_0": 64,
            "dppc_bilayer_-5": 64,
            "dppc_bilayer_5": 64,
            "dppc_monolayer_18": 223,
            "dppc_monolayer_40": 223,
            "dppc_monolayer_55": 223,
            "dlpc_bilayer_0": 62,
            "dmpc_bilayer_0": 63,
            "popc_bilayer_0": 65,
            "prpc_bulk_0": None
        }
        return finder["{}_{}_{}".format(
            self.ts.lipname, self.ts.system_type, self.ts.surface_tension
        )]

    @property
    def intgrt(self):
        return 'L'

    @property
    def barostat(self):
        if self.ts.system_type != 'bulk':
            return "MCM"
        else:
            return "MC"

    @property
    def zmode(self):
        if self.ts.system_type == "monolayer":
            return 1
        elif self.ts.system_type == "bilayer":
            return 0
        else:
            assert self.ts.system_type == "bulk"
            return None


# ----------------------------- Guess Functions ------------------------------
# ------------------------------- Delete Later -------------------------------

def make_guess_of_block_size(btype, prop_type):
    if btype == 0:
        if prop_type == 'rdf':
            return 1
        else:
            return 10
    elif btype == 1:
        return 10
    else:
        raise Exception("type {} not supported, use 0 or 1".format(btype))


def make_guess_of_trajectory_range(name):
    # TODO: could add a check for last trajectory (through energy_rel_dir)
    if 'scd' in name:
        return 0, 0
    elif 'area' in name and 'bilayer' in name:
        return 0, 0
    elif 'area' in name and 'mono' in name:
        return 0, 0
    else:
        return 0, 0


def make_guess_of_intervals(name):
    if 'scd' in name:
        return 10, 10
    elif 'area' in name:
        return 10, 10
    else:
        return 10, 1


def make_guess_of_scaling(name):
    if 'area' in name:
        # (nm^2)
        return 1 / 0.6
    elif 'peak' in name:
        return 1 / 1
    elif 'foot' in name:
        return 1 / 0.5
    elif 'scd' in name:
        return 1 / 0.15
    elif 'ka' in name:
        return 1 / 200
    elif 'delta_area' in name:
        return 1 / 0.02


def make_guess_of_layertype(name):
    if 'bilayer' in name:
        return 'bilayer'
    elif 'mono' in name:
        return 'monolayer'
    elif 'scd' in name:
        return 'bilayer'
    else:
        return 'bulk'


lipfinder = {
    'dppc': pc,
    'prpc': pc,
    'dmpc': pc,
    'dopc': pc,
    'popc': pc,
    'dlpc': pc
}



