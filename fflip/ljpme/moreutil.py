# -*- coding: utf-8 -*-

import time
import math
import glob
import pandas as pd

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


def get_avail_exp_prop_names(file_template = '/u/alanyu/c36ljpme/fflow/exp/*.exp'):
    exp_props = glob.glob(file_template)
    useful_props = []
    for prop in exp_props:
        name = prop.split('/')[-1].strip().split('.')[0]
        useful_props.append(name)
    return useful_props


def get_sim_scd_names(file_template = '/u/alanyu/c36ljpme/fflow/runner/scd_dppc/block_data/*'):
    scd_names = []
    file_names = glob.glob(file_template)
    for name in file_names:
        scd_name = name.strip().split('/')[-1].split('_')[0]
        if scd_name not in scd_names:
            scd_names.append(scd_name)
    return scd_names


def get_rdf_names_as_properties(file_template = '/u/alanyu/c36ljpme/fflow/runner/rdf/block_data/sparse*'):
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
    dict = {}
    for i, name in enumerate(names):
        dict[i] = name
    return dict


def order_peak_foot(name):
    splist = name.split('_')
    if splist[1] == 'peak':
        b = 0
    elif splist[1] == 'foot':
        b = 1
    else:
        raise OtherReweightError('Name of property () is not suitable'.format(name))
    a = int(splist[2])
    return 2*(a-1) + b


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def find_peak_position(x, first_n_peak = 2):
    local_max_bool = np.r_[True, x[1:] > x[:-1]] & np.r_[x[:-1] > x[1:], True]
    _where_max = np.where(local_max_bool == True)
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


def find_foot_position(x, first_n_foots = 1):
    local_min_bool = np.r_[True, x[1:] < x[:-1]] & np.r_[x[:-1] < x[1:], True]
    _where_min = np.where(local_min_bool == True)
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


def find_rdf_peaks_and_foots(r, rdf, first_n_peaks = 2, first_n_foots = 1, smooth_window_size = 1):
    """
    Args:
        r: numpy.array, (n,), radius of the rdf
        rdf: numpy.array, (m,n), the rdf
        first_n_peaks: integer, number of the nearest peaks to find ()
    """
    r = np.array(r)
    rdf = np.array(rdf)
    shape_original = np.shape(rdf)
    if len(shape_original) == 2:
        value_list = []
        r_list = []
        for i in range(int(shape_original[0])):
            copy_of_rdf = rdf[i, :]
            peak_indexes_smoothed = find_peak_position(smooth(copy_of_rdf, smooth_window_size), first_n_peaks)
            foot_indexes_smoothed = find_foot_position(smooth(copy_of_rdf, smooth_window_size), first_n_foots)
            peak_indexes = find_peak_position(copy_of_rdf, first_n_peaks)
            foot_indexes = find_foot_position(copy_of_rdf, first_n_foots)
            r_list.append(r[np.array(peak_indexes_smoothed + foot_indexes_smoothed)])
            value_list.append(copy_of_rdf[np.array(peak_indexes + foot_indexes)])
        return r_list, value_list
    else:
        peak_indexes_smoothed = find_peak_position(smooth(rdf, smooth_window_size), first_n_peaks)
        foot_indexes_smoothed = find_foot_position(smooth(rdf, smooth_window_size), first_n_foots)
        peak_indexes = find_peak_position(rdf, first_n_peaks)
        foot_indexes = find_foot_position(rdf, first_n_foots)
        return r[np.array(peak_indexes_smoothed + foot_indexes_smoothed)], rdf[np.array(peak_indexes + foot_indexes)]


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

# ----------------------------- Guess Functions ------------------------------

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


def make_guess_of_layertype(name):
    if 'bilayer' in name:
        return 'bilayer'
    elif 'mono' in name:
        return 'monolayer'
    elif 'scd' in name:
        return 'bilayer'
    else:
        return 'bulk'


# ------------------------------ Delete Later ----------------------------
def gen_dppc_bilayer_sim(
        iteration, surface_tension, change_para, solution_file=None,
        temperature=323.15, last_seqno=100,
        where='/gs-scratch/mbs/alanyu/c36ljpme',
        overwrite=False):
    md_options = {'surface_tension': surface_tension,
                  'boxx': 48,
                  'boxz': 68,
                  'temperature': temperature,
                  'zmode': 0,
                  'barostat': 'MCM',
                  'change_para': change_para
                  }
    job_dir = os.path.join(where, 'iter{}/dppc_bilayer_{}_{}'.format(iteration, surface_tension, temperature))
    # sim = dppc_bilayer_calc.gensim(job_dir,
    #                                last_seqno,
    #                                md_options=md_options,
    #                                solution_file=solution_file,
    #                                overwrite=overwrite,
    #                                template='sim_template')
    # return sim


def gen_dmpc_bilayer_sim(iteration, surface_tension, change_para,
                         solution_file=None, temperature=303.15, last_seqno=100,
                         where='/gs-scratch/mbs/alanyu/c36ljpme',
                         overwrite=False):
    # md_options = {'surface_tension': surface_tension,
    #               'boxx': 48,
    #               'boxz': 68,
    pass


def gen_dopc_bilayer_sim(iteration, surface_tension, change_para,
                         solution_file=None, temperature=303.1, last_seqno=100,
                         where='/gs-scratch/mbs/alanyu/c36ljpme',
                         overwrite=False):
    # md_options = {'surface_tension': surface_tension,
    #              'boxx': 50,
    #              'boxz': 76,
    pass


def gen_dppc_monolayer_sim(iteration, surface_tension, change_para,
                           solution_file=None, temperature=321.15,
                           last_seqno=120,
                           where='/gs-scratch/mbs/alanyu/c36ljpme',
                           overwrite=False):
    # md_options = {'surface_tension': surface_tension,
    #               'boxx': 43, 'zmode': 1?sss
    #               'boxz': 223,
    pass


def gen_prpc_sim(iteration, surface_tension, change_para, solution_file=None,
                 temperature=298.15, last_seqno=100,
                 where='/gs-scratch/mbs/alanyu/c36ljpme', overwrite=False):
    md_options = {'surface_tension': surface_tension,
                  'boxx': 43,
                  'temperature': temperature,
                  'barostat': 'MC',
                  'change_para': change_para
                  }
    job_dir = os.path.join(where, 'iter{}/prpc_bulk_{}_{}'.format(iteration, surface_tension, temperature))
    # sim = prpc_calc.gensim(job_dir,
    #                        last_seqno,
    #                        md_options=md_options,
    #                        solution_file=solution_file,
    #                        overwrite=overwrite,
    #                        template='sim_template')
    # return sim



lipfinder = {
    'dppc': dppc,
    'prpc': prpc,
    'dmpc': dmpc,
    'dopc': dopc,
    'dlpc': dppc
}


