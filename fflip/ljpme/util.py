#!/usr/bin/env python
# coding: utf-8

import os
import math
import glob
import numpy as np


class ReweightingError(Exception):
    pass


class NoSimDataError(ReweightingError):
    def __init__(self):
        print(
            'There is no simulated data, please run the reweight method first!'
        )


class NoRewDataError(ReweightingError):
    def __init__(self):
        print(
            'There is no reweighted data, please run the reweight method first!'
        )


class OtherReweightError(ReweightingError):
    def __int__(self, message):
        print(message)


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

def get_rdf_pf_names(
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


def get_rdf_rmsd_names(file_template):
    def not_in(a, bb):
        for b in bb:
            if a in b:
                return False
        return True
    rdf_rmsd_names = []
    file_names = glob.glob(file_template)
    for name in file_names:
        name_parts = name.strip().split('/')[-1].split('-')
        rdf_name = name_parts[1] + '-' + name_parts[2]
        if not_in(rdf_name, rdf_rmsd_names):  # and 'Os2' not in rdf_name:
            rdf_rmsd_names.append(rdf_name + '_rmsd')
    return rdf_rmsd_names


def rename_row_col(names):
    dictt = {}
    for i, name in enumerate(names):
        dictt[i] = name
    return dictt


def move_traj(fromdir, todir):
    # Use with care
    os.system("mv {}/dyn?.dcd {}".format(fromdir, todir))
    os.system("mv {}/dyn??.dcd {}".format(fromdir, todir))
    os.system("mv {}/dyn???.dcd {}".format(fromdir, todir))


def construct_indexes(inp, properties):
    if inp=='all' or inp=='All':
        return list(range(len(properties)))
    elif ',' in inp:
        indexes = []
        values = inp.split(',')
        for v in values:
            assert 0 <= int(v) < len(properties)
            indexes.append(int(v))
        return indexes
    elif '-' in inp:
        first_last = inp.split('-')
        assert len(first_last) == 2
        first = int(first_last[0])
        last = int(first_last[1])
        assert len(properties) > last > first >= 0
        return list(range(first, last + 1))
    else:
        try:
            index = int(inp)
            assert 0 <= index < len(properties)
            return [index]
        except:
            raise Exception('Invalid property index(es)!')


def parse_first_last(inp):
    if ',' in inp:
        indexes = []
        values = inp.split(',')
        for v in values:
            indexes.append(int(v))
        return indexes
    else:
        try:
            index = int(inp)
            return [index]
        except:
            raise Exception('Invalid trajectory index(es)!')


def dict_to_str(dictionary):
    string = ""
    for key, value in dictionary.items():
        string += "{0:10s}: {1:10s}\n".format(str(key), str(value))
    return string
