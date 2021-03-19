#!/usr/bin/env python
# coding: utf-8

import os
from fflip.omm.paragroup import *


def save_nb_change_log(gs, offsets, file="./gtcnp_order.txt"):
    if os.path.isfile(file):
        os.system("rm {}".format(file))
    for g, off in zip(gs, offsets):
        with open(file, 'a') as f:
            f.write(
                "{0:>5} {1:>9} {2:>10.4f}\n".format(
                    g.center_names[0], g.par_type, off
                )
            )


# new for DRUDE
def save_offsets(groups, offsets, file="./offsets.log"):
    if os.path.isfile(file):
        os.system("rm {}".format(file))
    for g, off in zip(groups, offsets):
        with open(file, 'a') as f:
            f.write(
                "{0:>9} {1:>6} {2:>9} {3:>12.7f}\n".format(
                    g.lipid, g.center_names[0], g.par_type, off
                )
            )


def get_parameter_set_and_offset_by_index(index_, lip_, percentage_):
    parameter_sets_ = lip_.parse_gtcnp()
    pset_ = [parameter_sets_[index_-1]]
    # keep this list format for the other function/class (parameterEnergy?)
    offset_ = [
        gen_sensitivity_offset(
            parameter_sets_[index_-1], percentage=percentage_
        )
    ]  # same here
    return pset_, offset_


def get_one_group_with_offset(index_, lip_, percentage_, id_allowed_):
    parameter_sets_ = lip_.parse_groups(id_allowed=id_allowed_)
    pset_ = [parameter_sets_[index_-1]]
    # keep this list format for the other function/class (parameterEnergy?)
    offset_ = [
        gen_sensitivity_offset(
            parameter_sets_[index_-1], percentage=percentage_
        )
    ]  # same here
    return pset_, offset_
