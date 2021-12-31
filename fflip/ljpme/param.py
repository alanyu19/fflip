#!/usr/bin/env python
# coding: utf-8

import os
from fflip.omm.playpara import *
from fflip.omm.paragroup import *


def get_sign(x):
    return math.copysign(1, x)


def gen_param_offset(param_group, amount, sign=1):
    """
    Function to generate the nonbonded parameter offset for a group of atoms.
    Args:
        param_group: the parameter group
        amount (float): the size of the offset
        sign (int): default to 1, the sign of the offset, can be 1/-1
    Returns:
        the offset of the parameter
    """
    pt = param_group.par_type
    if pt == 'charge':
        num_atom = len(param_group.center_names)
        return amount / num_atom
    if pt == 'alpha':
        return amount
    if pt == 'thole' or pt == 'nbthole':
        return amount
    if pt == 'epsilon' or pt == 'sigma':
        return amount


def save_nb_change_log(gs, offsets, file="./nbgroup_order.txt"):
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


def get_parameter_set_and_offset_by_index(index, lipid, amount):
    """
    Generate the parameter group and corresponding offset based on user
    provided index and perturbation size (amount)
    Args:
        index (int): index of the parameter group within the lipid.
        lipid: the lipid object
        amount (float): amount of perturbation
    Returns:
        the parameter group and the offset
    """
    parameter_sets = lipid.parse_groups()
    pset = [parameter_sets[index-1]]
    offset = [gen_param_offset(parameter_sets[index-1], amount)]
    return pset, offset


def get_one_group_with_offset(index, lipid, amount, id_allowed):
    """
    Upgraded function for get_parameter_set_and_offset_by_index.

    Generate the parameter group and corresponding offset based on user
    provided index and perturbation size (amount).

    Args:
        index (int): index of the parameter group within the lipid.
        lipid: the lipid object.
        amount (float): amount of perturbation.
        id_allowed: parameter indexes allowed for perturbation.
    Returns:
        the parameter group and the offset.
    """
    parameter_sets = lipid.parse_groups(id_allowed=id_allowed)
    pset = [parameter_sets[index-1]]
    offset = [
        gen_param_offset(
            parameter_sets[index-1], amount=amount
        )
    ]
    return pset, offset


def select_a_list_of_atoms(atom_names, topology):
    return_list = []
    for atom_Name in atom_names:
        return_list += list(topology.select("name {}".format(atom_Name)))
    return return_list


def find_a_paragroup(groups, atom_name, par_type):
    succeed = False
    for g in groups:
        if g.par_type == par_type and atom_name in g.center_names:
            succeed = True
            group = g
    assert succeed
    return group


def reassign_tail_lj_parameters(system, topology, new_parameters, sn1, sn2):
    """

    Args:
        system: the openmm system
        topology: topology of the system
        new_parameters (dict): the parameters
        sn1: size (length) of sn1 chain
        sn2: size of sn2 chain

    Returns:
        the new system
    """
    ch2 = ['C3{}'.format(i) for i in range(2, sn1)] + \
          ['C2{}'.format(i) for i in range(2, sn2)]
    h2 = ['H3{}X'.format(i) for i in range(2, sn1)] + \
         ['H3{}Y'.format(i) for i in range(2, sn1)] + \
         ['H2{}R'.format(i) for i in range(2, sn2)] + \
         ['H2{}S'.format(i) for i in range(2, sn2)]
    ch3 = ['C3{}'.format(sn1), 'C2{}'.format(sn2)]
    h3 = ['C3{}X'.format(sn1), 'C3{}Y'.format(sn1), 'C3{}Z'.format(sn1),
          'C2{}R'.format(sn2), 'C2{}R'.format(sn2), 'C2{}T'.format(sn2)]

    atom_sele_ch2 = select_a_list_of_atoms(ch2, topology)
    atom_sele_h2 = select_a_list_of_atoms(h2, topology)
    atom_sele_ch3 = select_a_list_of_atoms(ch3, topology)
    atom_sele_h3 = select_a_list_of_atoms(h3, topology)

    force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            break

    assert force is not None
    for atom_id in atom_sele_ch2:
        reassign_epsilon_by_id(
            force, int(atom_id), round(
                new_parameters["ch2_epsilon"] * 4.184, 4
            )
        )
        reassign_sigma_by_id(
            force, int(atom_id), round(
                (1/2)**(1/6) * new_parameters["ch2_half_r_min"] * 0.2, 5
            )
        )
    for atom_id in atom_sele_ch3:
        reassign_epsilon_by_id(
            force, int(atom_id),  round(
                new_parameters["ch3_epsilon"] * 4.184, 4
            )
        )
        reassign_sigma_by_id(
            force, int(atom_id), round(
                (1/2)**(1/6) * new_parameters["ch3_half_r_min"] * 0.2, 5
            )
        )
    for atom_id in atom_sele_h2:
        reassign_epsilon_by_id(
            force, int(atom_id), round(
                new_parameters["h2_epsilon"] * 4.184, 4
            )
        )
        reassign_sigma_by_id(
            force, int(atom_id), round(
                (1/2)**(1/6) * new_parameters["h2_half_r_min"] * 0.2, 5
            )
        )
    for atom_id in atom_sele_h3:
        reassign_epsilon_by_id(
            force, int(atom_id),  round(
                new_parameters["h3_epsilon"] * 4.184, 4
            )
        )
        reassign_sigma_by_id(
            force, int(atom_id), round(
                (1/2)**(1/6) * new_parameters["h3_half_r_min"] * 0.2, 5
            )
        )
    return system
