#!/usr/bin/env python
# coding: utf-8


from coffe.omm.playpara import *


def select_a_list_of_atoms(atom_names, topology):
    return_list = []
    for atom_Name in atom_names:
        return_list += list(topology.select("name {}".format(atom_Name)))
    return return_list


def reassign_tail_lj_parameters(system, topology, new_parameters, sn1, sn2):
    """

    Args:
        system:
        topology:
        new_parameters:
        sn1:
        sn2:

    Returns:

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
            force, atom_id, round(
                (1/2)**(1/6) * new_parameters["ch2_epsilon"] * 0.2, 5
            )
        )
        reassign_sigma_by_id(
            force, atom_id, round(new_parameters["ch2_half_r_min"] * 4.184, 4)
        )
    for atom_id in atom_sele_ch3:
        reassign_epsilon_by_id(
            force, atom_id,  round(
                (1/2)**(1/6) * new_parameters["ch3_epsilon"] * 0.2, 5
            )
        )
        reassign_sigma_by_id(
            force, atom_id, round(new_parameters["ch3_half_r_min"] * 4.184, 4)
        )
    for atom_id in atom_sele_h2:
        reassign_epsilon_by_id(
            force, atom_id, round(
                (1/2)**(1/6) * new_parameters["h2_epsilon"] * 0.2, 5
            )
        )
        reassign_sigma_by_id(
            force, atom_id, round(new_parameters["h2_half_r_min"] * 4.184, 4)
        )
    for atom_id in atom_sele_h3:
        reassign_epsilon_by_id(
            force, atom_id,  round(
                (1/2)**(1/6) * new_parameters["h3_epsilon"] * 0.2, 5
            )
        )
        reassign_sigma_by_id(
            force, atom_id, round(new_parameters["h3_half_r_min"] * 4.184, 4)
        )
