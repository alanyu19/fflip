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


def get_equil_data(data, step_size, nskip, block_size=None):
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


def calc_block_avg(data, timestep, blocksize, time_to_skip=0, verbose=True):
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
    if block_number == 0:
        raise Exception("No sufficient data. If you have tried forcing" + \
            " the calculation, that means you need longer simulation!")
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


def construct_scd_bond(topology):
    bond_dict_res = {}
    for bond in topology.bonds:
        atom1 = bond[0]
        atom2 = bond[1]
        residue = atom1.residue.name
        atomname1 = atom1.name
        atomname2 = atom2.name
        # loose condition, since 'C' and 'H' might be used as indexing other than element
        if 'C' in atomname1 and 'H' in atomname2 or \
        'H' in atomname1 and 'C' in atomname2:
            if 'C' in atomname1 and 'H' in atomname2:
                carbon = atomname1; hydrogen = atomname2
            else:
                carbon = atomname2; hydrogen = atomname1
            if not residue in bond_dict_res:
                bond_dict_res[residue] = [(carbon, hydrogen)]
            else:
                if not (atomname1, atomname2) in bond_dict_res[residue]:
                    bond_dict_res[residue].append((atomname1, atomname2))
    return bond_dict_res


def gen_scd_pairs(scd_res_dict, special_carbons_for_residues={}):
    return_dict = {}
    for res in scd_res_dict:
        res_list = []
        carbons = []
        if not res in special_carbons_for_residues:
            # first find all the carbons
            for bond in scd_res_dict[res]:
                if not bond[0] in carbons:
                    carbons.append(bond[0])
            # then generate the scd pairs
            for carbon in carbons:
                hydrogens_bond_to_the_carbon = []
                for bond in scd_res_dict[res]:
                    if bond[0] == carbon:
                        hydrogens_bond_to_the_carbon.append(bond[1])
                res_list.append((carbon, hydrogens_bond_to_the_carbon))
        else:
            # first find all the carbons
            for bond in scd_res_dict[res]:
                if not bond[0] in carbons and not bond[0] in special_carbons_for_residues[res]:
                    carbons.append(bond[0])
            for carbon in carbons:
                hydrogens_bond_to_the_carbon = []
                for bond in scd_res_dict[res]:
                    if bond[0] == carbon:
                        hydrogens_bond_to_the_carbon.append(bond[1])
                res_list.append((carbon, hydrogens_bond_to_the_carbon))
            for carbon in special_carbons_for_residues[res]:
                for bond in scd_res_dict[res]:
                    if bond[0] == carbon:
                        res_list.append((carbon, [bond[1]]))
        return_dict[res] = res_list
    return return_dict


def hard_up_low_atoms(atoms, recipe):
    if 'explicit_grouping' in recipe:
        groups = recipe['explicit_grouping']
        assert isinstance(groups, tuple)
        assert len(groups) == 2
        assert isinstance(groups[0], list)
        assert isinstance(groups[1], list)
        upper_atoms = []
        lower_atoms = []
        for atom in atoms:
            if atom in groups[0]:
                upper_atoms.append(atom)
            else:
                assert atom in groups[1]
                lower_atoms.append(atom)
        return upper_atoms, lower_atoms
    elif 'hard_border' in recipe:
        assert isinstance(recipe['hard_border'], int)
        upper_atoms = []
        lower_atoms = []
        for atom in atoms:
            if atom < recipe['hard_border']:
                upper_atoms.append(atom)
            else:
                lower_atoms.append(atom)
        return upper_atoms, lower_atoms
    else:
        pass  # being stupid for now


def manually_select_res_atom(topology_file, atom_selection):
    """
    Only support selection words like:
        resname * and (name * or name * or ... or name *) .OR.
        resname * and name *
    Args:
        topology_file: psf_file for CHARMM
        atom_selection: selection words
    Returns:
        0 based atom indexes as a numpy.array
    """
    if "\'" in atom_selection:
        # Do our own selection to avoid error in mdtraj for H3' and O3' in
        # sterols and possible future bad naming of atoms
        select_words = atom_selection.split()
        atom_names = list()
        for iw, w in enumerate(select_words):
            if 'resname' in w:
                res_name = select_words[iw + 1]
            if 'name' in w and 'resname' not in w:
                atom_names.append(select_words[iw + 1].strip(')'))
        atom_ids = []
        with open(topology_file, 'r') as topfile:
            lines = topfile.readlines()
        for line in lines:
            if res_name in line:
                for atom in atom_names:
                    if atom.upper() in line:
                        # mdtraj uses python indexing (starting from 0)
                        atom_ids.append(int(line.strip().split()[0]) - 1)
        return np.array(atom_ids)


def res_leader(resname):
    if resname in ['PLPC', 'DLIPE', 'PNPG', 'PLPI', 'DLIPS']:
        return 'P'
    elif 'CET' in resname:
        return 'NF'
    elif 'SITO' in resname:
        return 'O3'


def find_up_low_from_crd(crd_file):
    with open(crd_file, 'r') as f:
        lines = f.readlines()
    reading = False
    upper = []
    lower = []
    for l_number, l in enumerate(lines):
        if reading == True:
            elements = l.strip().split()
            atomindex = elements[0]
            resname = elements[2]
            atomname = elements[3]
            zcoor = float(elements[6])
            segname = elements[7]
            resid = int(elements[8])
            if atomname == res_leader(resname):
                if zcoor > 0:
                    if 'GLP' in segname:
                        upper.append(segname)
                    else:
                        upper.append(resname + ' ' + str(resid))
                else:
                    if 'GLP' in segname:
                        lower.append(segname)
                    else:
                        lower.append(resname + ' ' + str(resid))
        if 'EXT' in l:
            reading = True

    reading = False
    upper_atoms = []
    lower_atoms = []
    for l_number, l in enumerate(lines):
        if reading == True:
            elements = l.strip().split()
            atomindex = elements[0]
            resname = elements[2]
            atomname = elements[3]
            zcoor = float(elements[6])
            segname = elements[7]
            resid = int(elements[8])
            if 'GLP' in segname:
                if segname in upper:
                    upper_atoms.append(int(atomindex))
                elif segname in lower:
                    lower_atoms.append(int(atomindex))
            else:
                if resname + ' ' + str(resid) in upper:
                    upper_atoms.append(int(atomindex))
                elif resname + ' ' + str(resid) in lower:
                    lower_atoms.append(int(atomindex))
        if 'EXT' in l:
            reading = True
    return upper_atoms, lower_atoms


def select_atoms(topology_from, sel):
    """
    A short helper function to enable selection via atom ids or selection strings.

    Args:
        topology_from (mdtraj.Trajectory or mdtraj.Topology): The object defining the topology.
        sel: Either a selection string or a list of atom ids.

    Returns:
        list of int: Selected atom ids.
    """
    if hasattr(topology_from, "topology"):
        topology = topology_from.topology
    else:
        topology = topology_from
    if sel is None:
        return []
    elif isinstance(sel, str):
        return topology.select(sel)
    else:
        return sel


def manually_select_atoms(topology_file, atom_name):
    if "\'" in atom_name:
        # Do our own selection to avoid error in mdtraj for H3' and O3' in
        # sterols and possible future bad naming of atoms
        atmindxs = []
        with open(topology_file, 'r') as topfile:
            lines = topfile.readlines()
        for line in lines:
            if atom_name.upper() in line:
                atmindxs.append(int(line.strip().split()[0]) - 1)
        selection = atmindxs
        return selection
