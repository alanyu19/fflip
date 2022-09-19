# -*- coding: utf-8 -*-

"""Utility and helper functions for OpenMM"""

# my test change
# test change

from __future__ import absolute_import, division, print_function

import os
import copy
import numpy as np
import warnings

import simtk.unit as u
from simtk.openmm import Context, LangevinIntegrator
from simtk.openmm import NonbondedForce, CustomNonbondedForce
from simtk.openmm.app import CharmmPsfFile, CharmmParameterSet, CharmmPSFWarning
from simtk.openmm.app import LJPME, PME, HBonds
from simtk.openmm import Platform
import mdtraj as md

from fflip.omm import omm_vfswitch  # this should be from omm
from fflip.omm.playpara import (
    prepare_system, BrutalNonbondedParameter, change_nb_exceptions,
    change_drude_ff_parameters
)


# ---------------------------- File Handling ------------------------------

def copy_folder(sample, destination):
    des = os.path.abspath(destination)
    des_parent = '/'.join(des.split('/')[:-1])
    if not os.path.isdir(des_parent):
        os.system("mkdir -p " + des_parent)
    os.system("cp -r {} {}".format(sample, destination))


def copy_one_layer(sample, destination):
    os.system("cp {}/* {}/".format(sample, destination))


def check_and_make_dir(directory):
    if not os.path.isdir(directory):
        os.system("mkdir -p {}".format(directory))


# ---------------------------- Statistical & Thermo ------------------------

def beta_kjmol(temperature_kelvin):
    """Thermodynamic beta"""
    return 1.0/(temperature_kelvin * u.kelvin * u.MOLAR_GAS_CONSTANT_R).\
        value_in_unit(u.kilojoule_per_mole)


def beta_kcalmol(temperature_kelvin):
    return 1.0/(temperature_kelvin * u.kelvin * u.MOLAR_GAS_CONSTANT_R).\
        value_in_unit(u.kilocalories_per_mole)


# --------------------------- OpenMM Iteractions ---------------------------

def read_structure_parameter_files(psf, toppar):
    parameters = CharmmParameterSet(*toppar)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CharmmPSFWarning)
        psf = CharmmPsfFile(psf)
    topology = md.Topology.from_openmm(psf.topology)
    return psf, topology, parameters


def get_old_new_parameters(paragroups, paraoffsets):
    center_names = []
    par_types = []
    old_paras = []
    new_paras = []
    for g, offset in zip(paragroups, paraoffsets):
        # Attention: original_p is just float number, not u.Quantity
        center_names.append(g.center_names[0])
        par_types.append(g.par_type)
        old_paras.append(g.original_p)
        new_paras.append(g.original_p * (1 + offset))
    return center_names, par_types, old_paras, new_paras


def get_contexts_using_offset(
    system, topology, paragroups, paraoffsets, **kwargs
):
    prepare_system(system, topology, paragroups)
    contexts = []
    for g, offset in zip(paragroups, paraoffsets):
        if 'use_platform' in kwargs:
            platform = Platform.getPlatformByName(kwargs['use_platform'])
            context_g = Context(
                system, LangevinIntegrator(1., 1., 1.), platform
            )
        else:
            context_g = Context(
                system, LangevinIntegrator(1., 1., 1.)
            )
        if g.par_type:
            context_g.setParameter(g.force_names[0], offset)
            if g.par_type == "charge":
                for j, force_name in enumerate(g.force_names[1]):
                    context_g.setParameter(force_name, offset * g.roc[j])
                for k, force_name in enumerate(g.force_names[2]):
                    context_g.setParameter(force_name, offset * g.ron[k])
        contexts.append(context_g)
    return contexts


def get_contexts_using_offset_defult_platform(
    system, topology, paras_groups, paras_offsets
):
    """ use offset method """
    warnings.warn("This is a incomplete code, quitting ...")
    contexts = []
    for pgs, pos in zip(paras_groups, paras_offsets):
        sys = copy.deepcopy(system)
        prepare_system(sys, topology, pgs)
        context = Context(sys, LangevinIntegrator(1., 1., 1.))
        for g, offset in zip(pgs, pos):
            if g.par_type:
                context.setParameter(g.force_names[0], offset)
            if g.par_type == "charge":
                for j, force_name in enumerate(g.force_names[1]):
                    context.setParameter(force_name, offset * g.roc[j])
                for k, force_name in enumerate(g.force_names[2]):
                    context.setParameter(force_name, offset * g.ron[k])
        contexts.append(context)
    return contexts


def get_contexts(system_old, psf, paragroups, paraoffsets, **kwargs):
    """ use brutal method """
    contexts = []
    topology = md.Topology.from_openmm(psf.topology)
    for g, offset in zip(paragroups, paraoffsets):
        system = copy.deepcopy(system_old)
        if g.par_type:
            system, old_p, new_p = BrutalNonbondedParameter(
                system, topology, g, offset
            )
        change_nb_exceptions(psf, system)
        nonbonded_method = None
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                nonbonded_method = force.getNonbondedMethod()
                r_off = force.getCutoffDistance().value_in_unit(u.nanometer)
                # if PME
                if nonbonded_method == 4:
                    try:
                        r_on = force.getSwitchingDistance().value_in_unit(
                            u.nanometer
                        )
                        if r_on < 0.6 or r_off - r_on > 0.5:
                            r_on = r_off - 0.2
                    except:
                        r_on = r_off - 0.2
        if not nonbonded_method:
            pass  # TODO: it should raise some error ... (Yalun Yu)
        if nonbonded_method == 4:
            ''' Add CHARMM combination rule and switch functions '''
            system = omm_vfswitch.vfswitch(system, psf, r_on, r_off)
            ''' We don't want the default Long Range Correction, always '''
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setUseDispersionCorrection(False)
                if isinstance(force, CustomNonbondedForce):
                    force.setUseLongRangeCorrection(False)
        if 'use_platform' in kwargs:
            platform = Platform.getPlatformByName(kwargs['use_platform'])
            print('Using Platform {}'.format(kwargs['use_platform']))
            if kwargs['use_platform'] == 'CUDA':
                context_g = Context(
                    system, LangevinIntegrator(1., 1., 1.), platform,
                    {'Precision': 'mixed'}
                )
            else:
                context_g = Context(
                    system, LangevinIntegrator(1., 1., 1.), platform
                )
        else:
            context_g = Context(system, LangevinIntegrator(1., 1., 1.))
        contexts.append(context_g)
    return contexts


def get_contexts_default_platform(
        system_old, psf, paras_groups, paras_offsets
):
    contexts = []
    topology = md.Topology.from_openmm(psf.topology)
    for pgs, pos in zip(paras_groups, paras_offsets):
        system = copy.deepcopy(system_old)
        for g, offset in zip(pgs, pos):
            if g.par_type:
                system, old_p, new_p = BrutalNonbondedParameter(
                    system, topology, g, offset
                )
        change_nb_exceptions(psf, system)
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                nonbonded_method = force.getNonbondedMethod()
                r_off = force.getCutoffDistance().value_in_unit(u.nanometer)
                # if PME
                if nonbonded_method == 4:
                    try:
                        r_on = force.getSwitchingDistance().value_in_unit(
                            u.nanometer
                        )
                        if r_on < 0.6 or r_off - r_on > 0.5:
                            r_on = r_off - 0.2
                    except:
                        r_on = r_off - 0.2

        if not nonbonded_method:
            pass  # TODO: should raise some error, complete later ... (YYL)

        if nonbonded_method == 4:
            ''' Add CHARMM combination rule and switch functions '''
            system = omm_vfswitch.vfswitch(system, psf, r_on, r_off)
            ''' We don't want the default Long Range Correction, always '''
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setUseDispersionCorrection(False)
                if isinstance(force, CustomNonbondedForce):
                    force.setUseLongRangeCorrection(False)
        context_g = Context(system, LangevinIntegrator(1.0, 1.0, 1.0))
        contexts.append(context_g)

    return contexts


def create_drude_contexts(system_old, psf, paragroups, paraoffsets, **kwargs):
    contexts = []
    topology = md.Topology.from_openmm(psf.topology)
    for g, offset in zip(paragroups, paraoffsets):
        system = copy.deepcopy(system_old)
        if g.par_type:
            change_drude_ff_parameters(system, topology, g, offset)
        change_nb_exceptions(psf, system, isdrude=True)
        nonbonded_method = None
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                nonbonded_method = force.getNonbondedMethod()
                r_off = force.getCutoffDistance().value_in_unit(u.nanometer)
                # if PME
                if nonbonded_method == 4:
                    try:
                        r_on = force.getSwitchingDistance().value_in_unit(
                            u.nanometer
                        )
                        if r_on < 0.6 or r_off - r_on > 0.5:
                            r_on = r_off - 0.2
                    except:
                        r_on = r_off - 0.2
        if not nonbonded_method:
            pass  # TODO: it should raise some error ... (Yalun Yu)
        if nonbonded_method == 4:
            ''' Add CHARMM combination rule and switch functions '''
            system = omm_vfswitch.vfswitch(system, psf, r_on, r_off)
            ''' We don't want the default Long Range Correction, always '''
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setUseDispersionCorrection(False)
                if isinstance(force, CustomNonbondedForce):
                    force.setUseLongRangeCorrection(False)
        if 'use_platform' in kwargs:
            platform = Platform.getPlatformByName(kwargs['use_platform'])
            print('Using Platform {}'.format(kwargs['use_platform']))
            if kwargs['use_platform'] == 'CUDA':
                context_g = Context(
                    system, LangevinIntegrator(1., 1., 1.), platform,
                    {'Precision': 'mixed'}
                )
            else:
                context_g = Context(
                    system, LangevinIntegrator(1., 1., 1.), platform
                )
        else:
            context_g = Context(system, LangevinIntegrator(1., 1., 1.))
        contexts.append(context_g)
    return contexts


def setup_sim_system(psf, topology, parameters, paragroups, paraoffsets,
                     nonbonded_method, change_par=False, r_off=1.0, r_on=0.8):
    system = psf.createSystem(
        parameters,
        nonbondedMethod=nonbonded_method,
        constraints=HBonds,
        nonbondedCutoff=r_off * u.nanometer,
        ewaldErrorTolerance=0.0005
    )
    # add if necessary:
    # switchDistance=1.0 * u.nanometer
    if change_par:
        for g, offset in zip(paragroups, paraoffsets):
            # the old_para, new_para are dummy
            system, old_parameter, new_parameter = BrutalNonbondedParameter(
                system, topology, g, offset
            )

    ''' Add CHARMM combination rule and switch functions '''
    if nonbonded_method == PME:
        system = omm_vfswitch.vfswitch(system, psf, r_on, r_off)

    ''' We don't want the default Long Range Correction, always '''
    if nonbonded_method == PME:
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setUseDispersionCorrection(False)
            if isinstance(force, CustomNonbondedForce):
                force.setUseLongRangeCorrection(False)
    return system


# --------------------------- parsing inputs ----------------------------

def gen_omm_options_file(md_options, filename):
    """
    Write MD simulation inputs into a file
    Args:
        md_options: a dictionary contains simulation parameters
        filename: the file to write
    """
    with open(filename, 'w') as f:
        for item in md_options:
            f.write(item + '::' + str(md_options[item]))
            f.write('\n')


def get_md_options(filename):
    """
    Retrieve the MD simulation inputs from the .mdo file
    """
    md_options = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
    for l in lines:
        key = str(l.strip().split('::')[0])
        content = l.strip().split('::')[1]
        md_options[key] = content
    return md_options


def parse_md_options(md_options):
    """
    md_options are a dictionary containing input paramters
    """
    mdo = dict()
    for item in md_options:
        if item in [
            'psf', 'crd', 'lipid', 'ff', 'lipname', 'toppar_path',
            'change_param', 'sfile', 'integrator', 'barostat'
        ]:
            mdo[item] = str(md_options[item])
        elif item == 'surface_tension':
            # TODO: unify the naming
            mdo['surf_ts'] = 20 * float(md_options[item])
            mdo[item] = 20 * float(md_options[item])
        elif item in ['box', 'temperature',]:
            mdo[item] = float(md_options[item])
    if 'boxx' in md_options:
        mdo['xvec'] = float(md_options['boxx'])
        if 'boxy' in md_options:
            mdo['yvec']= float(md_options['boxy'])
        else:
            mdo['yvec'] = float(md_options['boxx'])
        if 'boxz' in md_options:
            mdo['zvec'] = float(md_options['boxz'])
        else:
            mdo['zvec'] = float(md_options['boxx'])
    if 'zmode' in md_options:
        mdo['zmode'] = int(md_options['zmode'])
    else:
        mdo['zmode'] = 0
    return mdo


def filter_solution(file_to_load='solution.txt', threshold=0.0001):
    data = np.loadtxt(file_to_load)
    sol = []
    for d in data:
        if np.abs(d) < threshold:
            sol.append(0)
        else:
            sol.append(d)
    return sol

def create_system_with_lj_offset(
    parameter_group, parameter_offset, psf_file, parameter_files,
    unitcell_lengths, nbmethod=LJPME, nbcutoff=1.0, change_14=True
):
    parameters = CharmmParameterSet(*parameter_files)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CharmmPSFWarning)
        psf = CharmmPsfFile(psf_file)
    topology = md.Topology.from_openmm(psf.topology)
    if parameter_group.par_type == 'sigma':
        par_type = 'rmin'
    elif parameter_group.par_type == 'epsilon':
        par_type = 'epsilon'
    else:
        raise Exception(
            "Parameter type {} not supported by this function".format(
                parameter_group.par_type
            )
        )
    # initialize atom_type
    atom_type = None
    for name in parameter_group.center_names:
        atoms = topology.select(
            "name {}".format(name)
             # changed because parameter_group.lipid.upper
             # is not necessarily resname
        )
        first_atom_index = int(atoms[0])
        assert psf.atom_list[first_atom_index].name == name
        if not atom_type:
            atom_type = psf.atom_list[first_atom_index].attype
        else:
            assert atom_type == psf.atom_list[first_atom_index].attype,\
                "Different atom types in one parameter (LJ) group!"
    if par_type == 'rmin':
        parameters.atom_types_str[atom_type].rmin = \
        parameters.atom_types_str[atom_type].rmin * (1 + parameter_offset)
        if change_14:
            parameters.atom_types_str[atom_type].rmin_14 = \
            parameters.atom_types_str[atom_type].rmin_14 * (1 + parameter_offset)
    elif par_type == 'epsilon':
        parameters.atom_types_str[atom_type].epsilon = \
        parameters.atom_types_str[atom_type].epsilon * (1 + parameter_offset)
        if change_14:
            parameters.atom_types_str[atom_type].epsilon_14 = \
            parameters.atom_types_str[atom_type].epsilon_14 * (1 + parameter_offset)
    psf.setBox(*unitcell_lengths)
    system_ = psf.createSystem(
        parameters,
        nonbondedMethod=nbmethod,
        constraints=HBonds,
        nonbondedCutoff=nbcutoff * u.nanometer,
        ewaldErrorTolerance=0.0001
    )
    return system_

# Function to change the LJ parameters compatible with NBFIX
def change_lj_param(psfworkflow,solution_file,lipid,cutoff_dist=10.0,change_14=True):
    sol = filter_solution(solution_file)
    parameter_sets = lipid.parse_groups()
    all_offsets = [gen_param_offset(ps, amount=sol[i]) \
    for i, ps in enumerate(parameter_sets)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CharmmPSFWarning)
        psf = CharmmPsfFile(psf_file)
    topology = md.Topology.from_openmm(psfworkflow.psf.topology)
    for g, offset in zip(parameter_sets, all_offsets):
        if g.par_type == 'sigma':
            par_type = 'rmin'
        elif g.par_type == 'epsilon':
            par_type = 'epsilon'
        atom_type = None
        for name in g.center_names:
            atoms = topology.select("name {}".format(name))
            first_atom_index = int(atoms[0])
            assert psf.atom_list[first_atom_index].name == name
            if not atom_type:
                atom_type = psf.atom_list[first_atom_index].atom_type
            else:
                assert atom_type == psf.atom_list[first_atom_index].attype, \
                "Different atom types in one parameter (LJ) group!"
        if par_type == 'rmin':
            psfworkflow.parameters.atom_types_str[atom_type].rmin *= \
            ( 1 + offset )
            if change_14:
                psfworkflow.parameters.atom_types_str[atom_type].rmin_14 *= \
                ( 1 + offset )
        elif par_type == 'epsilon':
            psfworkflow.parameters.atom_types_str[atom_type].epsilon *= \
            ( 1 + offset )
            if change_14:
                psfworkflow.parameters.atom_types_str[atom_type].epsilon_14 *= \
                ( 1 + offset )
    psfworkflow.create_system(cutoff_distance=cutoff_dist*u.angstrom)


# Function to change the charge parameters of the atoms
def change_charge_param(solution_file, lipid):
    sol = filter_solution(solution_file)
    parameter_sets = lipid.parse_groups()
    all_offsets = [gen_param_offset(ps, amount=sol[i]) \
    for i, ps in enumerate(parameter_sets)]
    for force in workflow.system.getForces():
        if isinstance(force, NonbondedForce):
            for g, offset in zip(parameter_sets, all_offsets):
                if g.par_type == 'charge':
                    # Initialize the center and neighbor atom lists.
                    atom_sele_neighb = []
                    atom_sele_center = []
                    for i in range(len(g.center_names)):
                        atom_sele_center.append(
                            topology.select("name {}".format(g.center_names[i]))
                        )
                    for i in range(len(g.neighbors)):
                        atom_sele_neighb.append(
                            topology.select("name {}".format(g.neighbors[i]))
                        )

                    # Change parameters for the center atoms.
                    for i in range(len(g.center_names)):
                        for atom in atom_sele_center[i]:
                            atom = int(atom)
                            charge, sigma, epsilon - force.getParticleParameters(atom)
                            charge_new = charge + u.Quantity(offset, unit=u.elementary_charge)
                            force.setParticleParameters(atom, charge_new, sigma, epsilon)

                    # Change parameters for the neighboring atoms.
                    for i in range(len(g.neighbors)):
                        for atom in atom_sele_neighb[i]:
                            atom = int(atom)
                            charge, sigma, epsilon = force.getParticleParameters(atom)
                            charge_new = charge + u.Quantity(offset*g.ron[i], unit=u.elementary_charge)
                            force.setParticleParameters(atom, charge_new, sigma, epsilon)
