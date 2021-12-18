#!/usr/bin/env python

from __future__ import division, print_function

import warnings
import numpy as np
import simtk.unit as u
from simtk.openmm import NonbondedForce, DrudeForce, CustomBondForce
import simtk.openmm as mm

from fflip.omm.exceptions import ParameterTypeNotExistentError,\
    GlobalParameterForceTypeError


def prepare_system(system, topology, groups):
    # Add(initialize) all the GlobalParameters
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            initialize_parameter_offset(topology, force, groups)


def add_global_parameter(force, name, value):
    """
    Args:
        force: 
        name: name of the parameter, name should indicates 
        nonbonded force subtype (charge, epsilon or sigma)
        value: default value
    Returns:
    """
    if not isinstance(force, NonbondedForce):
        raise GlobalParameterForceTypeError(force)

    index = force.addGlobalParameter(name, value)
    return index, name


def add_nonbonded_offset(force, name, particles, par_type):
    """
    This is assigning unit value, we might want something different ...
    Args:
        force: OpenMM forces, must be nonbonded force
        name: name of the force
        particles: int, particle (atom) index
        par_type: charge / epsilon / sigma

    Returns:
    """
    if not isinstance(force, NonbondedForce):
        raise GlobalParameterForceTypeError(force)
    if par_type == "epsilon":
        # name of the force, index of the particle
        for particle in particles:
            q, sig, eps = force.getParticleParameters(particle)
            # force.addParticleParameterOffset(name, particle, 0, 0, 1)
            force.addParticleParameterOffset(name, particle, 0, 0, eps)
    elif par_type == "sigma":
        for particle in particles:
            q, sig, eps = force.getParticleParameters(particle)
            # force.addParticleParameterOffset(name, particle, 0, 1, 0)
            force.addParticleParameterOffset(name, particle, 0, sig, 0)
    elif par_type == "charge":
        # Note: this is called by the initializer,
        # which takes care of subgroups of charge,
        # no worry here.
        for particle in particles:
            # q, sig, eps = force.getParticleParameters(particle)
            force.addParticleParameterOffset(name, particle, 1, 0, 0)
            # force.addParticleParameterOffset(name, particle, q, 0, 0)
    else:
        raise Exception("???")


def smart_select(topology, group_info):
    """

    :param topology:
    :param group_info:
    :return:
    """
    par_type = group_info.par_type
    if par_type == "charge":
        atom_sele_cooper = []
        atom_sele_neighb = []
        atom_sele_center = []
        for name in group_info.center_names:
            # several arrays like [array([...]), array([...])]
            atom_sele_center.append(topology.select(
                "name {}".format(name)))
        for name in group_info.cooperators:
            atom_sele_cooper.append(topology.select(
                "name {}".format(name)))
        for name in group_info.neighbors:
            atom_sele_neighb.append(topology.select(
                "name {}".format(name)))
        # these element themselves are list containing array(s)
        return [atom_sele_center, atom_sele_cooper, atom_sele_neighb]

    if par_type == "sigma" or par_type == "epsilon":
        atom_sele = []
        if group_info.center_names[0] == "all":
            sele = topology.select_atom_indices(selection='all')
            atom_sele.append(sele)
        elif group_info.center_names[0] == "no water":
            sele = topology.select("not water")
            atom_sele.append(sele)
        elif group_info.center_names[0] == "headgroup":
            sele = topology.select(
                "not ((name =~ 'C2[3-9]') or (name =~ 'C3[3-9]') \
                or (name =~ 'C21[0-9]') or (name =~ 'C31[0-9]') \
                or (name =~ 'H[3-9][R-T]') or (name =~ 'H[3-9][X-Z]') \
                or (name =~ 'H1[0-9][R-T]') or (name =~ 'H1[0-9][X-Z]') \
                or water)")
            atom_sele.append(sele)
        else:
            for name in group_info.center_names:
                atom_sele.append(topology.select("name {}".format(name)))
        return [atom_sele]


def convert_selection_to_int(atom_sele):
    # Hope this is what we thought...
    sele_int = []
    for item in atom_sele:
        arr_int = []
        for arr in item:
            sele = []
            for atom in arr:
                sele.append(int(atom))
            arr_int.append(sele)
        sele_int.append(arr_int)
    return sele_int  # [arr[sele1,sele2], array_of_ints_2, array_of_ints_3]


def reassign_global_parameter(context, names, values):
    """
    This alone is not very useful, should be used with
    the actual objective function.
    Args:
        context:
        names:
        values:

    Returns:
    """
    for name, value in zip(names, values):
        context.setParameter(name, value)


def initialize_parameter_offset(topology, force, groups):
    """
    :param topology:
    :param force:
    :param groups:
    :return:
    """
    if not isinstance(force, NonbondedForce):
        raise  GlobalParameterForceTypeError(force)

    for group in groups:
        if group.par_type:
            # add Global Parameter
            add_global_parameter(force, group.force_names[0], 0.0)
            if group.par_type is "charge":
                for name_list in group.force_names[1:]:
                    for name in name_list:
                        # for now, just say the default for
                        # all the global parameter is 0.0, actually, it gets
                        # forgotten once the optimizer is called.
                        add_global_parameter(force, name, 0.0)

            atom_sele = smart_select(topology, group)
            if group.par_type == "charge":
                center, cooperator, neighbor = convert_selection_to_int(
                    atom_sele
                )
                for p_group in center:
                    # uniform force_name for center atoms
                    add_nonbonded_offset(
                        force, group.force_names[0], p_group, group.par_type
                    )
                for p_group, f_name in zip(cooperator, group.force_names[1]):
                    # p_group is group of particle indexes.
                    add_nonbonded_offset(force, f_name, p_group, group.par_type)
                for p_group, f_name in zip(neighbor, group.force_names[2]):
                    add_nonbonded_offset(force, f_name, p_group, group.par_type)
            else:
                selection = convert_selection_to_int(atom_sele)
                # Double check here !!!
                for p_group in selection[0]:
                    add_nonbonded_offset(
                        force, group.force_names[0], p_group, group.par_type
                    )


def BrutalNonbondedParameter(system, topology, paragroup, paraoffset):
    return brutal_nonbonded_parameter(system, topology, paragroup, paraoffset)


def brutal_nonbonded_parameter(system, topology, paragroup, paraoffset):
    """
    Change some nonbonded parameters of the system.
    """
    par_type = paragroup.par_type
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            # CHARGE
            if par_type == 'charge':
                atom_sele_cooper = []
                atom_sele_neighb = []
                atom_sele_center = []
                # Fill in atom selections for center atoms and others
                for i in range(len(paragroup.center_names)):
                    atom_sele_center.append(
                        topology.select(
                            "name {}".format(paragroup.center_names[i])
                        )
                    )
                for i in range(len(paragroup.cooperators)):
                    atom_sele_cooper.append(
                        topology.select(
                            "name {}".format(paragroup.cooperators[i])
                        )
                    )
                for i in range(len(paragroup.neighbors)):
                    atom_sele_neighb.append(
                        topology.select(
                            "name {}".format(paragroup.neighbors[i])
                        )
                    )
                # Change charge parameters for the center atoms
                for i in range(len(paragroup.center_names)):
                    for atom in atom_sele_center[i]:
                        atom = int(atom)
                        charge, sigma, epsilon = \
                            force.getParticleParameters(atom)
                        charge_new = charge + u.Quantity(
                            paraoffset, unit=u.elementary_charge
                        )
                        force.setParticleParameters(
                            atom, charge_new, sigma, epsilon
                        )
                # Sotre the charge information in 'old' and 'new'
                # respectively after loop
                old = charge._value
                new = charge_new._value

                ''' Change charge parameters for other associated atoms '''
                for i in range(len(paragroup.cooperators)):
                    for atom in atom_sele_cooper[i]:
                        atom = int(atom)
                        charge, sigma, epsilon = \
                            force.getParticleParameters(atom)
                        charge_new = charge + u.Quantity(
                            paraoffset*paragroup.roc[i],
                            unit=u.elementary_charge
                        )
                        force.setParticleParameters(
                            atom, charge_new, sigma, epsilon
                        )
                for i in range(len(paragroup.neighbors)):
                    for atom in atom_sele_neighb[i]:
                        atom = int(atom)
                        charge, sigma, epsilon = \
                            force.getParticleParameters(atom)
                        charge_new = charge + u.Quantity(
                            paraoffset*paragroup.ron[i],
                            unit=u.elementary_charge
                        )
                        force.setParticleParameters(
                            atom, charge_new, sigma, epsilon
                        )
            # SIGMA
            elif par_type == 'sigma':
                atom_sele = []
                # If change for all
                if paragroup.center_names[0] == 'all':
                    atom_sele = topology.select_atom_indices(selection='all')
                    for atom in atom_sele:
                        atom = int(atom)
                        charge, sigma, epsilon = \
                            force.getParticleParameters(atom)
                        # in this case, the paraoffset should be fraction
                        sigma_new = sigma * (1 + paraoffset)
                        force.setParticleParameters(
                            atom, charge, sigma_new, epsilon
                        )
                    old = 1
                    new = 1 + paraoffset
                # If change for all non water atoms
                elif paragroup.center_names[0] == 'no water':
                    atom_not_sele = topology.select('water')
                    atom_all = topology.select_atom_indices(selection='all')
                    for atom in atom_all:
                        if atom not in atom_not_sele:
                            charge, sigma, epsilon = \
                                force.getParticleParameters(atom)
                            sigma_new = sigma*(1 + paraoffset)
                            force.setParticleParameters(
                                atom, charge, sigma_new, epsilon
                            )
                    old = 1
                    new = 1 + paraoffset
                # Most generic situation
                else:
                    for sister in range(len(paragroup.center_names)):
                        atom_sele.append(
                            topology.select(
                                "name {}".format(
                                    paragroup.center_names[sister])
                            )
                        )
                        for atom in atom_sele[sister]:
                            atom = int(atom)
                            charge, sigma, epsilon = \
                                force.getParticleParameters(atom)
                            sigma_new = sigma * (1 + paraoffset)
                            force.setParticleParameters(
                                atom, charge, sigma_new, epsilon
                            )
                    old = sigma._value
                    new = sigma_new._value
            # EPSILON
            elif par_type == 'epsilon':
                atom_sele = []
                # If change for all
                if paragroup.center_names[0] == 'all':
                    atom_sele = topology.select_atom_indices(selection='all')
                    for atom in atom_sele:
                        atom = int(atom)
                        charge, sigma, epsilon = \
                            force.getParticleParameters(atom)
                        epsilon_new = epsilon * (1 + paraoffset)
                        force.setParticleParameters(
                            atom, charge, sigma, epsilon_new
                        )
                    old = 1
                    new = 1 + paraoffset
                # If change for all non water atoms
                elif paragroup.center_names[0] == 'no water':
                    atom_not_sele = topology.select('water')
                    atom_all = topology.select_atom_indices(selection='all')
                    for atom in atom_all:
                        if atom not in atom_not_sele:
                            atom = int(atom)
                            charge, sigma, epsilon = \
                                force.getParticleParameters(atom)
                            epsilon_new = epsilon * (1 + paraoffset)
                            force.setParticleParameters(
                                atom, charge, sigma, epsilon_new
                            )
                    old = 1
                    new = 1 + paraoffset
                # Most generic situation
                else:
                    for sister in range(len(paragroup.center_names)):
                        atom_sele.append(topology.select(
                            "name {}".format(paragroup.center_names[sister])))
                        for atom in atom_sele[sister]:
                            atom = int(atom)
                            charge, sigma, epsilon = \
                                force.getParticleParameters(atom)
                            epsilon_new = epsilon * (1 + paraoffset)
                            force.setParticleParameters(
                                atom, charge, sigma, epsilon_new
                            )
                    # sotre information about epsilon
                    old = epsilon._value
                    new = epsilon_new._value
            # Wrong parameter type
            else:
                raise ParameterTypeNotExistentError(par_type)
    return system, old, new


def change_nb_exceptions(psf, system, isdrude=False):
    force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            break
    assert force is not None
    from math import sqrt
    sigma_scale = 2**(-1/6)
    ene_conv = u.kilocalorie_per_mole.conversion_factor_to(u.kilojoule_per_mole)
    length_conv = u.angstrom.conversion_factor_to(u.nanometer)
    num_exceptions = force.getNumExceptions()
    for i in range(num_exceptions):
        index1, index2, cprod, _, _ = force.getExceptionParameters(i)
        if cprod._value == 0:
            continue
        particle1_parameters = force.getParticleParameters(index1)
        particle2_parameters = force.getParticleParameters(index2)
        chg1 = particle1_parameters[0]._value
        sig1 = particle1_parameters[1]._value
        eps1 = particle1_parameters[2]._value
        chg2 = particle2_parameters[0]._value
        sig2 = particle2_parameters[1]._value
        eps2 = particle2_parameters[2]._value
        rmin_14_1 = psf.atom_list[index1].type.rmin_14
        rmin_1_org = psf.atom_list[index1].type.rmin
        rmin_14_2 = psf.atom_list[index2].type.rmin_14
        rmin_2_org = psf.atom_list[index2].type.rmin
        epsilon_14_1 = psf.atom_list[index1].type.epsilon_14
        epsilon_1_org = psf.atom_list[index1].type.epsilon
        epsilon_14_2 = psf.atom_list[index2].type.epsilon_14
        epsilon_2_org = psf.atom_list[index2].type.epsilon
        # charge
        charge_prod = (chg1 * chg2)
        # epsilon
        if epsilon_14_1 != epsilon_1_org or isdrude:
            eps_1 = abs(epsilon_14_1) * ene_conv
        else:
            eps_1 = eps1
        if epsilon_14_2 != epsilon_2_org or isdrude:
            eps_2 = abs(epsilon_14_2) * ene_conv
        else:
            eps_2 = eps2
        epsilon = (sqrt(eps_1 * eps_2))
        # sigma
        if rmin_14_1 != rmin_1_org or isdrude:
            sig_1 = rmin_14_1 * (length_conv * sigma_scale)
        else:
            sig_1 = sig1 / 2
        if rmin_14_2 != rmin_2_org or isdrude:
            sig_2 = rmin_14_2 * (length_conv * sigma_scale)
        else:
            sig_2 = sig2 / 2
        sigma = sig_1 + sig_2
        force.setExceptionParameters(
            i, index1, index2, charge_prod, sigma, epsilon
        )


def find_nb_parameter(system, topology, names, par_types):
    values = []
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            for name, par_type in zip(names, par_types):
                atoms = topology.select("name {}".format(name))
                atom = int(atoms[0])
                charge, sigma, epsilon = force.getParticleParameters(atom)
                if par_type == 'charge':
                    values.append(charge._value)
                elif par_type == 'sigma':
                    values.append(sigma._value)
                elif par_type == 'epsilon':
                    values.append(epsilon._value)
    return values


def find_drude_nb_forces(system):
    # get the forces we need
    drudeforce = None
    nbforce = None
    for force in system.getForces():
        if isinstance(force, DrudeForce):
            drudeforce = force
        if isinstance(force, NonbondedForce):
            nbforce = force
    return drudeforce, nbforce


def change_polarizability(drudeforce, nbforce, drude_id, change_of_alpha):
    parameters = drudeforce.getParticleParameters(drude_id)
    # print(parameters)
    drude_particle_id = parameters[0]
    attached_id = parameters[1]
    qd_value = parameters[5]._value
    alpha_value = parameters[6]._value
    new_alpha_value = alpha_value * (1 + change_of_alpha)
    new_qd_value = np.round(
        -np.sqrt((qd_value**2/alpha_value)*new_alpha_value), 4
    )
    change_of_qd = new_qd_value - qd_value
    chg, sig, eps = nbforce.getParticleParameters(attached_id)
    dchg, dsig, deps = nbforce.getParticleParameters(drude_particle_id)
    new_chg = chg + u.Quantity(-change_of_qd, unit=u.elementary_charge)
    new_dchg = dchg + u.Quantity(change_of_qd, unit=u.elementary_charge)
    parameters[5]._value = new_qd_value
    parameters[6]._value = new_alpha_value
    drudeforce.setParticleParameters(drude_id, *parameters)
    nbforce.setParticleParameters(attached_id, new_chg, sig, eps)
    nbforce.setParticleParameters(drude_particle_id, new_dchg, dsig, deps)


def change_thole(drudeforce, thole_id, change_of_thole):
    parameters = drudeforce.getScreenedPairParameters(thole_id)
    particle_1 = parameters[0]
    particle_2 = parameters[1]
    thole = parameters[2] + change_of_thole
    drudeforce.setScreenedPairParameters(thole_id, particle_1, particle_2, thole)


def change_charge(nbforce, particle_id, change_of_charge):
    charge, sigma, epsilon = nbforce.getParticleParameters(particle_id)
    new_charge = charge + u.Quantity(change_of_charge, unit=u.elementary_charge)
    nbforce.setParticleParameters(particle_id, new_charge, sigma, epsilon)


def reassign_charge_by_id(nbforce, particle_id, new_charge):
    charge, sigma, epsilon = nbforce.getParticleParameters(particle_id)
    new_c = u.Quantity(new_charge, unit=u.elementary_charge)
    nbforce.setParticleParameters(particle_id, new_c, sigma, epsilon)


def change_sigma(nbforce, particle_id, change_of_sigma):
    """
    change_of_sigma must be portion of the original
    """
    charge, sigma, epsilon = nbforce.getParticleParameters(particle_id)
    new_sigma = sigma * (1 + change_of_sigma)
    nbforce.setParticleParameters(particle_id, charge, new_sigma, epsilon)


def reassign_sigma_by_id(nbforce, particle_id, new_sigma):
    """
    change_of_sigma must be portion of the original
    """
    charge, sigma, epsilon = nbforce.getParticleParameters(particle_id)
    new_sig = u.Quantity(new_sigma, unit=u.nanometer)
    nbforce.setParticleParameters(particle_id, charge, new_sig, epsilon)


def change_epsilon(nbforce, particle_id, change_of_epsilon):
    """
    change_of_epsilon must be portion of the original
    """
    charge, sigma, epsilon = nbforce.getParticleParameters(particle_id)
    new_epsilon = epsilon * (1 + change_of_epsilon)
    nbforce.setParticleParameters(particle_id, charge, sigma, new_epsilon)


def reassign_epsilon_by_id(nbforce, particle_id, new_epsilon):
    """
    change_of_epsilon must be portion of the original
    """
    charge, sigma, epsilon = nbforce.getParticleParameters(particle_id)
    new_eps = u.Quantity(new_epsilon, unit=u.kilojoule_per_mole)
    nbforce.setParticleParameters(particle_id, charge, sigma, new_eps)


def change_drude_ff_parameters(
    system, topology, parameter, paraoffset, psf=None
):
    """
    Change some nonbonded parameters of the system.
    """
    par_type = parameter.par_type
    drudeforce, nbforce = find_drude_nb_forces(system)
    particleMap = {}
    particleMapR = {}
    for i in range(drudeforce.getNumParticles()):
        particleMap[drudeforce.getParticleParameters(i)[0]] = i
        particleMapR[i] = drudeforce.getParticleParameters(i)[0]
    # custombondforce = None
    # for force in system.getForces():
    #     if isinstance(force, CustomBondForce):
    #         custombonforce = force
    #         break
    # nbt14force_new = mm.CustomBondForce(
    #     '-138.935456*charge_prod*(1.0+0.5*screen*r)*exp(-1.0*screen*r)/r'
    # )
    # nbt14force_new.addPerBondParameter("charge_prod")
    # nbt14force_new.addPerBondParameter("screen")
    # nbt14force_new.setForceGroup(psf.NONBONDED_FORCE_GROUP)
    # # new nbthole (currently changing existing NBTHOLE is not supported!)
    if par_type is 'nbthole':
        pass
    #     num_nbt14_new = 0
    #     atype1 = parameter.center_names[0]
    #     atype2 = parameter.center_names[1]
    #     for dih in psf.dihedral_list:
    #         a1, a4 = dih.atom1, dih.atom4
    #         if sorted([a1.attype, a4.attype]) == sorted([atype1, atype2]):
    #             # These (d1_idx, d4_idx) are for the drudeforce only
    #             d1_idx = particleMap[a1.idx + 1]
    #             d4_idx = particleMap[a4.idx + 1]
    #             d1_parameters = drudeforce.getParticleParameters(d1_idx)
    #             d4_parameters = drudeforce.getParticleParameters(d4_idx)
    #             # q1, q4 = d1_parameters[5]._value, d4_parameters[5]
    #             # # 1000 is for unit conversion
    #             # alpha1 = pow(1000 * d1_parameters[6]._value, -1./6.)
    #             # alpha4 = pow(1000 * d4_parameters[6]._value, -1./6.)
    #             # screen = parameter.original_p + paraoffset
    #             # screen = screen * alpha1 * alpha4 * 10.0
    #             # q_prod = q1 * q4
    #             # nbt14force_new.addBond(d1_idx, d4_idx, [q_prod, screen])
    #             q1, q4 = a1.charge, a4.charge
    #             # 1000 is for unit conversion
    #             alpha1 = pow(-1 * psf.drudeconsts_list[a1.idx][0],-1./6.)
    #                      # pow(1000 * d1_parameters[6]._value, -1./6.)
    #             alpha4 = pow(-1 * psf.drudeconsts_list[a4.idx][0],-1./6.)
    #                      # pow(1000 * d4_parameters[6]._value, -1./6.)
    #             nbt_value = parameter.original_p + paraoffset
    #             screen = nbt_value * alpha1 * alpha4 * 10
    #             q_prod = q1 * q4
    #             nbt14force_new.addBond(
    #                 a1.idx, a4.idx, [q_prod, screen])
    #             num_nbt14_new += 1
    #     if num_nbt14_new > 0:
    #         system.addForce(nbt14force_new)

    if par_type is 'alpha':
        particles = []
        for i, dp in enumerate(parameter.drude_particles):
            particles += list(
                topology.select(
                    "resname {} and name {}".format(
                        parameter.lipid.upper(), dp
                    )
                )
            )
        drude_indexes = []
        for index_drude in range(drudeforce.getNumParticles()):
            parameters = drudeforce.getParticleParameters(index_drude)
            if parameters[0] in particles:
                drude_indexes.append(index_drude)
        if not drude_indexes:
            warnings.warn("No Drude Particle Selected!")
        for di in drude_indexes:
            change_polarizability(drudeforce, nbforce, di, paraoffset)
    if par_type is 'thole':
        particles = []
        for i, dp in enumerate(parameter.drude_particles):
            particles += list(
                topology.select(
                    "resname {} and name {}".format(
                        parameter.lipid.upper(), dp
                    )
                )
            )
        d_particles = []
        for index_d in range(drudeforce.getNumParticles()):
            dd_param = drudeforce.getParticleParameters(index_d)
            if dd_param[0] in particles:
                d_particles.append(index_d)
        thole_indexes = []
        for index_t in range(drudeforce.getNumScreenedPairs()):
            sp_param = drudeforce.getScreenedPairParameters(index_t)
            if sp_param[0] in d_particles or sp_param[1] in d_particles:
                # name1 = psf.atom_list[particleMapR[sp_param[0]]].name
                # name2 = psf.atom_list[particleMapR[sp_param[1]]].name
                # if 'DC2' in [name1, name2] or ('DC3' in [name1, name2] and 'DOP2' in [name1, name2]):
                #     print(name1, name2)
                thole_indexes.append(index_t)
                if sp_param[0] in d_particles and sp_param[1] in d_particles:
                    # repeat to include both changes
                    thole_indexes.append(index_t)
        if not thole_indexes:
            warnings.warn("No Thole Pairs Selected!")
        for ti in thole_indexes:
            change_thole(drudeforce, ti, paraoffset)
    elif par_type is 'charge':
        particles = []
        # Fill in atom selections for center atoms and neighbors
        for i, cn in enumerate(parameter.center_names):
            particles += list(
                topology.select(
                    "resname {} and name {}".format(
                        parameter.lipid.upper(), cn
                    )
                )
            )
        for pid in particles:
            change_charge(nbforce, int(pid), paraoffset)
        for i, nb in enumerate(parameter.neighbors):
            neighbors = list(
                topology.select(
                    "resname {} and name {}".format(
                        parameter.lipid.upper(), nb
                    )
                )
            )
            for nid in neighbors:
                change_charge(nbforce, int(nid), paraoffset * parameter.ron[i])
    elif par_type is 'sigma' or par_type is 'epsilon':
        warnings.warn(
            "Might not be effective if there is any NBFIX in the system!")
        particles = []
        for i, pn in enumerate(parameter.center_names):
            particles += list(
                topology.select(
                    "resname {} and name {}".format(
                        parameter.lipid.upper(), pn
                    )
                )
            )
        for particle in particles:
            if par_type is 'sigma':
                change_sigma(nbforce, int(particle), paraoffset)
            else:
                change_epsilon(nbforce, int(particle), paraoffset)


def find_drude_force_parameter(system, topology, name, plrzb=True, dname=None):
    """
    dname: name of the drude particle, if not provided, 
    the program will guess based on the name of the heavy atom
    """
    force = None
    for force in system.getForces():
        if isinstance(force, DrudeForce):
            ddforce = force
            break
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nbforce = force
            break
    assert isinstance(ddforce, DrudeForce)
    assert isinstance(nbforce, NonbondedForce)
    atoms = topology.select("name {}".format(name))
    atom = int(atoms[0])
    assert atom is not None
    if plrzb:
        if dname is None:
            dname = 'D' + name
        # datom is the indexing of the Drude Particle
        # in the full partical space
        datoms = topology.select("name {}".format(dname))
        datom = int(datoms[0])
        assert datom is not None
        for dindex in range(ddforce.getNumParticles()):
            dd_param = ddforce.getParticleParameters(dindex)
            if dd_param[0] == datom:
                break
    # charge
    nbparameters = nbforce.getParticleParameters(atom)
    charge = round(nbparameters[0]._value, 4)
    # alpha
    if plrzb:
        for index_drude in range(ddforce.getNumParticles()):
            parameters = ddforce.getParticleParameters(index_drude)
            if parameters[0] == datom:
                alpha = round(parameters[6]._value * 1000, 3)
    else:
        alpha = 0
    # thole - since I can only get pair values, this is lengthy
    if plrzb:
        neighbor_list = []
        for ti in range(ddforce.getNumScreenedPairs()):
            sp_param = ddforce.getScreenedPairParameters(ti)
            if sp_param[0] == dindex:
                neighbor_list.append(sp_param[1])
            if sp_param[1] == dindex:
                neighbor_list.append(sp_param[0])
        for ti in range(ddforce.getNumScreenedPairs()):
            sp_param = ddforce.getScreenedPairParameters(ti)
            if sp_param[0] in neighbor_list and sp_param[1] in neighbor_list:
                neighbors = [sp_param[0], sp_param[1]]
                neighbor_sum = sp_param[2]
                break
        bigsum = 0
        for ti in range(ddforce.getNumScreenedPairs()):
            sp_param = ddforce.getScreenedPairParameters(ti)
            if sp_param[0] == dindex and sp_param[1] in neighbors:
                bigsum += sp_param[2]
            if sp_param[1] == dindex and sp_param[0] in neighbors:
                bigsum += sp_param[2]
        thole = round((bigsum - neighbor_sum)/2, 3)
    else:
        thole = 0
    return charge, alpha, thole
