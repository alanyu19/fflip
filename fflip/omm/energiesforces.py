# -*- coding: utf-8 -*-

"""Utility and helper functions for OpenMM"""

from __future__ import absolute_import, division, print_function

from mdtraj.utils import lengths_and_angles_to_box_vectors
from fflip.omm.util import *
from fflip.omm.paragroup import *


class EnergyDecomposition(object):
    """Energy decompostion for CHARMM force field
    """
    def __init__(self, system):
        self.system = system
        self.context = Context(system, LangevinIntegrator(1., 1., 1.))

    def get_force_group(self, forces):
        self.forcegroups = {}
        for i in range(self.system.getNumForces()):
            force = self.system.getForce(i)
            force.setForceGroup(i)
            # make the force object to str and human readable
            f = str(force).split("::")[1].split("*")[0].strip()
            if f in forces or forces == "all":
                self.forcegroups[f] = i
        for i in range(self.system.getNumForces()):
            force = self.system.getForce(i)
            if isinstance(force, NonbondedForce):
                force.setReciprocalSpaceForceGroup(i+1)
                self.forcegroups["Reciprocal"] = i

    def energy_decomposition(self, energies):
        for f, i in self.forcegroups.items():
            energies[f].append(
                self.context.getState(
                    getEnergy=True, groups={i}
                ).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
            )

    def __call__(self, traj, **kwargs):
        """
        Args:
        - traj: the trajectory to get energies from
        - force_to_return (kwargs):
        HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce, CustomTorsionForce,
        CMAPTorsionForce, NonbondedForce, CMMotionRemover, and MORE TO COME ...
        - n_frames (kwargs): int, only the first n frames of the trajectory
        are used to calculate energies if this is provided
        """
        if "n_frames" not in kwargs:
            n_frames = traj.n_frames
        else:
            n_frames = kwargs["n_frames"]
        if "forces_to_return" not in kwargs:
            forces_to_return = "all"
        else:
            forces_to_return = kwargs["forces_to_return"]
        self.get_force_group(forces_to_return)
        energies = {}
        for f, i in self.forcegroups.items():
            # Initialize energy dictionary
            energies[f] = []
        for frame in range(n_frames):
            box_vectors = lengths_and_angles_to_box_vectors(
                *traj.unitcell_lengths[frame], *traj.unitcell_angles[frame]
            )
            self.context.setPeriodicBoxVectors(*box_vectors)
            self.context.setPositions(traj.xyz[frame])
            self.energy_decomposition(energies)
        energy_list = []
        for f, e in energies.items():
            energy_list.append(e)
        energy_array = np.swapaxes(np.array(energy_list), 0, 1)
        if energy_array.shape[1] == 1:
            return np.reshape(energy_array, -1)
        else:
            return energy_array


class ParameterEnergy(object):
    """Evaluate potential energies of trajectory frames in OpenMM.
    """
    def __init__(self, system, psf, paragroups=[], paraoffsets=[], **kwargs):
        """
        Args:
        - system (openmm.System): defines all forces in OpenMM
        - psf (openmm.app.CharmmPsfFile): OpenMM representation of a psf file
        - paragroups (list of nbgroup): Instances that define changed nonbonded parameters.
          If paragroups=[] (the default), the energies are evaluated for the original
          force field parameters.
        - paraoffsets(list of float): Offsets for nonbonded parameters (see OpenMM documentation
          of NonbondedForce)
        """
        topology = md.Topology.from_openmm(psf.topology)
        self.psf = psf
        self.system = system
        self.topology = topology
        self.name = "Potential Energy (kJ/mole)"
        self.paragroups = paragroups
        if self.paragroups:
            assert paraoffsets
            self.paraoffsets = paraoffsets
        else:
            # 'empty_nbgroup' is defined in paragroup
            self.paragroups = [empty_nbgroup()]
            self.paraoffsets = [0]
        if "use_new_method" in kwargs:
            use_new_method = kwargs["use_new_method"]
        else:
            use_new_method = False
        if use_new_method:
            self.contexts = get_contexts_using_offset(
                self.system, self.topology,
                self.paragroups, self.paraoffsets, **kwargs
            )
        else:
            self.contexts = get_contexts(
                self.system, self.psf, self.paragroups,
                self.paraoffsets, **kwargs
            )

    def __call__(self, traj, **kwargs):
        """
        Args:
        - traj (mdtraj.Trajectory): the trajectory object to get energies from
        - n_frames (kwargs): the first n frames of the trajectory to calculate energies
        """
        all_energies = []

        # create empty list for energies data
        for _ in range(len(self.paragroups)):
            all_energies.append([])
        if "n_frames" not in kwargs:
            n_frames = traj.n_frames
        else:
            n_frames = kwargs["n_frames"]
        for frame in range(n_frames):
            box_vectors = lengths_and_angles_to_box_vectors(
                *traj.unitcell_lengths[frame], *traj.unitcell_angles[frame]
            )
            for g in range(len(self.paragroups)):
                self.contexts[g].setPeriodicBoxVectors(*box_vectors)
                self.contexts[g].setPositions(traj.xyz[frame])
                all_energies[g].append(
                    self.contexts[g].getState(
                        getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                )
        # attention: the energies are returned, but the parameter change is
        # not recorded here
        if len(all_energies) == 1:
            return np.reshape(np.array(all_energies), (-1, 1))
        else:
            energy_array = np.swapaxes(np.array(all_energies), 0, 1)
            return energy_array


class DrudeEnergy(object):
    """Evaluate potential energies of trajectory frames in OpenMM.
    """
    def __init__(self, system, psf, paragroups=[], paraoffsets=[], **kwargs):
        """
        Args:
        - system (openmm.System): defines all forces in OpenMM
        - psf (openmm.app.CharmmPsfFile): OpenMM representation of a psf file
        - paragroups (list of nbgroup): Instances that define changed nonbonded parameters.
          If paragroups=[] (the default), the energies are evaluated for the original
          force field parameters.
        - paraoffsets(list of float): Offsets for nonbonded parameters (see OpenMM documentation
          of NonbondedForce)
        """
        topology = md.Topology.from_openmm(psf.topology)
        self.psf = psf
        self.system = system
        self.topology = topology
        self.name = "Potential Energy (kJ/mole)"
        self.paragroups = paragroups
        if self.paragroups:
            assert paraoffsets
            self.paraoffsets = paraoffsets
        else:
            # 'empty_drude_parameter' is defined in paragroup
            self.paragroups = [empty_drude_parameter()]
            self.paraoffsets = [0]
        self.contexts = create_drude_contexts(
            self.system, self.psf, self.paragroups,
            self.paraoffsets, **kwargs
        )

    def __call__(self, traj, **kwargs):
        """
        Args:
        - traj (mdtraj.Trajectory): the trajectory object to get energies from
        - n_frames (kwargs): the first n frames of the trajectory to calculate energies
        """
        all_energies = []

        # create empty list for energies data
        for _ in range(len(self.paragroups)):
            all_energies.append([])
        if "n_frames" not in kwargs:
            n_frames = traj.n_frames
        else:
            n_frames = kwargs["n_frames"]
        for frame in range(n_frames):
            box_vectors = lengths_and_angles_to_box_vectors(
                *traj.unitcell_lengths[frame], *traj.unitcell_angles[frame]
            )
            for g in range(len(self.paragroups)):
                self.contexts[g].setPeriodicBoxVectors(*box_vectors)
                self.contexts[g].setPositions(traj.xyz[frame])
                all_energies[g].append(
                    self.contexts[g].getState(
                        getEnergy=True
                    ).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                )
        # attention: the energies are returned, but the parameter change is
        # not recorded here
        if len(all_energies) == 1:
            return np.reshape(np.array(all_energies), (-1, 1))
        else:
            energy_array = np.swapaxes(np.array(all_energies), 0, 1)
            return energy_array


class ParameterForce(object):
    """Evaluate forces of trajectory frames in OpenMM.
    """
    def __init__(self, system, psf, paragroups=[], paraoffsets=[], **kwargs):
        """
        Args:
        - system (openmm.System): defines all forces in OpenMM
        - psf (openmm.app.CharmmPsfFile): OpenMM representation of a psf file
        - paragroups (list of nbgroup): Instances that define changed nonbonded parameters.
          If paragroups=[] (the default), the energies are evaluated for the original
          force field parameters.
        - paraoffsets(list of float): Offsets for nonbonded parameters (see OpenMM documentation
          of NonbondedForce)
        """
        topology = md.Topology.from_openmm(psf.topology)
        self.system = system
        self.topology = topology
        self.name = "Forces on Particles (kJ/mole/A)"
        self.paragroups = paragroups
        if self.paragroups:
            assert paraoffsets
            self.paraoffsets = paraoffsets
        else:
            # 'empty_nbgroup' is defined in paragroup
            self.paragroups = [empty_nbgroup()]
            self.paraoffsets = [0]

    def __call__(self, traj, **kwargs):

        all_forces = []

        if "use_new_method" in kwargs:
            use_new_method = kwargs["use_new_method"]
        else:
            use_new_method = False
        if use_new_method:
            contexts = get_contexts_using_offset(
                self.system, self.topology,
                self.paragroups, self.paraoffsets, **kwargs
            )
        else:
            contexts = get_contexts(
                self.system, self.psf, self.paragroups,
                self.paraoffsets, **kwargs
            )

        # create empty list for energies data
        for _ in range(len(self.paragroups)):
            all_forces.append([])
        if "n_frames" not in kwargs:
            # Umm, can we use warning? Do we really need this? --YYL
            print("Warning: trajectory might be too long for forces storage!")
            n_frames = traj.n_frames
        else:
            n_frames = kwargs["n_frames"]
        for frame in range(n_frames):
            box_vectors = lengths_and_angles_to_box_vectors(
                *traj.unitcell_lengths[frame], *traj.unitcell_angles[frame]
            )
            for g in range(len(self.paragroups)):
                contexts[g].setPeriodicBoxVectors(*box_vectors)
                contexts[g].setPositions(traj.xyz[frame])
                all_forces[g].append(
                    contexts[g].getState(getForces=True).getForces().
                        value_in_unit(u.kilojoule_per_mole/u.angstrom)
                )
        # Attention: the energies are returned, but the parameter change is
        # not recorded here
        return all_forces


class ParametersEnergy(object):
    """Evaluate potential energies of trajectory frames in OpenMM.
    """
    def __init__(self, system, psf, paras_groups=[[]], paras_offsets=[[]], **kwargs):
        """
        Args:
        - system (openmm.System): defines all forces in OpenMM
        - psf (openmm.app.CharmmPsfFile): OpenMM representation of a psf file
        - paragroups (list of nbgroup): Instances that define changed nonbonded parameters.
          If paragroups=[] (the default), the energies are evaluated for the original
          force field parameters.
        - paraoffsets(list of float): Offsets for nonbonded parameters (see OpenMM documentation
          of NonbondedForce)
        """
        topology = md.Topology.from_openmm(psf.topology)
        self.psf = psf
        self.system = system
        self.topology = topology
        self.name = "Potential Energy (kJ/mole)"
        self.paras_groups = paras_groups
        if self.paras_groups[0]:
            assert paras_offsets
            self.paras_offsets = paras_offsets
        else:
            # 'empty_nbgroup' is defined in paragroup
            self.paras_groups = [[empty_nbgroup()]]
            self.paras_offsets = [[0]]
        if "use_new_method" in kwargs:
            use_new_method = kwargs["use_new_method"]
        else:
            use_new_method = False
        if use_new_method:
            self.contexts = get_contexts_using_offset_defult_platform(
                self.system, self.topology, self.paras_groups,
                self.paras_offsets
            )
        else:
            self.contexts = get_contexts_default_platform(
                self.system, self.psf, self.paras_groups, self.paras_offsets
            )

    def __call__(self, traj, **kwargs):
        """
        Args:
        - traj (mdtraj.Trajectory): the trajectory object to get energies from
        - n_frames (kwargs): the first n frames of the trajectory to calculate energies
        """
        all_energies = []

        # create empty list for energies data
        for _ in range(len(self.paras_groups)):
            all_energies.append([])
        if "n_frames" not in kwargs:
            n_frames = traj.n_frames
        else:
            n_frames = kwargs["n_frames"]
        for frame in range(n_frames):
            box_vectors = lengths_and_angles_to_box_vectors(
                *traj.unitcell_lengths[frame], *traj.unitcell_angles[frame]
            )
            for g in range(len(self.paras_groups)):
                self.contexts[g].setPeriodicBoxVectors(*box_vectors)
                self.contexts[g].setPositions(traj.xyz[frame])
                all_energies[g].append(
                    self.contexts[g].getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                )
        # attention: the energies are returned, but the parameter change is not recorded here
        if len(all_energies)==1:
            return np.reshape(np.array(all_energies), (-1,1))
        else:
            energy_array = np.swapaxes(np.array(all_energies), 0, 1)
            return energy_array


class SimEnergiesOpenmm():
    def __init__(self):
        self.name = "Simulation energies from openmm output (kj/mol)"
    
    def __call__(self, out_file):
        data = np.genfromtxt(out_file, skip_header = 1)
        data = data[:,2]
        return data


class PressureVolume(object):
    def __init__(self, pressure=1.0*u.atmosphere):
        self.name = "pV (units)"
        self.pressure = pressure

    def __call__(self, traj):
        return (self.pressure *
                (traj.unitcell_lengths[:,0]*traj.unitcell_lengths[:,1]*traj.unitcell_lengths[:,2]
                )*(u.nanometer**3) * u.AVOGADRO_CONSTANT_NA
               ).value_in_unit(u.kilojoules_per_mole)


class MdtrajRdf(object):
    def __init__(self, psf, name, atom_selections, r_range=[0, 1]):
        topology = md.Topology.from_openmm(psf.topology)
        self.topology = topology
        self.r_range = r_range
        self.name = name
        self.sele_words = atom_selections
        self.atom1 = self.topology.select(atom_selections[0])
        self.atom2 = self.topology.select(atom_selections[1])
        grid1, grid2 = np.meshgrid(self.atom1, self.atom2)
        self.pairs = []
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                self.pairs.append([grid1[i][j], grid2[i][j]])
        self.count_traj = 0

    def __call__(self, traj, **kwargs):
        self.count_traj += 1
        print('Calculating <{}> rdf for trajectory {} ...'.format(
            self.name, self.count_traj))
        a, b = md.compute_rdf(traj, self.pairs, self.r_range)
        self.radius = a
        if self.count_traj == 1:
            self.rdf = b
        else:
            self.rdf = (self.rdf * (self.count_traj - 1) + b) / self.count_traj
        return a, b

    def plot(self, cut, color='blue'):
        from matplotlib import pyplot as plt
        plt.plot(self.radius[:int(cut / 0.005)], self.rdf[:int(cut / 0.005)], color=color)
        plt.title(self.name)
        plt.show()
