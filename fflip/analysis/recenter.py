import mdtraj as md
import numpy as np
import os
from rflow.trajectory import TrajectoryIterator, center_of_mass_of_selection
from rflow.openmm.app import CharmmPsfFile


class TrajectoryRecenter():
    def __init__(self, psf_file, input_template, first, last):
        self.input_template = input_template
        self.first = first
        self.last = last
        self.psf_file = psf_file
        self.psf = CharmmPsfFile(psf_file)
        self.topology = md.Topology.from_openmm(self.psf.topology)
        for i in range(self.first, self.last+1):
            assert os.path.isfile(input_template.format(i)), \
            "missing file " + input_template.format(i)
        molecules = self.topology.find_molecules()
        # assume all lipid molecules have at least 50 atoms
        self.lipid_moles = [mole for mole in molecules if len(mole) > 50]
        self.other_moles = [mole for mole in molecules if mole not in self.lipid_moles]
        self.lipid_molecules_atom_indices = [
            np.fromiter((a.index for a in mole), dtype=np.int32) for \
            mole in self.lipid_moles]
        self.other_molecules_atom_indices  = [
            np.fromiter((a.index for a in mole), dtype=np.int32) for \
            mole in self.other_moles]
    
    def check_membrane_at_center(self, trajectory, lipid_head_atoms,
        lipid_tail_atoms):
        # check if the membrane is put at the center already or not
        # trajectory must be a md.Trajectory
        head_atoms = self.topology.select(lipid_head_atoms)
        tail_atoms = self.topology.select(lipid_tail_atoms)
        # xyz np.ndarray, shape=(n_frames, n_atoms, 3)
        system_com_z = trajectory.xyz[0,:,2].mean()
        head_disp_abs = np.abs(
            trajectory.xyz[0, head_atoms, 2] - system_com_z)
        tail_disp_abs = np.abs(
            trajectory.xyz[0, tail_atoms, 2] - system_com_z)
        if head_disp_abs.mean() > tail_disp_abs.mean():
            return True
        else:
            return False

    def recenter_apart(self, trajectory):
        # unitcell_vectors = trajectory.unitcell_vectors
        unitcell_lengths = trajectory.unitcell_lengths
        # 2 is for z axis
        center_z = center_of_mass_of_selection(trajectory, "all", 2)
        # normalize z coordinates for faster computation
        z_normalized = trajectory.xyz[:, :, 2].transpose() - center_z.transpose()
        z_normalized /= trajectory.unitcell_lengths[:, 2].transpose()
        z_offsets_normalized = np.zeros(z_normalized.shape)
        for frame in range(trajectory.n_frames):
            for lipid_molecule in self.lipid_molecules_atom_indices:
                if z_normalized[lipid_molecule, frame].mean() < 0:
                    z_offsets_normalized[lipid_molecule, frame] += 1
            for other_molecule in self.other_molecules_atom_indices:
                if z_normalized[other_molecule, frame].mean() < 0:
                    z_offsets_normalized[other_molecule, frame] += 1
        z_normalized += z_offsets_normalized
        z = z_normalized * trajectory.unitcell_lengths[:, 2].transpose()
        z = z.transpose()
        trajectory.xyz[:, :, 2] = z
        return trajectory

    def recenter_intact(self, trajectory):
        # unitcell_vectors = trajectory.unitcell_vectors
        unitcell_lengths = trajectory.unitcell_lengths
        all_lipid_atom_indices = []
        for lipid in self.lipid_molecules_atom_indices:
            all_lipid_atom_indices += list(lipid)
        all_lipid_atom_indices = np.array(all_lipid_atom_indices)
        # 2 is for z axis
        center_z_lipid = trajectory.xyz[:, all_lipid_atom_indices, 2].mean(axis=1)
        # normalize z coordinates for faster computation
        z_normalized = trajectory.xyz[:, :, 2].transpose() - center_z_lipid.transpose()
        z_normalized /= trajectory.unitcell_lengths[:, 2].transpose()
        # system_com_z = trajectory.xyz[:,:,1].mean(axis=1)
        z_offsets_normalized = np.zeros(z_normalized.shape)
        for frame in range(trajectory.n_frames):
            for other_molecule in self.other_molecules_atom_indices:
                if z_normalized[other_molecule, frame].mean() < -0.5:
                    z_offsets_normalized[other_molecule, frame] += 1
                elif z_normalized[other_molecule, frame].mean() > 0.5:
                    z_offsets_normalized[other_molecule, frame] -= 1
        z_normalized += z_offsets_normalized
        # final_offset = z_normalized[all_lipid_molecule_atom_indices, :]
        z = z_normalized * trajectory.unitcell_lengths[:, 2].transpose()
        z = z.transpose()
        trajectory.xyz[:, :, 2] = z
        return trajectory
    
    def recenter_all_trajectories(self, output_format, lipid_head_atoms='name C2',
        lipid_tail_atoms='name C212'):
        for i in range(self.first, self.last+1):
            trajectory = md.load_dcd(
                self.input_template.format(i), top=self.psf_file)
            at_center_already = self.check_membrane_at_center(
                trajectory, lipid_head_atoms=lipid_head_atoms,
                lipid_tail_atoms=lipid_tail_atoms)
            if at_center_already:
                new_traj = self.recenter_intact(trajectory)
            else:
                new_traj = self.recenter_apart(trajectory)
            new_traj.save_dcd(output_format.format(i))


class TrajectoryRecenterNAMD():
    pass