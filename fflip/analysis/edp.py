# -*- coding: utf-8 -*-

import numpy as np
import os
import mdtraj as md
from rflow.edp import *
from fflip.analysis.edp_util import *


class ElectronDensityFactory:
    """
    Do production calculation and grouping
    """
    def __init__(self, psf_file, traj_template,
                 com_selection=None, box_size=10,
                 save_to_folder='.'):
        """
        Args:
            psf_file: str, the psf file
            traj_template: str, the trajectory template
            com_selection: str, the center of mass selection words (ref to
            mdtraj for this)
            save_to_folder: str, the folder to save all the atom data,
            usually don't need to change it.
        """
        self.psf_file = psf_file
        self.traj_template = traj_template
        self.com_selection = com_selection
        self.top = make_topology(psf_file)
        self.residues = find_resnames_from_psf(psf_file)
        if not os.path.isdir(save_to_folder):
            os.mkdir(save_to_folder)
        self.save_to_folder = save_to_folder
        self.box_size = box_size

    def __call__(self, first, last):
        """
        Args:
            first: int, the first trajectory index
            last: int, the last trajectory index
        Returns:
            None
        """
        trajs = TrajectoryIterator(
            first_sequence=first, last_sequence=last,
            filename_template=self.traj_template,
            topology_file=self.psf_file,
            atom_selection="all", load_function=md.load_dcd
        )
        for res in self.residues:
            subdir = os.path.join(
                self.save_to_folder,
                './atoms_{}'.format(res.lower())
            )
            if not os.path.isdir(subdir):
                os.mkdir(subdir)
            blocks_dir = os.path.join(subdir, 'blocks')
            if not os.path.isdir(blocks_dir):
                os.mkdir(blocks_dir)
            atoms = find_atoms_from_psf(self.psf_file, res)
            for atom in atoms:
                edc = ElectronDensityCalculator(
                    atom_selection="name {}".format(atom.upper()),
                    com_selection=self.com_selection,
                    box_length_fixed=self.box_size,
                    topology_file=self.psf_file
                )
                for traj in trajs:
                    edc(traj)
                np.savetxt(
                    os.path.join(blocks_dir, "{}_{}.dat".format(atom, first)),
                    edc.adp_angstrom
                )
                np.savetxt(
                    os.path.join(blocks_dir, "{}_{}_edp.dat".format(atom, first)),
                    edc.edp_angstrom
                )
            z_bin_file = os.path.join(subdir, 'z.txt')
            if not os.path.isfile(z_bin_file):
                np.savetxt(z_bin_file, np.linspace(
                    -edc.box_length_fixed / 2,
                    edc.box_length_fixed / 2,
                    edc.num_bins
                )
                       )