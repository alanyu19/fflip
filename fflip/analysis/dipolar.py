# -*- coding: utf-8 -*-

import numpy as np
from fflip.analysis.util import hard_up_low_atoms
from fflip.analysis.angle import angle, unit_vector


class DipolarCouplingCalculator:
    def __init__(self, topology, atom_pair, sep_leaflet=False, recipe=None):
        self.atom1 = atom_pair[0]
        self.atom2 = atom_pair[1]
        self.topology = topology
        self.sep_leaflet = sep_leaflet
        self.recipe = recipe
        self.sele1 = self.topology.select(self.atom1)
        self.sele2 = self.topology.select(self.atom2)
        self.traj_values = []
        
    def __call__(self, traj):
        if not self.sep_leaflet:
            atom1 = traj.xyz[:, self.sele1]
            atom2 = traj.xyz[:, self.sele2]
            vector1 = atom1 - atom2
            vector1copied = vector1
            vector2 = np.tile(
                [0, 0, 1],
                (vector1.shape[0], vector1.shape[1], 1)
            )  # duplicated by number of frames and residues
            # normalize all the vectors and do the inner product
            vector1 = unit_vector(vector1)
            vector2 = unit_vector(vector2)
            angles = angle(vector1, vector2)
            # angles_normalized = angles / np.pi
            distances = np.linalg.norm(vector1copied, axis=2) * 10  # to angstrom
            coupling = 12236.5 * np.mean(
                (3 * np.cos(angles)**2 - 1) * distances**(-3) / 2
            )
            self.traj_values.append(np.abs(coupling))
            self.avg = np.mean(self.traj_values)
            self.se = np.std(self.traj_values) / np.sqrt(len(self.traj_values))

