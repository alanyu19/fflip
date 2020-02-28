# -*- coding: utf-8 -*-

import numpy as np
from fflip.analysis.util import hard_up_low_atoms


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    one = np.linalg.norm(vector, axis=2)
    stacked = np.stack((one, one, one), axis=2)
    return vector / stacked


def angle(v1, v2):
    assert v1.shape == v2.shape
    product = []
    for frame in range(v1.shape[0]):
        product.append([])
        for residue in range(v2.shape[1]):
            dot = np.dot(v1[frame, residue, :], v2[frame, residue, :])
            dot = np.minimum(1.0, dot)
            dot = np.maximum(-1.0, dot)
            product[frame].append(dot)
    product = np.array(product)
    return np.arccos(product)


class AngleDistribution(object):
    """
    A tool that generate order parameter calculations automatically for a
    given topology
    """
    def __init__(self, topology, first_vector, constant_second_vector=True,
                 second_vector=[0, 0, 1], sep_leaflet=False, recipe=None,
                 nbins=90):
        """
        Args:
            topology: the mdtraj topology
            first_vector: the atom pair of the vector
            second_vector: can be a constant vector or another atom pair
        """
        self.nbins = nbins
        self.atom1 = first_vector[0]
        self.atom2 = first_vector[1]
        self.constant_second_vector = constant_second_vector
        if not constant_second_vector:
            self.atom3 = second_vector[0]
            self.atom4 = second_vector[1]
        else:
            self.constant_vector = second_vector
        self.n_frames = 0.0
        self.topology = topology
        self.sep_leaflet = sep_leaflet
        self.recipe = recipe
        self.sele1 = self.topology.select(self.atom1)
        self.sele2 = self.topology.select(self.atom2)
        if not self.sep_leaflet:
            self.distribution = 0.0
        else:
            self.distribution_up = 0.0
            self.distribution_low = 0.0
            self.sele1_up, self.sele1_low = hard_up_low_atoms(
                self.sele1, self.recipe
            )
            self.sele2_up, self.sele2_low = hard_up_low_atoms(
                self.sele2, self.recipe
            )
        if not constant_second_vector:
            self.sele3 = self.topology.select(self.atom3)
            self.sele4 = self.topology.select(self.atom4)
            if not self.sep_leaflet:
                pass
            else:
                self.sele3_up, self.sele3_low = hard_up_low_atoms(
                    self.sele3, self.recipe
                )
                self.sele4_up, self.sele4_low = hard_up_low_atoms(
                    self.sele4, self.recipe
                )

    def __call__(self, traj):
        if not self.sep_leaflet:
            atom1 = traj.xyz[:, self.sele1]
            atom2 = traj.xyz[:, self.sele2]
            vector1 = atom1 - atom2
            if self.constant_second_vector:
                vector2 = np.tile(
                    self.constant_vector, (vector1.shape[0], vector1.shape[1])
                )  # duplicated by number of frames and residues
            else:
                atom3 = traj.xyz[:, self.sele3]
                atom4 = traj.xyz[:, self.sele4]
                vector2 = atom3 - atom4
            # normalize all the vectors and do the inner product
            vector1 = unit_vector(vector1)
            vector2 = unit_vector(vector2)
            angles = angle(vector1, vector2)
            distribution, edges = np.histogram(
                angles, bins=self.nbins, range=(0, np.pi)
            )
            self.edges = edges
            self.distribution = self.n_frames * self.distribution + \
                distribution * traj.n_frames
            self.n_frames += traj.n_frames
            self.distribution /= self.n_frames
        else:
            atom1_up = traj.xyz[:, self.sele1_up]
            atom2_up = traj.xyz[:, self.sele2_up]
            atom1_low = traj.xyz[:, self.sele1_low]
            atom2_low = traj.xyz[:, self.sele2_low]
            vector1_up = atom1_up - atom2_up
            vector1_low = atom1_low - atom2_low
            if self.constant_second_vector:
                vector2_up = np.tile(
                    self.constant_second_vector, (
                        vector1_up.shape[0], vector1_up.shape[1]
                    )  # duplicated by number of frames and residues
                )
                vector2_low = np.tile(
                    -self.constant_second_vector, (
                        vector1_low.shape[0], vector1_low.shape[1]
                    )  # duplicated by number of frames and residues
                )
            else:
                atom3_up = traj.xyz[:, self.sele3_up]
                atom4_up = traj.xyz[:, self.sele4_up]
                atom3_low = traj.xyz[:, self.sele3_low]
                atom4_low = traj.xyz[:, self.sele4_low]
                vector2_up = atom3_up - atom4_up
                vector2_low = atom3_low - atom4_low
            # normalize all the vectors and do the inner product
            vector1_up = unit_vector(vector1_up)
            vector2_up = unit_vector(vector2_up)
            vector1_low = unit_vector(vector1_low)
            vector2_low = unit_vector(vector2_low)
            angles_up = angle(vector1_up, vector2_up)
            angles_low = angle(vector1_low, vector2_low)
            distribution_up, edges = np.histogram(
                angles_up, bins=self.nbins, range=(0, np.pi)
            )
            distribution_low, edges = np.histogram(
                angles_low, bins=self.nbins, range=(0, np.pi)
            )
            self.edges = edges
            self.distribution_up = self.n_frames * self.distribution_up + \
                distribution_up * traj.n_frames
            self.distribution_low = self.n_frames * self.distribution_low + \
                distribution_low * traj.n_frames
            self.n_frames += traj.n_frames
            self.distribution_up /= self.n_frames
            self.distribution_low /= self.n_frames
