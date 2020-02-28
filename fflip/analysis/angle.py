# -*- coding: utf-8 -*-

import numpy as np
from fflip.analysis.util import hard_up_low_atoms


class AngleDistribution(object):
    """
    A tool that generate order parameter calculations automatically for a
    given topology
    """
    def __init__(self, topology, first_vector, constant_second_vector=True,
                 second_vector=[0, 0, 1], sep_leaflet=False, recipe=None):
        """
        Args:
            topology: the mdtraj topology
            first_vector: the atom pair of the vector
            second_vector: can be a constant vector or another atom pair
        """
        self.atom1 = first_vector[0]
        self.atom2 = first_vector[1]
        self.constant_second_vector = constant_second_vector
        if not constant_second_vector:
            self.atom3 = second_vector[0]
            self.atom4 = second_vector[1]
        else:
            self.constant_vector = second_vector
        self.n_frames = 0.0
        self.distribution = 0.0
        self.topology = topology
        self.sep_leaflet = sep_leaflet
        self.recipe = recipe
        self.sele1 = self.topology.select(self.atom1)
        self.sele2 = self.topology.select(self.atom2)
        if not self.sep_leaflet:
            pass
        else:
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
                vector2 = self.constant_vector  # duplicated by number of frames
            else:
                atom3 = traj.xyz[:, self.sele3]
                atom4 = traj.xyz[:, self.sele4]
                vector2 = atom3 - atom4
            return vector1, vector2
            # normalize all the vectors and do the inner product
            #self.distribution = self.n_frames * self.distribution + \
            #           distribution * traj.n_frames
            #self.n_frames += traj.n_frames
            #self.distribution /= self.n_frames
        else:
            atom1_up = traj.xyz[:, self.sele1_up]
            atom2_up = traj.xyz[:, self.sele2_up]
            atom1_low = traj.xyz[:, self.sele1_low]
            atom2_low = traj.xyz[:, self.sele2_low]
            vector1_up = atom1_up - atom2_up
            vector1_low = atom1_low - atom2_low
            if self.constant_second_vector:
                vector2_up = self.constant_second_vector
                vector2_low = self.constant_second_vector
            # normalize and do inner product
