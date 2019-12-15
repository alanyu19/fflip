# -*- coding: utf-8 -*-

import numpy as np
from fflip.analysis.util import construct_scd_bond, gen_scd_pairs


class OrderParameterFactory(object):
    """
    A tool that generate order parameter calculations automatically for a
    given topology
    """
    def __init__(self, topology, special_carbons_for_splitting={}):
        """
        Args:
            topology: the mdtraj topology
            special_carbons_for_splitting: carbon names for splitting
        """
        self.topology = topology
        self.scd_res_dict = construct_scd_bond(topology)
        self.bonds_for_residues = construct_scd_bond(topology)
        self.scd_pairs_for_residues = gen_scd_pairs(
            self.bonds_for_residues,
            special_carbons_for_residues=special_carbons_for_splitting
        )

    def __call__(self):
        """
        Returns: all OrderParameterCalculation(s)
        """
        calculations = []
        i = 0
        for residue in self.scd_pairs_for_residues:
            for atom_pair in self.scd_pairs_for_residues[residue]:
                i += 1
                opc = OrderParameterCalculation(
                    self.topology, residue, atom_pair[0], atom_pair[1]
                )
                calculations.append(opc)
        return calculations

    def __iter__(self):
        """
        Returns: OrderParameterCalculation(s)
        """
        i = 0
        for residue in self.scd_pairs_for_residues:
            for atom_pair in self.scd_pairs_for_residues[residue]:
                i += 1
                opc = OrderParameterCalculation(
                    self.topology, residue, atom_pair[0], atom_pair[1]
                )
                opc.i = i
                yield opc

    def __len__(self):
        l = 0
        for residue in self.scd_pairs_for_residues:
            for _ in self.scd_pairs_for_residues[residue]:
                l += 1
        return l


class OrderParameterCalculation(object):
    """
    Calculate the order parameter (average over all molecules) of a (group
    of) chemical bond
    TODO: THE BILAYER NORMAL DIRECTION
    """
    def __init__(self, topology, residue, atom1, atom2, ):
        """
        Args:
            topology: a topology object of mdtraj
            residue: str, the residue name
            atom1: str, the center atom name
            atom2: the list of hydrogens
        """
        self.n_frames = 0.0
        self.scd = 0.0
        self.topology = topology
        self.residue = residue
        self.atom1 = atom1
        self.atom2 = atom2

    def __call__(self, traj):
        sele1 = self.topology.select(
            "resname {} and name {}".format(self.residue, self.atom1)
        )
        center_atom = traj.xyz[:, sele1]
        for i, atom in enumerate(self.atom2):
            sele2 = self.topology.select(
                "resname {} and name {}".format(self.residue, atom)
            )
            the_other_atom = traj.xyz[:, sele2]
            vectors = the_other_atom - center_atom
            if i==0:
                scd = vectors[:,:,2] * vectors[:,:,2] / (
                    vectors[:,:,0] * vectors[:,:,0] +
                    vectors[:,:,1] * vectors[:,:,1] +
                    vectors[:,:,2] * vectors[:,:,2]
                )
            else:
                scd += vectors[:,:,2] * vectors[:,:,2] / (
                    vectors[:,:,0] * vectors[:,:,0]
                    + vectors[:,:,1] * vectors[:,:,1]
                    + vectors[:,:,2] * vectors[:,:,2]
                )
        scd = scd/(i+1)

        # update the average scd
        self.scd = self.n_frames * self.scd + (-1.5 * np.mean(scd) + 0.5) * \
        traj.n_frames
        self.n_frames += traj.n_frames
        self.scd /= self.n_frames

        # return the average scd for every frame
        return -1.5 * np.mean(scd, axis=1) + 0.5
