# -*- coding: utf-8 -*-

import numpy as np
from fflip.analysis.util import (
    construct_scd_bond, gen_scd_pairs, hard_up_low_atoms
)


class OrderParameterFactory(object):
    """
    A tool that generate order parameter calculations automatically for a
    given topology
    """
    def __init__(self, topology, special_carbons_for_splitting={},
                 skip=None, sep_leaflet=False, recipe=None):
        """
        Args:
            topology: the mdtraj topology
            special_carbons_for_splitting: carbon names for splitting
        """
        self.topology = topology
        scd_res_dict = construct_scd_bond(topology)
        if skip is not None:
            for resskip in skip:
                scd_res_dict.pop(resskip)
        self.bonds_for_residues = scd_res_dict
        self.scd_pairs_for_residues = gen_scd_pairs(
            self.bonds_for_residues,
            special_carbons_for_residues=special_carbons_for_splitting
        )
        self.sep_leaflet = sep_leaflet
        self.recipe = recipe

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
                    self.topology, residue, atom_pair[0], atom_pair[1],
                    sep_leaflet=self.sep_leaflet, recipe=self.recipe
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
        length = 0
        for residue in self.scd_pairs_for_residues:
            for _ in self.scd_pairs_for_residues[residue]:
                length += 1
        return length


class OrderParameterCalculation(object):
    """
    Calculate the order parameter (average over all molecules) of a (group
    of) chemical bond
    TODO: THE BILAYER NORMAL DIRECTION
    """
    def __init__(self, topology, residue, atom1, atom2,
                 sep_leaflet=False, recipe=None):
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
        self.sep_leaflet = sep_leaflet
        self.recipe = recipe
        if self.sep_leaflet:
            self.scd_up = 0.0
            self.scd_low = 0.0

    def __call__(self, traj, skip=1, per_mol=False):
        if not self.sep_leaflet:
            sele1 = self.topology.select(
                "resname {} and name {}".format(self.residue, self.atom1)
            )
            center_atom = traj.xyz[::skip, sele1]
            for i, atom in enumerate(self.atom2):
                sele2 = self.topology.select(
                    "resname {} and name {}".format(self.residue, atom)
                )
                the_other_atom = traj.xyz[::skip, sele2]
                vectors = the_other_atom - center_atom
                if i == 0:
                    scd = vectors[:, :, 2] * vectors[:, :, 2] / (
                        vectors[:, :, 0] * vectors[:, :, 0] +
                        vectors[:, :, 1] * vectors[:, :, 1] +
                        vectors[:, :, 2] * vectors[:, :, 2]
                    )
                else:
                    scd += vectors[:, :, 2] * vectors[:, :, 2] / (
                        vectors[:, :, 0] * vectors[:, :, 0]
                        + vectors[:, :, 1] * vectors[:, :, 1]
                        + vectors[:, :, 2] * vectors[:, :, 2]
                    )
            scd = scd / (i + 1)

            # update the average scd
            self.scd = self.n_frames * self.scd + \
                (-1.5 * np.mean(scd) + 0.5) * traj.n_frames / skip
            self.n_frames += traj.n_frames / skip
            self.scd /= self.n_frames

            # return the scd for every frame
            if per_mol:
                return (-1.5 * scd + 0.5).flatten()
            else:
                return -1.5 * np.mean(scd, axis=1) + 0.5
        else:
            assert self.recipe is not None
            sele1 = self.topology.select(
                "resname {} and name {}".format(self.residue, self.atom1)
            )
            sele1_up, sele1_low = hard_up_low_atoms(sele1, self.recipe)
            center_atom_up = traj.xyz[::skip, sele1_up]
            center_atom_low = traj.xyz[::skip, sele1_low]
            for i, atom in enumerate(self.atom2):
                sele2 = self.topology.select(
                    "resname {} and name {}".format(self.residue, atom)
                )
                sele2_up, sele2_low = hard_up_low_atoms(sele2, self.recipe)
                other_atom_up = traj.xyz[::skip, sele2_up]
                other_atom_low = traj.xyz[::skip, sele2_low]
                vectors_up = other_atom_up - center_atom_up
                vectors_low = other_atom_low - center_atom_low
                if i == 0:
                    scd_up = vectors_up[:, :, 2] * vectors_up[:, :, 2] / (
                        vectors_up[:, :, 0] * vectors_up[:, :, 0] +
                        vectors_up[:, :, 1] * vectors_up[:, :, 1] +
                        vectors_up[:, :, 2] * vectors_up[:, :, 2]
                    )
                    scd_low = vectors_low[:, :, 2] * vectors_low[:, :, 2] / (
                        vectors_low[:, :, 0] * vectors_low[:, :, 0] +
                        vectors_low[:, :, 1] * vectors_low[:, :, 1] +
                        vectors_low[:, :, 2] * vectors_low[:, :, 2]
                    )
                else:
                    scd_up += vectors_up[:, :, 2] * vectors_up[:, :, 2] / (
                            vectors_up[:, :, 0] * vectors_up[:, :, 0] +
                            vectors_up[:, :, 1] * vectors_up[:, :, 1] +
                            vectors_up[:, :, 2] * vectors_up[:, :, 2]
                    )
                    scd_low += vectors_low[:, :, 2] * vectors_low[:, :, 2] / (
                            vectors_low[:, :, 0] * vectors_low[:, :, 0] +
                            vectors_low[:, :, 1] * vectors_low[:, :, 1] +
                            vectors_low[:, :, 2] * vectors_low[:, :, 2]
                    )
            scd_up = scd_up / (i + 1)
            scd_low = scd_low / (i + 1)

            # update the average scd
            self.scd_up = self.n_frames * self.scd_up + \
                (-1.5 * np.mean(scd_up) + 0.5) * traj.n_frames / skip
            self.scd_low = self.n_frames * self.scd_low + \
                (-1.5 * np.mean(scd_low) + 0.5) * traj.n_frames / skip
            self.n_frames += traj.n_frames / skip
            self.scd_up /= self.n_frames
            self.scd_low /= self.n_frames

            # return the scd for every frame
            if per_mol:
                return -1.5 * scd_up + 0.5, -1.5 * scd_low + 0.5
            else:
                return -1.5 * np.mean(scd_up, axis=1) + 0.5, -1.5 * np.mean(
                    scd_low, axis=1) + 0.5
