# -*- coding: utf-8 -*-

from coffe.omm.paragroup import *
from fflip.drude.exceptions import *


def add_a_new_group(existing_groups, new_group):
    # make sure that there is no repeated definition
    for g in existing_groups:
        if g.par_type == new_group.par_type and \
                g.center_names == new_group.center_names:
            raise GroupRepeatedError(new_group)
    existing_groups.append(new_group)


class DrudeChargeGroup:
    def __init__(self, id_, atom_groups, drude_particles,
                 charges, alphas, tholes, add_group, neighbors,
                 atoms_same_charge, atoms_same_alpha, atoms_same_thole=None,
                 add_alpha=None, add_thole=None, exclusion=None):
        """
        The class for defining drude FF charge group.
        Args:
            id: int, the id within the lipid
            atom_groups: list of list, each sub-list contains atoms sharing the
            same charge.
            drude_particles: list of the Drude Particle names
            charges: list, each corresponds to the sub-list in atoms.
            alphas: list, each corresponds to the sub-list in atoms.
            tholes: list, each corresponds to the sub-list in atoms.
            add_group: bool, if parse this group.
            neighbors: list of list, the atoms to exchange charge.
            atoms_same_charge: list of list, atoms in other groups that
            add_alpha:
            atoms_same_alpha:
            must use the same charge.
            exclusion: bool, if ignore the whole group
        Example:
            dcg = DrudeChargeGroup()
        """
        self.id = id_
        self.atom_groups = atom_groups  # should be list inside list structure
        self.drude_particles = drude_particles
        self.charges = charges
        self.alphas = alphas
        self.tholes = tholes
        self.add_group = add_group
        self.neighbors = neighbors
        self.atoms_same_charge = atoms_same_charge
        self.atoms_same_alpha = atoms_same_alpha
        if atoms_same_thole is None:
            self.atoms_same_thole = self.atoms_same_alpha
        else:
            self.atoms_same_thole = atoms_same_thole
        if add_alpha is None:
            self.add_alpha = add_group
        else:
            self.add_alpha = add_alpha
        if add_thole is None:
            self.add_thole = self.add_alpha
        else:
            self.add_thole = add_thole
        if exclusion is None:
            self.exclusion = []
        else:
            self.exclusion = exclusion


class LJGroup:
    def __init__(self, _id, atom_type_dict, half_r_mins=None, epsilons=None):
        """
        The class for defining all parametrizable LJ interactions, NBFIX is
        not included in this class and can be set separately using the
        DrudeNBFIX class
        Args:
            id: int, the id within the lipid
            atom_type_dict: a dict contains all atom types and their related
            atom names.
            half_r_mins: can be the original parameters for half_r_mins (A)
            epsilons: can be the original parameters for epsilon (kcal/mol)
        """
        self.id = _id
        self.atom_type_dict = atom_type_dict
        self.half_r_mins = half_r_mins
        self.epsilons = epsilons


class DrudeNBFIX:
    pass


class DrudeLipid:
    def __init__(self, lipname, charge_groups=None, lj_groups=None):
        """
        DRUDE Lipid, NBFIX to come
        Args:
            lipname: the name of the lipid / lipid family
            charge_groups: as named, the DrudeChargeGroup class
            lj_groups: the LJGroup class
        """
        self.lipname = lipname
        self.charge_groups = charge_groups
        self.lj_groups = lj_groups
        if charge_groups:
            self.num_charge_groups = len(self.charge_groups)
        if lj_groups:
            self.num_lj_groups = len(self.lj_groups)
            
    @staticmethod
    def level_print(content, level, print_level):
        if print_level >= level:
            print(content)
        else:
            pass

    def parse_groups(self, print_level=0, chg_offset=0.2, alpha_offset=0.02,
                     thole_offset=0.02, id_allowed='all'):
        gs = []
        for counter, chggp in enumerate(self.charge_groups):
            if id_allowed is not 'all' and chggp.id not in id_allowed:
                continue
            self.level_print("Parsing {} charge group {} ...".format(
                self.lipname, counter + 1
            ), 1, print_level)
            # used to record the indexing of the parameter within the group
            internal_id = 0
            for i in range(len(chggp.atom_groups)):
                if chggp.add_group[i]:
                    center_names = []
                    nbr_names = []
                    # Create and fill in the atoms
                    for atom in chggp.atom_groups[i]:
                        center_names.append(atom)
                    for nbr_index in chggp.neighbors[i]:
                        for atom in chggp.atom_groups[nbr_index]:
                            nbr_names.append(atom)
                    # Find any associated charges in other charge groups
                    if (chggp.atoms_same_charge is not None) and \
                            (i not in chggp.exclusion):
                        for atom in chggp.atoms_same_charge[i]:
                            center_names.append(atom)
                        for nbr_index in chggp.neighbors[i]:
                            for atom in chggp.atoms_same_charge[nbr_index]:
                                nbr_names.append(atom)
                    # Now we have the name lists, create the 'gtcnp' for charge:
                    ron = []
                    # I think this decision will influence the speed of the
                    # optimization, think twice before using ...
                    for n in range(len(nbr_names)):
                        ron.append(round(-len(center_names)/len(nbr_names), 5))
                    # for r in range(len(roc)):
                    #     roc.append(1)
                    internal_id += 1
                    add_a_new_group(
                        gs, DrudeParameter(
                            lipidname=self.lipname.lower(),
                            cgid=chggp.id,
                            internal_id=internal_id,
                            par_type="charge",
                            center_names=center_names,
                            original_p=chggp.charges[i],
                            targeted_range=[
                                round(chggp.charges[i] - chg_offset, 4),
                                round(chggp.charges[i] + chg_offset, 4)
                            ],
                            neighbors=nbr_names,
                            drude_particles=None,
                            ron=ron
                        )
                    )
                else:
                    self.level_print(
                        "Skipping charge(s) for {} ...".format(
                            chggp.atom_groups[i]
                        ), 2, print_level
                    )
                if chggp.add_alpha[i]:
                    assert chggp.drude_particles[i] is not None
                    drude_atoms = []
                    heavy_atoms = []
                    assert len(chggp.atom_groups[i]) == \
                        len(chggp.drude_particles[i])
                    for hatom, datom in zip(
                            chggp.atom_groups[i], chggp.drude_particles[i]
                    ):
                        drude_atoms.append(datom)
                        heavy_atoms.append(hatom)
                    if chggp.atoms_same_alpha is not None:
                        # atom_combo is a set like (drude_atom, heavy_atom)
                        for atom_combo in chggp.atoms_same_alpha[i]:
                            drude_atoms.append(atom_combo[0])
                            heavy_atoms.append(atom_combo[1])
                    internal_id += 1
                    add_a_new_group(
                        gs, DrudeParameter(
                            lipidname=self.lipname.lower(),
                            cgid=chggp.id,
                            internal_id=internal_id,
                            par_type="alpha",
                            center_names=heavy_atoms,
                            original_p=chggp.alphas[i],
                            targeted_range=[
                                round(chggp.alphas[i] - alpha_offset, 9),
                                round(chggp.alphas[i] + alpha_offset, 9)
                            ],
                            drude_particles=drude_atoms
                        )
                    )
                else:
                    self.level_print(
                        "Skipping alpha(s) for {} ...".format(
                            chggp.atom_groups[i]
                        ), 2, print_level
                    )
                if chggp.add_thole[i]:
                    assert chggp.drude_particles[i] is not None
                    drude_atoms = []
                    heavy_atoms = []
                    assert len(chggp.atom_groups[i]) == \
                        len(chggp.drude_particles[i])
                    for hatom, datom in zip(
                            chggp.atom_groups[i], chggp.drude_particles[i]
                    ):
                        drude_atoms.append(datom)
                        heavy_atoms.append(hatom)
                    if chggp.atoms_same_thole is not None:
                        # atom_combo is a set like (drude_atom, heavy_atom)
                        for atom_combo in chggp.atoms_same_thole[i]:
                            drude_atoms.append(atom_combo[0])
                            heavy_atoms.append(atom_combo[1])
                    internal_id += 1
                    add_a_new_group(
                        gs, DrudeParameter(
                            lipidname=self.lipname.lower(),
                            cgid=chggp.id,
                            internal_id=internal_id,
                            par_type="thole",
                            center_names=heavy_atoms,
                            original_p=chggp.tholes[i],
                            targeted_range=[
                                round(chggp.tholes[i] - thole_offset, 9),
                                round(chggp.tholes[i] + thole_offset, 9)
                            ],
                            drude_particles=drude_atoms
                        )
                    )
                else:
                    self.level_print(
                        "Skipping thole(s) for {} ...".format(
                            chggp.atom_groups[i]
                        ), 2, print_level
                    )
            self.level_print("", 1, print_level)
        self.level_print(
            "Total {} DrudeParameters created for {}\n".format(
                len(gs), self.lipname
            ), 1, print_level
        )
        return gs
