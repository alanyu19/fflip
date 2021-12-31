# -*- coding: utf-8 -*-

import warnings
from fflip.omm.paragroup import *


def add_a_new_group(existing_groups, new_group):
    # make sure that there is no repeated definition
    for g in existing_groups:
        if g.par_type == new_group.par_type and \
                g.center_names == new_group.center_names:
            raise Exception(new_group)
    existing_groups.append(new_group)


class CharmmGroup(object):
    def __init__(self, num_atom_category, atoms, charges, half_r_mins,
                 epsilons, add_charge_group, atoms_same_charge,
                 neighbors, add_lj_group, atoms_same_lj, **kwargs):
        self.num_atom_category = num_atom_category
        self.atoms = atoms  # should be list inside list structure
        self.charges = charges
        self.half_r_mins = half_r_mins
        self.epsilons = epsilons
        # Charge
        self.add_charge_group = add_charge_group
        self.atoms_same_charge = atoms_same_charge
        self.neighbors = neighbors
        # LJ
        self.add_lj_group = add_lj_group
        self.atoms_same_lj = atoms_same_lj
        if "exclusion" in kwargs:
            self.exclusion = kwargs["exclusion"]
        else:
            self.exclusion = []


class Lipid(object):
    def __init__(self, charmm_group_list, name, **kwargs):
        self.charmm_group_list = charmm_group_list
        self.num_charmm_groups = len(self.charmm_group_list)
        self.name = name
        if "lipname" in kwargs:
            warnings.warn(
                "'lipname' will not be supported in future versions, please use 'name'!"
            )
            self.name = kwargs["lipname"]
        self.cgroups = []
        for group in self.charmm_group_list:
            self.cgroups.append(group)
            
    def level_print(self, content, level, print_level):
        if print_level >= level:
            print(content)
        else:
            pass

    def parse_groups(self, id_allowed='all', print_level=0):
        gs = []
        for counter, chm_gp in enumerate(self.cgroups):
            if id_allowed is not 'all':
                assert isinstance(id_allowed, list)
                if counter not in id_allowed:
                    continue
            self.level_print(
                "Creating 'NonbondedGroup's for {} group {} ... ".format(
                    self.name, counter+1
                ), 1, print_level
            )
            for i in range(chm_gp.num_atom_category):
                # LJ
                if chm_gp.add_lj_group[i]:
                    atom_list=[]
                    for atom in chm_gp.atoms[i]:
                        atom_list.append(atom)
                    for atom in chm_gp.atoms_same_lj[i]:
                        atom_list.append(atom)
                    add_a_new_group(
                        gs, NonbondedGroup(
                            lipid_name=self.name,
                            par_type="sigma", center_names=atom_list,
                            original_p=round(
                                (1/2)**(1/6) * chm_gp.half_r_mins[i] * 0.2, 5
                            ),
                            targeted_range=[0.8, 1.2]
                        )
                    )
                    add_a_new_group(
                        gs, NonbondedGroup(
                            lipid_name=self.name.lower(),
                            par_type="epsilon", center_names=atom_list,
                            original_p=round(chm_gp.epsilons[i] * 4.184, 4),
                            targeted_range=[0.8, 1.2]
                        )
                    )
                else:
                    self.level_print(
                        "Skipping LJ parameters for {} ...".format(
                        chm_gp.atoms[i]
                        ), 2, print_level
                    )
                # Charge
                if chm_gp.add_charge_group[i]:
                    center_names=[]; nb_names=[]
                    # Create and fill in the atoms with the current CHARMM group
                    for atom in chm_gp.atoms[i]:
                        center_names.append(atom)
                    for nb_index in chm_gp.neighbors[i]:
                        for atom in chm_gp.atoms[nb_index]:
                            nb_names.append(atom)
                    # Find any associated charges in other groups
                    if (chm_gp.atoms_same_charge is not None) and \
                            (i not in chm_gp.exclusion):
                        for atom in chm_gp.atoms_same_charge[i]:
                            center_names.append(atom)
                        for nb_index in chm_gp.neighbors[i]:
                            for atom in chm_gp.atoms_same_charge[nb_index]:
                                nb_names.append(atom)
                    # Now we have the name lists, create the 'NonbondedGroup' for charge:
                    roc = []
                    ron = []
                    # I think this decision will influence the speed of the
                    # optimization, think twice before using ...
                    for n in range(len(nb_names)):
                        ron.append(round(-len(center_names)/len(nb_names), 5))
                    # for r in range(len(roc)):
                    #     roc.append(1)
                    add_a_new_group(
                        gs, NonbondedGroup(
                            lipid_name=self.name.lower(),
                            par_type="charge", center_names=center_names,
                            neighbors=nb_names,
                            original_p=chm_gp.charges[i],
                            targeted_range=[round(chm_gp.charges[i] - 0.2, 2),
                                            round(chm_gp.charges[i] + 0.2, 2)],
                            roc=roc, ron=ron
                        )
                    )
                else:
                    self.level_print(
                        "Skipping charge(s) for {} ...".format(
                            chm_gp.atoms[i]
                        ), 2, print_level
                    )
            self.level_print("", 1, print_level)
        self.level_print(
            "Total {} NonbondedGroups created for {}\n".format(len(gs), self.name),
            1, print_level
        )
        return gs


