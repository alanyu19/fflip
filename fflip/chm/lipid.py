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


class CharmmGroup:
    def __init__(self, **kwargs):
        self.num_atom_category = kwargs["num_atom_category"]
        self.atoms = kwargs["atoms"]  # should be list inside list structure
        self.charges = kwargs["charges"]
        self.half_r_mins = kwargs["half_r_mins"]
        self.epsilons = kwargs["epsilons"]
        self.add_lj_nbgroup = kwargs["add_lj_nbgroup"]
        self.atoms_same_lj = kwargs["atoms_same_lj"]
        self.add_charge_nbgroup = kwargs["add_charge_nbgroup"]
        self.neighbors = kwargs["neighbors"]
        self.cooperators = kwargs["cooperators"]
        self.atoms_same_charge = kwargs["atoms_same_charge"]
        if "exclusion" in kwargs:
            self.exclusion = kwargs["exclusion"]
    

class Lipid:
    def __init__(self, **kwargs):
        assert "charmm_group_list" in kwargs
        self.charmm_group_list = kwargs["charmm_group_list"]
        self.num_charmm_groups = len(self.charmm_group_list)
        assert "lipname" in kwargs or "name" in kwargs
        if "name" in kwargs:
            self.name = kwargs["name"]
        elif "lipname" in kwargs:
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

    def parse_nbgroups(self, groups='all', print_level=0):
        gs = []
        for counter, chm_gp in enumerate(self.cgroups):
            if groups is not 'all':
                assert isinstance(groups, list)
                if counter not in groups:
                    continue
            self.level_print(
                "Creating 'nbgroup's for {} group {} ... ".format(
                    self.name, counter+1
                ), 1, print_level
            )
            for i in range(chm_gp.num_atom_category):
                # LJ
                if chm_gp.add_lj_nbgroup[i]:
                    atom_list=[]
                    for atom in chm_gp.atoms[i]:
                        atom_list.append(atom)
                    for atom in chm_gp.atoms_same_lj[i]:
                        atom_list.append(atom)
                    add_a_new_group(
                        gs, nbgroup(
                            par_type="sigma", center_names=atom_list,
                            original_p=round(
                                (1/2)**(1/6) * chm_gp.half_r_mins[i] * 0.2, 5
                            ),
                            targeted_range=[0.8, 1.2]
                        )
                    )
                    add_a_new_group(
                        gs, nbgroup(
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
                if chm_gp.add_charge_nbgroup[i]:
                    center_names=[]; nb_names=[]; coop_names=[]
                    # Create and fill in the atoms with the current CHARMM group
                    for atom in chm_gp.atoms[i]:
                        center_names.append(atom)
                    for nb_index in chm_gp.neighbors[i]:
                        for atom in chm_gp.atoms[nb_index]:
                            nb_names.append(atom)
                    if chm_gp.cooperators is not None:
                        for coop_index in chm_gp.cooperators[i]:
                            for atom in chm_gp.atoms[coop_index]:
                                coop_names.append(atom)
                    # Find any associated charges in other groups
                    if (chm_gp.atoms_same_charge is not None) and \
                            (i not in chm_gp.exclusion):
                        for atom in chm_gp.atoms_same_charge[i]:
                            center_names.append(atom)
                        for nb_index in chm_gp.neighbors[i]:
                            for atom in chm_gp.atoms_same_charge[nb_index]:
                                nb_names.append(atom)
                        if chm_gp.cooperators is not None:
                            for coop_index in chm_gp.cooperators[i]:
                                for atom in chm_gp.atoms_same_charge[
                                    coop_index
                                ]:
                                    coop_names.append(atom)
                    # Now we have the name lists, create the 'nbgroup' for charge:
                    roc = []
                    ron = []
                    # I think this decision will influence the speed of the
                    # optimization, think twice before using ...
                    for n in range(len(nb_names)):
                        ron.append(round(-len(center_names)/len(nb_names), 5))
                    # for r in range(len(roc)):
                    #     roc.append(1)
                    add_a_new_group(
                        gs, nbgroup(
                            par_type="charge", center_names=center_names,
                            cooperators=coop_names, neighbors=nb_names,
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
            "Total {} nbgroups created for {}\n".format(len(gs), self.name),
            1, print_level
        )
        return gs


