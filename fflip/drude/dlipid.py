# -*- coding: utf-8 -*-

from coffe.omm.paragroup import *
from fflip.drude.exceptions import *


def add_a_new_group(existing_groups, new_group):
    # make sure that there is no repeated definition
    for g in existing_groups:
        if g.par_type == new_group.par_type and g.center_names == new_group.center_names:
            raise GroupRepeatedError(new_group)
    existing_groups.append(new_group)


class DrudeChargeGroup:
    def __init__(self, atoms, charges, alphas, tholes, add_group, neighbors,
                 atoms_same_charge, exclusion=False):
        """
        The class for defining drude FF charge group.
        Args:
            atoms: list of list, each sub-list contains atoms sharing the same
            charge.
            charges: list, each corresponds to the sub-list in atoms.
            alphas: list, each corresponds to the sub-list in atoms.
            tholes: list, each corresponds to the sub-list in atoms.
            add_group: bool, if parse this group.
            neighbors: list of list, the atoms to exchange charge.
            atoms_same_charge: list of list, atoms in other groups that
            must use the same charge.
            exclusion: bool, if ignore the whole group
        Example:
            dcg = DrudeChargeGroup()
        """
        self.atoms = atoms  # should be list inside list structure
        self.charges = charges
        self.alphas = alphas
        self.tholes = tholes
        self.add_group = add_group
        self.neighbors = neighbors
        self.atoms_same_charge = atoms_same_charge
        self.exclusion = exclusion
    

class LJGroup:
    def __init__(self, atoms, half_r_mins, epsilons):
        self.atoms = atoms
        self.half_r_mins = half_r_mins
        self.epsilons = epsilons


class DrudeLipid:
    def __init__(self, lipname, charge_groups, lj_groups):
        self.lipname = lipname
        self.charge_groups = charge_groups
        self.lj_groups = lj_groups
        self.num_charge_groups = len(self.charge_groups)
        self.num_lj_groups = len(self.lj_groups)
            
    @staticmethod
    def level_print(content, level, print_level):
        if print_level >= level:
            print(content)
        else:
            pass

    def parse_groups(self, print_level=0):
        gs = []
        for counter, chm_gp in enumerate(self.charge_groups):
            self.level_print("Parsing {} charge group {} ...".format(
                self.lipname, counter + 1
            ), 1, print_level)
            for i in range(chm_gp.num_atom_category):
                # LJ
                if chm_gp.add_lj_gtcnp[i]:
                    atom_list=[]
                    for atom in chm_gp.atoms[i]:
                        atom_list.append(atom)
                    for atom in chm_gp.atoms_same_lj[i]:
                        atom_list.append(atom)
                    add_a_new_group(gs, gtcnp(par_type="sigma", center_names=atom_list,
                                              original_p=round((1/2)**(1/6) * chm_gp.half_r_mins[i] * 0.2, 5),
                                              targeted_range=[0.8, 1.2]
                                             )
                                   )
                    add_a_new_group(gs, gtcnp(par_type = "epsilon", center_names=atom_list,
                                              original_p=round(chm_gp.epsilons[i] * 4.184, 4), 
                                              targeted_range=[0.8, 1.2]
                                             )
                                   )
                else:
                    self.level_print("Skipping LJ parameters for {} ...".format(chm_gp.atoms[i]), 2, print_level)
                # Charge
                if chm_gp.add_charge_gtcnp[i]:
                    center_names=[]; nb_names=[]; coop_names=[]
                    # Create and fill in the atoms with the current CHARMM group
                    for atom in chm_gp.atoms[i]:
                        center_names.append(atom)
                    for nb_index in chm_gp.neighbors[i]:
                        for atom in chm_gp.atoms[nb_index]:
                            nb_names.append(atom)
                    if chm_gp.cooperators!=None:
                        for coop_index in chm_gp.cooperators[i]:
                            for atom in chm_gp.atoms[coop_index]:
                                coop_names.append(atom)
                    # Find any associated charges in other groups
                    if (chm_gp.atoms_same_charge!=None) and (i not in chm_gp.exclusion):
                        for atom in chm_gp.atoms_same_charge[i]:
                            center_names.append(atom)
                        for nb_index in chm_gp.neighbors[i]:
                            for atom in chm_gp.atoms_same_charge[nb_index]:
                                nb_names.append(atom)
                        if chm_gp.cooperators!=None:
                            for coop_index in chm_gp.cooperators[i]:
                                for atom in chm_gp.atoms_same_charge[coop_index]:
                                    coop_names.append(atom)
                    # Now we have the name lists, create the 'gtcnp' for charge:
                    roc=[]; ron=[]
                    # I think this decision will influence the speed of optimization,
                    # think twice before using ...
                    for n in range(len(nb_names)):
                        ron.append(round(-len(center_names)/len(nb_names),5))
                    #for r in range(len(roc)):
                    #    roc.append(1)
                    add_a_new_group(gs, gtcnp(par_type="charge", center_names=center_names, 
                                              cooperators=coop_names, neighbors=nb_names,
                                              original_p=chm_gp.charges[i], 
                                              targeted_range=[round(chm_gp.charges[i] - 0.2, 2), 
                                                              round(chm_gp.charges[i] + 0.2, 2)],
                                              roc=roc, ron=ron)
                                   )
                else:
                    self.level_print("Skipping charge(s) for {} ...".format(chm_gp.atoms[i]), 2, print_level)
            self.level_print("", 1, print_level)
        self.level_print("Total {} gtcnps created for {}\n".format(len(gs), self.lipname), 1, print_level)
        return gs


