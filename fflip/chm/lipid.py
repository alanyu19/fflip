import os
import sys
#import omm_vfswitch
#import copy

import numpy as np
import simtk.unit as u
from simtk.openmm import Platform, MonteCarloMembraneBarostat, DrudeLangevinIntegrator
from simtk.openmm import Context, LangevinIntegrator, NonbondedForce, CustomNonbondedForce
from simtk.openmm.app import ForceField, Simulation
from simtk.openmm.app import CharmmPsfFile, CharmmParameterSet, CharmmCrdFile
from simtk.openmm.app import LJPME, PME, HBonds
import mdtraj as md
from openmmtools.testsystems import LennardJonesFluid
import scipy.optimize

from coffe.omm.playpara import *
from coffe.omm.paragroup import *


def add_a_new_group(existing_groups, new_group):
    # make sure that there is no repeated definition
    for g in existing_groups:
        if g.par_type == new_group.par_type and g.center_names == new_group.center_names:
            raise GroupRepeatedError(new_group)
    existing_groups.append(new_group)


class charmm_group():
    def __init__(self, **kwargs):
        self.num_atom_category = kwargs["num_atom_category"]
        self.atoms = kwargs["atoms"] # should be list inside list structure
        self.charges = kwargs["charges"]
        self.half_r_mins = kwargs["half_r_mins"]
        self.epsilons = kwargs["epsilons"]
        self.add_lj_gtcnp = kwargs["add_lj_gtcnp"]
        self.atoms_same_lj = kwargs["atoms_same_lj"]
        self.add_charge_gtcnp = kwargs["add_charge_gtcnp"]
        self.neighbors = kwargs["neighbors"]
        self.cooperators = kwargs["cooperators"]
        self.atoms_same_charge = kwargs["atoms_same_charge"]
        if "exclusion" in kwargs:
            self.exclusion = kwargs["exclusion"]
    

class lipid():

    def __init__(self, **kwargs):
        assert "charmm_group_list" in kwargs
        self.charmm_group_list = kwargs["charmm_group_list"]
        self.num_charmm_groups = len(self.charmm_group_list)
        assert "lipname" in kwargs
        self.lipname= kwargs["lipname"]
        self.cgroups = []
        for group in self.charmm_group_list:
            self.cgroups.append(group)
            

    def parse_gtcnp(self):    
        gs = []
        for counter, chm_gp in enumerate(self.cgroups):
            print("Creating 'gtcnp's for {} group {} ... ".format(self.lipname, counter+1))
            for i in range(chm_gp.num_atom_category):
                # LJ
                if chm_gp.add_lj_gtcnp[i]:
                    atom_list=[]
                    for atom in chm_gp.atoms[i]:
                        atom_list.append(atom)
                    for atom in chm_gp.atoms_same_lj[i]:
                        atom_list.append(atom)
                    add_a_new_group(gs, gtcnp(par_type = "sigma", center_names=atom_list, 
                                              original_p = round((1/2)**(1/6) * chm_gp.half_r_mins[i] * 0.2, 5), 
                                              targeted_range=[0.8, 1.2]
                                             )
                                   )
                    add_a_new_group(gs, gtcnp(par_type = "epsilon", center_names=atom_list,
                                              original_p=round(chm_gp.epsilons[i] * 4.184, 4), 
                                              targeted_range=[0.8, 1.2]
                                             )
                                   )
                else:
                    print("Skipping LJ parameters for {} ...".format(chm_gp.atoms[i]))
                #
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
                    print("Skipping charge(s) for {} ...".format(chm_gp.atoms[i]))
            print("")
        print("Total {} gtcnps created for {}\n".format(len(gs), self.lipname))
        return gs

# ***************************************************************************************************************
