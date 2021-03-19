# -*- coding: utf-8 -*-

import os
import sys
import math
import glob
import numpy as np
import argparse

import simtk.unit as u
from simtk.openmm import Platform, MonteCarloMembraneBarostat, DrudeLangevinIntegrator
from simtk.openmm import Context, LangevinIntegrator, NonbondedForce, CustomNonbondedForce
from simtk.openmm.app import ForceField, Simulation
from simtk.openmm.app import CharmmPsfFile, CharmmParameterSet, CharmmCrdFile, PDBFile
from simtk.openmm.app import LJPME, PME, HBonds
from openmmtools.testsystems import LennardJonesFluid
import mdtraj as md

import coffe.omm.energiesforces as ef
from coffe.omm.util import *
from rflow.observables import *
from rflow.trajectory import *

parser = argparse.ArgumentParser()

parser.add_argument('--psf', type = str, help = "psf file")
parser.add_argument('-t', '--trajectory_template', type = str,
                    help = "template of trajectory files")
parser.add_argument('-s', '--starting_trajectory', type = int, help = "the index of the first trajectory in a block")
args = parser.parse_args()
trj_template = args.trajectory_template

starting_traj = args.starting_trajectory

psf_file = args.psf
parameter_files = ["/u/alanyu/coffe/tests/omm/data/par_all36_lipid.prm",
                   "/u/alanyu/coffe/tests/omm/data/top_all36_lipid.rtf",
                   "/u/alanyu/coffe/tests/omm/data/toppar_water_ions.str",
                   "/u/alanyu/top_yalun/for_ua/c36ua.str" ]
psf, topology, parameters = ReadStructureParameterFiles(psf_file, parameter_files)

# The following traj is only used to set the unitcell lengths
print(trj_template)
traj = md.load_dcd(trj_template + "/dyn{}.dcd".format(starting_traj), top=topology)
psf.setBox(*traj.unitcell_lengths[0])
system = psf.createSystem(
    parameters, 
    nonbondedMethod=LJPME,
    constraints=HBonds, 
    nonbondedCutoff=1.0 * u.nanometer,
    ewaldErrorTolerance=0.0005
)

trajs = CharmmTrajectoryIterator(
    first_sequence = starting_traj,
    last_sequence=starting_traj + 9,
    filename_template = os.path.join(trj_template, 'dyn{}.dcd'),
    topology_file = psf_file,
    atom_selection = "all",
    load_function=md.load_dcd)

# this number of lipids should not be constant in the future
sa_evaluator = ef.AreaPerLipid(36)
to_calc = TimeSeries(evaluator=sa_evaluator, filename = "./block_data/area_{}.dat".format(starting_traj))

for traj in trajs:
    to_calc(traj)

