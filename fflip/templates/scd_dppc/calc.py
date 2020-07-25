# -*- coding: utf-8 -*-

import argparse
import numpy as np

import simtk.unit as u
from simtk.openmm import Platform, MonteCarloMembraneBarostat, DrudeLangevinIntegrator
from simtk.openmm import Context, LangevinIntegrator, NonbondedForce, CustomNonbondedForce
from simtk.openmm.app import ForceField, Simulation
from simtk.openmm.app import CharmmPsfFile, CharmmParameterSet, CharmmCrdFile, PDBFile
from simtk.openmm.app import LJPME, PME, HBonds
import mdtraj as md
from openmmtools.testsystems import LennardJonesFluid

import coffe.omm.energiesforces as ef
from coffe.omm.playpara import *
from coffe.omm.paragroup import *
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
                   "/u/alanyu/coffe/tests/omm/data/toppar_water_ions.str"]

psf, topology, parameters = ReadStructureParameterFiles(psf_file, parameter_files)

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


scd_evaluators = []
for i in range(3, 16):
    scd_evaluators.append(ef.order_parameter(topology, 'C2{}'.format(i), ['H{}R'.format(i), 'H{}S'.format(i)]))
    scd_evaluators.append(ef.order_parameter(topology, 'C3{}'.format(i), ['H{}X'.format(i), 'H{}Y'.format(i)]))

scd_evaluators_extra = []
scd_evaluators_extra.append(ef.order_parameter(topology, 'C22', 'H2R'))
scd_evaluators_extra.append(ef.order_parameter(topology, 'C22', 'H2S'))
scd_evaluators.append(ef.order_parameter(topology, 'C32', ['H2X', 'H2Y']))

# Headgroup (compare to Toward Atomistic Resolution Structure of Phosphoatidylcholine Headgroup and Glycerol Backbone at Different Ambient Conditions. JPCB 2015, 119, 15075-15088)
scd_evaluators_head = []
scd_evaluators_head.append(ef.order_parameter(topology, 'C11', ['H11A', 'H11B']))
scd_evaluators_head.append(ef.order_parameter(topology, 'C12', ['H12A', 'H12B']))
# Glycerol Backbone
scd_evaluators_gly = []
scd_evaluators_gly.append(ef.order_parameter(topology, 'C1', 'HA'))
scd_evaluators_gly.append(ef.order_parameter(topology, 'C1', 'HB'))
scd_evaluators_gly.append(ef.order_parameter(topology, 'C2', 'HS'))
scd_evaluators_gly.append(ef.order_parameter(topology, 'C3', 'HX'))
scd_evaluators_gly.append(ef.order_parameter(topology, 'C3', 'HY'))

scd = []
for scd_evaluator in scd_evaluators:
    scd.append(TimeSeries(evaluator=scd_evaluator,
        filename = "./block_data/scd-{}_{}.dat".format(scd_evaluator.atom1.lower(), starting_traj))
        )

for scd_evaluator in scd_evaluators_extra:
    scd.append(TimeSeries(evaluator=scd_evaluator, 
        filename = "./block_data/scd-{}-{}_{}.dat".format(scd_evaluator.atom1.lower(), scd_evaluator.atom2.lower(), starting_traj))
        )

for scd_evaluator in scd_evaluators_head:
    scd.append(TimeSeries(evaluator=scd_evaluator,
        filename = "./block_data/scd-{}_{}.dat".format(scd_evaluator.atom1.lower(), starting_traj))
        )

for scd_evaluator in scd_evaluators_gly:
    scd.append(TimeSeries(evaluator=scd_evaluator, 
        filename = "./block_data/scd-{}-{}_{}.dat".format(scd_evaluator.atom1.lower(), scd_evaluator.atom2.lower(), starting_traj))
        )

for traj in trajs:
    for s in scd:
        s(traj)

