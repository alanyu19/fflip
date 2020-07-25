# -*- coding: utf-8 -*-

import os
import sys
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
from coffe.omm.playpara import *
from coffe.omm.paragroup import *
from coffe.omm.util import *
from rflow.observables import *
from rflow.trajectory import *
from fflip.chm import *

parser = argparse.ArgumentParser()

parser.add_argument('--psf', type = str, help = "psf file")
parser.add_argument('-t', '--trajectory_template', type = str,
                    help = "template of trajectory files")
parser.add_argument('-p', '--percentage', type = float, default = 1, 
                    help = "perturbed percentage of parameter (for charges, 1% is assumed to be 0.004)")
parser.add_argument('-s', '--starting_trajectory', type = int, help = "the index of the first trajectory in a block")
args = parser.parse_args()
percentage = args.percentage
trj_template = args.trajectory_template
starting_traj = args.starting_trajectory
psf_file = args.psf

with open('lipname.txt', 'r') as f:
    lip = f.readline().strip()

lip = parse_lipid(lip)

index = int(os.environ["SLURM_ARRAY_TASK_ID"])

# keep changes ~ 1% * percentage of the average (over all the atoms in the lipid)
# of the nonbonded type (charge/epsilon/sigma)

def save_nb_change_log(gs, offsets, file_template = "./gtcnp_order.txt"):
    if os.path.isfile(file_template):
        os.system("rm {}".format(file_template))
    for g, off in zip(gs, offsets):
        with open(file_template, 'a') as f:
            f.write("{0:>5} {1:>9} {2:>10.4f}\n".format(g.center_names[0], g.par_type, off))


def lets_get_it(index, lip):
    parameter_sets=lip.parse_gtcnp()
    all_offsets = [gen_sensitivity_offset(ps, percentage = percentage) for ps in parameter_sets]
    pset = [parameter_sets[index-1]]  # keep this list format for the other function/class (parameterEnergy?)
    offset = [gen_sensitivity_offset(parameter_sets[index-1], percentage = percentage)]  # same here
    if index == 1: #TODO: and starting_traj == 31 (only save for the first trajectory is enough
        save_nb_change_log(parameter_sets, all_offsets, "./groups_{}.nbg".format(starting_traj))
    return pset, offset

number_of_paragroups = len(lip.parse_gtcnp())

assert 0 <= index <= number_of_paragroups, "the index of parameter group is out of range (ignored), \
        please check your slurmq array indexes :)"

if index > 0:
    gss, offsets = lets_get_it(index, lip)

# TODO
if lip.lipname.lower() == 'prpc':
    parameter_files = glob.glob("/u/alanyu/top_yalun/for_ljpme/prpc/original/*")
else:
    parameter_files = glob.glob("/u/alanyu/top_yalun/for_ljpme/original/*")

psf, topology, parameters = ReadStructureParameterFiles(psf_file, parameter_files)

# The following traj is only used to set the unitcell lengths
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

if os.path.isfile('solution.txt'):
    sol = filter_solution('solution.txt')
    parameter_sets=lip.parse_gtcnp()
    all_offsets = [gen_sensitivity_offset(ps, percentage = sol[i]) for i, ps in enumerate(parameter_sets)]
    for g, offset in zip(parameter_sets, all_offsets):
        system, oldp, newp = BrutalNonbondedParameter(system, topology, g, offset)


if index == 0:
    energy_evaluator = ef.parameterEnergy(system, psf, paragroups = [], paraoffsets = [], use_new_method=False, use_platform='CPU')
    to_calc = TimeSeries(evaluator=energy_evaluator, filename = "./block_data/original_{}.dat".format(starting_traj))
else:
    energy_evaluator = ef.parameterEnergy(system, psf, paragroups = gss, paraoffsets = offsets, use_new_method=False)
    to_calc = TimeSeries(evaluator=energy_evaluator, filename = "./block_data/perturbed_{}_{}.dat".format(starting_traj, index)) 
    # TODO: parameter set name would be better than just index


for traj in trajs:
    to_calc(traj)

