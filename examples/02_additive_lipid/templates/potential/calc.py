#!/usr/bin/env python
# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=4G

import os
import sys
import time
sys.path.append(os.getcwd())

from simtk.openmm.app import LJPME

import fflip.omm.energiesforces as ef
from fflip.omm.playpara import *
from fflip.omm.util import *
from fflip.omm.util import get_md_options as gmd
from rflow.observables import *
from rflow.trajectory import *
from fflip.chm import *
from def_lipids import *

option_file_name = "potcalc.inp"

opts = gmd(option_file_name)

blk_size = int(opts['block_size'])
percentage = float(opts['percentage'])
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
lipname = str(opts['lipname'])
first_dcd = int(opts['first_trj'])
solution = opts['solution']
torfix = opts['torfix']

starting_traj = int(sys.argv[1])

lip = match_lipid[lipname]

index = int(os.environ["SLURM_ARRAY_TASK_ID"])


def save_nb_change_log(gs, offsets, file_template = "./nbgroups_order.txt"):
    if os.path.isfile(file_template):
        os.system("rm {}".format(file_template))
    for g, off in zip(gs, offsets):
        with open(file_template, 'a') as f:
            f.write(
                "{0:>5} {1:>9} {2:>10.4f}\n".format(
                    g.center_names[0], g.par_type, off
                )
            )


def get_parameter_set_and_offset_by_index(index, lip):
    parameter_sets = lip.parse_nbgroups()
    pset = [parameter_sets[index-1]]
    # keep this list format for the other function/class (parameterEnergy?)
    offset = [
        gen_sensitivity_offset(parameter_sets[index-1], percentage=percentage)
    ]  # same here
    return pset, offset


if starting_traj == first_dcd and index == 0:
    parameter_sets = lip.parse_nbgroups()
    all_offsets = [
        gen_sensitivity_offset(ps, percentage=percentage)
        for ps in parameter_sets
    ]
    save_nb_change_log(parameter_sets, all_offsets, "./param_change.txt")

assert 0 <= index <= len(lip.parse_nbgroups()), \
    "the index of parameter group is out of range (ignored), " \
    "please check your slurmq array indexes :)"

if index > 0:
    gss, offsets = get_parameter_set_and_offset_by_index(index, lip)

# Update the toppar
parameter_files = glob.glob("./toppar/*")

psf, topology, parameters = read_structure_parameter_files(
    psf_file, parameter_files
)

# The following traj is only used to set the unitcell lengths
traj = md.load_dcd(traj_template.format(starting_traj), top=topology)

psf.setBox(*traj.unitcell_lengths[0])

system = psf.createSystem(
    parameters, 
    nonbondedMethod=LJPME,
    constraints=HBonds, 
    nonbondedCutoff=1.0 * u.nanometer,
    ewaldErrorTolerance=0.0001
)

trajs = TrajectoryIterator(
    first_sequence=starting_traj,
    last_sequence=starting_traj + blk_size - 1,
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd)

if os.path.isfile('solution.txt'):
    time.sleep(2)
    sol = filter_solution('solution.txt')
    parameter_sets = lip.parse_nbgroups()
    all_offsets = [
        gen_sensitivity_offset(ps, percentage=sol[i])
        for i, ps in enumerate(parameter_sets)
    ]
    for g, offset in zip(parameter_sets, all_offsets):
        system, oldp, newp = BrutalNonbondedParameter(
            system, topology, g, offset
        )
    change_nb_exceptions(psf, system)
    from torfix import *
    system = BrutalTorsionParameters(system, psf, tfixes)

if index == 0:
    energy_evaluator = ef.ParameterEnergy(
        system, psf, paragroups=[], paraoffsets=[],
        # Use 'CPU' if not running this on a gpu platform
        # If using CUDA, load the version matching your openmm
        # before calling fflip obspot
        use_new_method=False, use_platform='CPU'
    )
    to_calc = TimeSeries(
        evaluator=energy_evaluator,
        filename="./block_data/original_{}.dat".format(starting_traj)
    )
else:
    energy_evaluator = ef.ParameterEnergy(
        system, psf, paragroups=gss, paraoffsets=offsets,
        use_new_method=False, use_platform='CPU'
    )
    to_calc = TimeSeries(
        evaluator=energy_evaluator,
        filename="./block_data/perturbed_{}_{}.dat".format(starting_traj, index)
    )

for traj in trajs:
    to_calc(traj)

