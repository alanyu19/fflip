#!/usr/bin/env python
# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A.out
#SBATCH --error=./log/master_%A.err
#SBATCH --time=12:00:00
#SBATCH --partition=v100,pascal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB

import os
import sys
import time
sys.path.append(os.getcwd())

from simtk.openmm.app import LJPME, HBonds

import fflip.omm.energiesforces as ef
from rflow.observables import *
from rflow.trajectory import *
from fflip.omm.util import *
from fflip.omm.util import get_md_options as gmd
from fflip.chm import *
from def_lipids import *
from fflip.ljpme.param import gen_param_offset, save_offsets, \
    get_one_group_with_offset


# ----------------------------- Reading Inputs -----------------------------

option_file_name = "potcalc.inp"

opts = gmd(option_file_name)

blk_size = int(opts['block_size'])
amount = float(opts['amount'])
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
lipid_name = str(opts['name'])
first_dcd = int(opts['first_trj'])
solution = opts['solution']
torfix = opts['torfix']
toppar_path = opts['toppar_path']

starting_traj = int(sys.argv[1])

lipid = match_lipid[lipid_name]

# ---------------------------- Record the Offsets --------------------------
if starting_traj == first_dcd:
    parameter_groups = lipid.parse_groups(id_allowed='all')
    offsets = [
        gen_param_offset(ps, amount=amount)
        for ps in parameter_groups
    ]
    save_offsets(parameter_groups, offsets, "./offsets.log")
# --------------------------------------------------------------------------

parameter_files = glob.glob(os.path.join(toppar_path, "*"))

psf, topology, parameters = read_structure_parameter_files(
    psf_file, parameter_files
)

# --------------------------- Set up the System ----------------------------
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

if os.path.isfile('solution.txt'):
    time.sleep(2)
    sol = filter_solution('solution.txt')
    parameter_sets = lip.parse_groups()
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

trajs = TrajectoryIterator(
    first_sequence=starting_traj,
    last_sequence=starting_traj + blk_size - 1,
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd)

for index in range(len(lipid.parse_groups()) + 1): 
    if index > 0:
        # currently no forbidden ids
        pgroup, offset = get_one_group_with_offset(
            index, lipid, amount, id_allowed='all'
        )
    if index == 0:
        energy_evaluator = ef.ParameterEnergy(
            system, psf, paragroups=[], paraoffsets=[],
            # Use 'CPU' if not running this on a gpu platform
            # If using CUDA, load the version matching your openmm
            # before calling fflip obspot
            use_new_method=False, use_platform='CUDA'
        )
        to_calc = TimeSeries(
            evaluator=energy_evaluator,
            filename="./block_data/original_{}.dat".format(starting_traj)
        )
    elif index > 0:
        energy_evaluator = ef.ParameterEnergy(
            system, psf, paragroups=pgroup, paraoffsets=offset,
            use_new_method=False, use_platform='CUDA'
        )
        to_calc = TimeSeries(
            evaluator=energy_evaluator,
            filename="./block_data/perturbed_{}_{}.dat".format(starting_traj, index)
        )

    for traj in trajs:
        to_calc(traj)

