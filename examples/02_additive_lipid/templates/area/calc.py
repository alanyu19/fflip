#!/usr/bin/env python

# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --partition=ivy
#SBATCH --ntasks=1
#SBATCH -J area

from simtk.openmm.app import LJPME

from fflip.omm.util import get_md_options as gmd

from fflip.omm.util import *
from rflow.observables import *
from rflow.trajectory import *


option_file_name = "obscalc.inp"

opts = gmd(option_file_name)
blk_size = int(opts['block_size'])  # number of traj files each block
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
nlipids = int(opts['nlip'])
lipname = str(opts['name'])

starting_traj = int(os.environ['SLURM_ARRAY_TASK_ID'])


parameter_files = glob.glob("./toppar/*")

psf, topology, parameters = read_structure_parameter_files(
    psf_file, parameter_files
)


# The following traj is only used to set the unitcell lengths
traj = md.load_dcd(
    traj_template.format(starting_traj), top=topology
)
psf.setBox(*traj.unitcell_lengths[0])
system = psf.createSystem(
    parameters, 
    nonbondedMethod=LJPME,
    constraints=HBonds, 
    nonbondedCutoff=1.0 * u.nanometer,
    ewaldErrorTolerance=0.0005
)


trajs = TrajectoryIterator(
    first_sequence=starting_traj,
    last_sequence=starting_traj + blk_size - 1,
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd
)


# this number of lipids should not be constant in the future
sa_evaluator = AreaPerLipid(int(nlipids))
to_calc = TimeSeries(
    evaluator=sa_evaluator, filename="./block_data/area_{}.dat".format(
        starting_traj
    )
)

for traj in trajs:
    to_calc(traj)

