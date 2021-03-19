#!/usr/bin/env python

# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=ivy,sbr,hwell
#SBATCH --ntasks=1
#SBATCH -J thickness


from coffe.omm.util import get_md_options as gmd

from rflow.edp import *
from rflow.observables import *
from rflow.trajectory import *


option_file_name = "obscalc.inp"

opts = gmd(option_file_name)
blk_size = int(opts['block_size'])  # number of traj files each block
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
nlipids = int(opts['nlip'])
lipname = str(opts['lipname'])

starting_traj = int(os.environ['SLURM_ARRAY_TASK_ID'])

psf = CharmmPsfFile(psf_file)

trajs = TrajectoryIterator(
    first_sequence=starting_traj,
    last_sequence=starting_traj+blk_size-1,
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd
)

dhh_evaluator = MembraneDhh(nbins=500, box_length_fixed=10)

to_calc = TimeSeries(
    evaluator=dhh_evaluator, filename="./block_data/dhh_{}.dat".format(
        starting_traj
    )
)

for traj in trajs:
    to_calc(traj)

