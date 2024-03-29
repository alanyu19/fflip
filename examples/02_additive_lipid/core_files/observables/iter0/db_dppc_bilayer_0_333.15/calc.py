#!/usr/bin/env python

# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=ivy
#SBATCH --ntasks=1
#SBATCH -J thickness


from rflow.observables import *
from rflow.trajectory import *
from fflip.omm.util import get_md_options as gmd
from fflip.analysis.edp import *


option_file_name = "obscalc.inp"

opts = gmd(option_file_name)
blk_size = int(opts['block_size'])  # number of traj files each block
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
nlipids = int(opts['nlip'])

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


db_evaluator = MembraneThickness(nbins=500, edge_bins=130, box_length_fixed=10)

to_calc = TimeSeries(
    evaluator=db_evaluator, filename="./block_data/db_{}.dat".format(
        starting_traj
    )
)

for traj in trajs:
    to_calc(traj)

