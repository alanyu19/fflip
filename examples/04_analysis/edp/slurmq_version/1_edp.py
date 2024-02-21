#!/usr/bin/env python
#SBATCH --output=./log/trj_%a.out
#SBATCH --error=./log/trj_%a.err
#SBATCH --ntasks=1
#SBTACH --exclusive
#SBATCH --cpus-per-task=20
#SBATCH --array=101-105
#SBATCH --job-name=edp
#SBATCH -p norm
#SBATCH -t 12:00:00


import os
import time
import warnings
import mdtraj as md
import numpy as np
import os
from rflow.trajectory import TrajectoryIterator, center_of_mass_of_selection
from rflow.openmm.app import CharmmPsfFile
from fflip.analysis.recenter import TrajectoryRecenter
from fflip.analysis.edp import *


## update system info ##
PSF_FILE = "../../sim/plant.psf"
TRAJ_TEMPLATE = "../../sim/trj/dyn{}.dcd"
## end of update ##

## NO update needed below ##
FIRST = int(os.environ['SLURM_ARRAY_TASK_ID'])
LAST = FIRST
recentered_trajs_dir = "recenter"
if os.path.isdir(recentered_trajs_dir):
    warnings.warn(f"Directory {recentered_trajs_dir} exists, overwriting it...")
else:
    os.mkdir(recentered_trajs_dir)
output_template = os.path.join(recentered_trajs_dir, "dyn{}-recenter.dcd")
traj_recenter = TrajectoryRecenter(psf_file=PSF_FILE, input_template=TRAJ_TEMPLATE, first=FIRST, last=LAST)
traj_recenter.recenter_all_trajectories(output_format=output_template, lipid_head_atoms='name C2', lipid_tail_atoms='name C212')

edp_factory = ElectronDensityFactory(PSF_FILE, output_template, com_selection='not water')
edp_factory(FIRST, LAST)
