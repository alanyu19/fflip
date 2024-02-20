import os
import time
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
FIRST = 101
LAST = 101
## end of update ##


## NO update needed below ##
recentered_trajs_dir = "recenter"
if os.path.isdir(recentered_trajs_dir):
    raise Exception(f"Directory {recentered_trajs_dir} exists, please check and delete/rename it to continue.")
os.mkdir(recentered_trajs_dir)
output_template = os.path.join(recentered_trajs_dir, "dyn{}-recenter.dcd")
traj_recenter = TrajectoryRecenter(psf_file=PSF_FILE, input_template=TRAJ_TEMPLATE, first=FIRST, last=LAST)
traj_recenter.recenter_all_trajectories(output_format=output_template, lipid_head_atoms='name C2', lipid_tail_atoms='name C212')

edp_factory = ElectronDensityFactory(PSF_FILE, output_template, com_selection='not water')
edp_factory(FIRST, LAST)
