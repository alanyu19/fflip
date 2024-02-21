from rflow.observables import *
from rflow.trajectory import *
from fflip.omm.util import get_md_options as gmd
from fflip.analysis.recenter import TrajectoryRecenter
from fflip.analysis.edp import *


## update system info ##
PSF_FILE = "../sim/plant.psf"
TRAJ_TEMPLATE = "../sim/trj/dyn{}.dcd"
FIRST = 101
LAST = 101
## end of update ##

## special updates ##
EDGE_BINS=130
## end of special upates ##

## no updates needed below ##
recentered_trajs_dir = "recenter"
if os.path.isdir(recentered_trajs_dir):
    raise Exception(f"Directory {recentered_trajs_dir} exists, please check and delete/rename it to continue.")
os.mkdir(recentered_trajs_dir)
output_template = os.path.join(recentered_trajs_dir, "dyn{}-recenter.dcd")
traj_recenter = TrajectoryRecenter(psf_file=PSF_FILE, input_template=TRAJ_TEMPLATE, first=FIRST, last=LAST)
traj_recenter.recenter_all_trajectories(output_format=output_template, lipid_head_atoms='name C2', lipid_tail_atoms='name C212')

trajs = TrajectoryIterator(
    first_sequence=FIRST,
    last_sequence=LAST,
    filename_template=output_template,
    topology_file=PSF_FILE,
    atom_selection="all",
    load_function=md.load_dcd
)

db_evaluator = MembraneThickness(nbins=500, edge_bins=EDGE_BINS, box_length_fixed=10)
db_timeseries = TimeSeries(
    evaluator=db_evaluator, filename="db_{}_{}.dat".format(FIRST, LAST)
)

for traj in trajs:
    db_timeseries(traj)
