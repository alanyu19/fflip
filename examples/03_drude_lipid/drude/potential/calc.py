#!/usr/bin/env python
# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --partition=ivy,sbr,hwell,spodall,hpodall
#SBATCH --ntasks=1

# import os
import sys
import time
import fflip.omm.energiesforces as ef

from simtk.openmm.app import LJPME

from rflow.observables import *
from rflow.trajectory import *
from fflip.omm.util import *
from fflip.omm.util import get_md_options as gmd
from fflip.drude import *
from fflip.ljpme.calcfuncs import save_offsets, get_one_group_with_offset


# ----------------------------- Reading Inputs -----------------------------

sys.path.append(os.getcwd())

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

lip = parse_lipid(lipname)

index = int(os.environ["SLURM_ARRAY_TASK_ID"])

assert 0 <= index <= len(lip.parse_groups()), \
    "the index of parameter group is out of range (ignored), " \
    "please check your slurmq array indexes :)"

# --------------------------------------------------------------------------

# ---------------------------- Record the Offsets --------------------------
if starting_traj == first_dcd and index == 0:
    parameter_groups = lip.parse_groups(id_allowed='all')
    offsets = [
        gen_sensitivity_offset(ps, percentage=percentage)
        for ps in parameter_groups
    ]
    save_offsets(parameter_groups, offsets, "./offsets.log")
# --------------------------------------------------------------------------

parameter_files = glob.glob("/u/alanyu/top_yalun/for_drude/2013/*.str")

psf, topology, parameters = read_structure_parameter_files(
    psf_file, parameter_files
)

if os.path.isfile('solution.txt'):
    # should do something to the LJ
    time.sleep(2)
    pass

if index > 0:
    # currently no forbidden ids
    pgroup, offset = get_one_group_with_offset(
        index, lip, percentage, id_allowed_='all'
    )
    # TODO: finish this part
    if pgroup.partype in ['sigma', 'epsilon']:
        pass

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
    parameter_sets = lip.parse_gtcnp()
    all_offsets = [
        gen_sensitivity_offset(ps, percentage=sol[i])
        for i, ps in enumerate(parameter_sets)
    ]
    for pset, offset in zip(parameter_sets, all_offsets):
        # the LJ parameters should be skipped in the
        # change_drude_ff_parameters
        change_drude_ff_parameters(
            system, topology, pset, offset
        )
    # change_nb_exceptions(psf, system), might not be needed in DRUDE
    # from torfix import *
    # system = BrutalTorsionParameters(system, psf, tfixes)

trajs = TrajectoryIterator(
    first_sequence=starting_traj,
    last_sequence=starting_traj + blk_size - 1,
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd)

if index == 0:
    energy_evaluator = ef.DrudeEnergy(
        system, psf, paragroups=[], paraoffsets=[],
        use_new_method=False, use_platform='CPU'
    )
    to_calc = TimeSeries(
        evaluator=energy_evaluator,
        filename="./block_data/original_{}.dat".format(starting_traj)
    )

elif index > 0:
    energy_evaluator = ef.DrudeEnergy(
        system, psf, paragroups=pgroup, paraoffsets=offset,
        use_new_method=False
    )
    to_calc = TimeSeries(
        evaluator=energy_evaluator,
        filename="./block_data/perturbed_{}_{}_{}.dat".format(
            starting_traj, pgroup.cgid, pgroup.internal_id
        )
    )
    # TODO: parameter set name would be better than just index?

for traj in trajs:
    to_calc(traj)
