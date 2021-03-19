# -*- coding: utf-8 -*-

import numpy as np
import argparse
from simtk.openmm.app import CharmmPsfFile #CharmmParameterSet, CharmmCrdFile, PDBFile
import mdtraj as md
from rflow.trajectory import *
from fflip.analysis.rdf import *

parser = argparse.ArgumentParser()

parser.add_argument('--psf', type = str, help = "psf file")
parser.add_argument('-t', '--trajectory_template', type = str,
                    help = "template of trajectory files")
parser.add_argument('-s', '--starting_trajectory', type = int, help = "the index of the first trajectory in a block")
args = parser.parse_args()
trj_template = args.trajectory_template
starting_traj = args.starting_trajectory
psf_file = args.psf

index = int(os.environ["SLURM_ARRAY_TASK_ID"])

psf = CharmmPsfFile(psf_file)

phos_pairs = ["Os1-HW", "O2-HW", "Os1-OW", "O2-OW"]
phos_selections = [["name O11 or name O12", "name H1 or name H2"], ["name O13 or name O14", "name H1 or name H2"],
                   ["name O11 or name O12", "name OH2"], ["name O13 or name O14", "name OH2"]]
es_pairs = ["Ob-OW", "Os2-OW", "Ob-HW", "Os2-HW"]
es_selections = [["name O22 or name O32", "name OH2"], ["name O21 or name O31", "name OH2"],
                 ["name O22 or name O32", "name H1 or name H2"], ["name O21 or name O31", "name H1 or name H2"]]

pairs = phos_pairs + es_pairs
selections = phos_selections + es_selections

def get_rdf_mdtraj(name, atom_selection, psf, trajs, start_traj = starting_traj):
    rdf_mdtraj = rdf(psf, name, atom_selection, r_range = [0, 1.0], dimension = 3, method = 'mdtraj')
    for trj in trajs:
        rdf_mdtraj(trj, start_traj, save_sparse = True, sparse_bin_width = 0.003, sparse_subfolder = './block_data')
    rdf_mdtraj.save_to_file()

trajs_ljpme = CharmmTrajectoryIterator(
    first_sequence = starting_traj, last_sequence = starting_traj, # this is currently taking only one trj
    filename_template = os.path.join(trj_template, 'dyn{}.dcd'), 
    topology_file = psf_file, 
    atom_selection = "all",
    load_function=md.load_dcd)

# the following index is succeeded from the array job submission script 'rdf_calc.sh'
name = pairs[index-1]
atom_selection = selections[index-1]
get_rdf_mdtraj(name, atom_selection, psf, trajs_ljpme)

