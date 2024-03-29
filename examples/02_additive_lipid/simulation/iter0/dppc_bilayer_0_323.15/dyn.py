#! /usr/bin/env python

import os
import sys
import warnings
import glob
import numpy as np

import simtk.unit as u
from simtk.openmm.app import LJPME
from simtk.openmm import MonteCarloMembraneBarostat, MonteCarloBarostat, LangevinIntegrator

from rflow.integrators import NoseHooverChainVelocityVerletIntegrator
from rflow import RickFlow, FreeEnergyCosineSeries, RelativePartialCenterOfMassRestraint
from fflip.omm.playpara import *
from fflip.omm.paragroup import *
from fflip.ljpme.param import *
from fflip.omm.util import *
from fflip.chm import *

from def_lipids import *


options = get_md_options('input.mdo')
mdo = parse_md_options(options)

psf_file = mdo['psf'] 
crd_file = mdo['crd']
surface_tension = mdo['surface_tension'] 
xvec = mdo['xvec']
yvec = mdo['yvec']
zvec = mdo['zvec']
lipname = mdo['lipid']
lipid = match_lipid[lipname]
temp = mdo['temperature'] 
baro = mdo['barostat']
intgrt = mdo['integrator']
zmode = mdo['zmode']
change_para = mdo['change_param']
if change_para.lower() == 'yes':
    try:
        os.system("cp {} solution.txt".format(mdo['sfile']))
        os.system("cp {} torfix.py".format(mdo['tfile']))
    except:
        raise Exception("copy solution file and torfix file failed")

#
#  ========= Setup System ==========
#

# must update!
# keep it as original paramters in all iterations!
top_files = glob.glob("./toppar/*")

workflow = RickFlow(
    toppar=top_files,
    psf=psf_file,
    crd=crd_file,
    box_dimensions=[xvec, yvec, zvec], 
    gpu_id=0,
    cutoff_distance=10*u.angstrom,
    vdw_switching="openmm",
    nonbonded_method=LJPME,
    center_around='not water',
    dcd_output_interval=500,
    table_output_interval=500,
    steps_per_sequence=500000,
    misc_psf_create_system_kwargs={"ewaldErrorTolerance": 0.0001}
)

temperature = temp * u.kelvin

#
#  ========= Define Integrator and Barostat ==========
#

if intgrt == 'L': 
    integrator = LangevinIntegrator(
        temperature, 1/u.picosecond, 0.002*u.picosecond
    )
else:
    integrator = NoseHooverChainVelocityVerletIntegrator(
        workflow.system, temperature, 50.0 / u.picosecond, 2.0 * u.femtosecond,
        chain_length=3, num_mts=3, num_yoshidasuzuki=3
    )

if baro == 'MCM':
    barostat = MonteCarloMembraneBarostat(
        1.0 * u.atmosphere, surface_tension, 
        integrator.getTemperature(), 
        0, zmode, 50
    )

elif baro == 'MC':
    barostat = MonteCarloBarostat(
        1.0 * u.atmosphere, 
        integrator.getTemperature()
    )

#
#  ========= Add parameter offsets =========
#
if change_para.lower() == 'yes':
    # don't change this line, solution.txt is copied over by the program
    sol = filter_solution('solution.txt')
    parameter_sets = lipid.parse_groups()
    all_offsets = [gen_param_offset(ps, amount=sol[i]) for i, ps in enumerate(parameter_sets)]
    topology = md.Topology.from_openmm(workflow.psf.topology)
    for g, offset in zip(parameter_sets, all_offsets):
        workflow._system, oldp, newp = BrutalNonbondedParameter(workflow.system, topology, g, offset)
    change_nb_exceptions(workflow.psf, workflow.system)

if change_para.lower() == 'yes':
    from torfix import *
    workflow._system = BrutalTorsionParameters(workflow.system, workflow.psf, tfixes)

workflow.prepareSimulation(integrator, barostat)

#
#  ========= Run Simulation  ==========
#

if workflow.next_seqno == 1:
    workflow.simulation.minimizeEnergy()
if __name__ == "__main__":
    workflow.run()

