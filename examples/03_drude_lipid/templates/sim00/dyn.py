#! /usr/bin/env python

import os
import sys
import warnings
import glob
import numpy as np

import simtk.unit as u
from simtk.openmm.app import LJPME
from simtk.openmm import MonteCarloMembraneBarostat, DrudeNoseHooverIntegrator, DrudeLangevinIntegrator, MonteCarloBarostat

from rflow.integrators import NoseHooverChainVelocityVerletIntegrator
from rflow import RickFlow, FreeEnergyCosineSeries, RelativePartialCenterOfMassRestraint
from fflip.omm.playpara import *
from fflip.omm.paragroup import *
from fflip.omm.util import *
from fflip.drude import *


options = get_md_options('input.mdo')
mdo = parse_md_options(options)

psf_file = mdo['psf'] 
crd_file = mdo['crd']
surface_tension = mdo['surf_ts'] 
xvec = mdo['xvec']
yvec = mdo['yvec']
zvec = mdo['zvec']
lipname = mdo['lipid']
lipid = parse_lipid(lipname)
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

top_files = glob.glob("toppar/*.str")

workflow = RickFlow(
    toppar=top_files,
    psf=psf_file,
    crd=crd_file,
    box_dimensions=[xvec, yvec, zvec], 
    # box_dimensions=None, 
    gpu_id=0,
    switch_distance=0*u.angstrom,
    cutoff_distance=10*u.angstrom,
    vdw_switching="openmm",
    nonbonded_method=LJPME,
    center_around='not water',
    dcd_output_interval=1000,
    table_output_interval=1000,
    steps_per_sequence=1000000,
    initialize_velocities=False,
    misc_psf_create_system_kwargs={"ewaldErrorTolerance": 0.0001}
)

temperature = temp * u.kelvin

#
#  ========= Define Integrator and Barostat ==========
#

drude_temperature = 1.0 * u.kelvin
timestep = 0.001 * u.picosecond
if intgrt == 'N':
    friction = 10.0 / u.picosecond # That's how I interpret tau = 0.1 ps from Lamoureux and Roux
    drude_friction = 1 #/ u. picosecond  
    integrator = DrudeNoseHooverIntegrator(
           temperature, friction, drude_temperature, drude_friction, timestep,
           3,3,3
    )
elif intgrt == 'L':
    integrator = DrudeLangevinIntegrator(
        temperature, 5/u.picosecond, drude_temperature,
        20/u.picosecond, 0.001 * u.picoseconds
    )
integrator.setMaxDrudeDistance(0.02)  # hard wall constraint; recommendation from Jing Huang

if baro == 'MCM':
    barostat = MonteCarloMembraneBarostat(
        1.0 * u.atmosphere, surface_tension, 
        integrator.getTemperature(), 
        0, zmode, 25
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
    sol = filter_solution('solution.txt')
    parameter_sets = lipid.parse_groups()
    all_offsets = [
        gen_sensitivity_offset(
            ps, percentage=sol[i]
        ) for i, ps in enumerate(parameter_sets)
    ]
    topology = md.Topology.from_openmm(workflow.psf.topology)
    for pset, offset in zip(parameter_sets, all_offsets):
        change_drude_ff_parameters(
            workflow._system, topology, pset, offset
        )
    change_nb_exceptions(workflow.psf, workflow.system, isdrude=True)

if change_para.lower() == 'yes':
    from torfix import *
    # workflow._system = BrutalTorsionParameters(workflow.system, workflow.psf, tfixes)

#
#  ========= Run Simulation  ==========
#

newarea = False
if workflow.next_seqno == 1:
    newarea = True
    integrator.setStepSize(0.0001 * u.picoseconds)
    workflow.prepareSimulation(integrator, barostat)
    # workflow.simulation.minimizeEnergy()
else:
    workflow.prepareSimulation(integrator, barostat)

workflow.context.applyConstraints(1e-7)
workflow.context.computeVirtualSites()

if __name__ == "__main__":
    workflow.run()
    if newarea:
        os.system("area.py --new")
    else:
        os.system("area.py")

