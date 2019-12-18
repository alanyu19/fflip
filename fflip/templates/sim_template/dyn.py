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
from coffe.omm.playpara import *
from coffe.omm.paragroup import *
from coffe.omm.util import *
from fflip.chm import *

md_options = get_md_options('input.mdo')

psf_file = str(md_options['psf'])

crd_file = str(md_options['crd'])

surface_tension = 20 * float(md_options['surface_tension'])

xvec = float(md_options['boxx'])

if 'boxy' in md_options:
    yvec = float(md_options['boxy'])
else:
    yvec = float(md_options['boxx'])

if 'boxz' in md_options:
    zvec = float(md_options['boxz'])
else:
    zvec = float(md_options['boxx'])

if 'lipid' in md_options:
    lipname = str(md_options['lipid'])
else:
    lipname = 'dppc' # maybe 'pc' is better?
    lipid = parse_lipid(lipname)

temp = float(md_options['temperature']) 

baro = str(md_options['barostat'])

if 'zmode' in md_options:
    zmode = int(md_options['zmode'])
else:
    zmode = 0
    if baro == 'MCM':
        warnings.warn('MonteCarloMembraneBarostat is used, but no zmode is specified, using default zmode(0)')

change_para = str(md_options['change_para'])

def simplify_solution(file_to_load = 'solution.txt', threshold = 0.1):
    data = np.loadtxt(file_to_load)
    sol = []
    for d in data:
        if np.abs(d) < threshold:
            sol.append(0)
        else:
            sol.append(d)
    return sol

#
#  ========= Setup System ==========
#

top_files = glob.glob("/u/alanyu/top_yalun/for_ljpme/original/*")

workflow = RickFlow(
    #toppar=["top_all36_lipid.rtf", "par_all36_lipid.prm", "toppar_water_ions.str"], # TODO
    toppar=top_files,
    psf=psf_file,
    crd=crd_file,
    box_dimensions=[xvec, yvec, zvec], # TODO
    gpu_id=0,
    #switch_distance=12*u.angstrom,
    cutoff_distance=10*u.angstrom,
    nonbonded_method=LJPME,
    dcd_output_interval=500,
    table_output_interval=500,
    steps_per_sequence=500000,
    use_vdw_force_switch=False,
    use_only_xml_restarts=True
)

temperature = temp * u.kelvin

#
#  ========= Define Integrator and Barostat ==========
#

integrator = NoseHooverChainVelocityVerletIntegrator(
        workflow.system, temperature, 50.0 / u.picosecond, 2.0 * u.femtosecond,
        chain_length=3, num_mts=3, num_yoshidasuzuki=3
)

integrator = LangevinIntegrator(
        temperature, 1/u.picosecond, 0.002*u.picosecond
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
if change_para.lower() == 'yes':
    sol = simplify_solution('solution.txt')
    parameter_sets=lipid.parse_gtcnp()
    all_offsets = [gen_sensitivity_offset(ps, percentage = sol[i]) for i, ps in enumerate(parameter_sets)]
    workflow.add_parameter_offsets(parameter_sets, all_offsets)

workflow.prepareSimulation(integrator, barostat)


#
#  ========= Run Simulation  ==========
#

if workflow.next_seqno == 1:
    workflow.simulation.minimizeEnergy()
if __name__ == "__main__":
    workflow.run()
