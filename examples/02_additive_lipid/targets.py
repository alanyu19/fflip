#!/usr/bin/env python

from fflip.ljpme.targetprop import *
from fflip.ljpme.optimizers import *
from def_lipids import *
import os

### Use absolute paths for these:
# the path for outputs (except simulation)
root = os.path.realpath('./core_files')
# the path of psf/crd files
psf_dir = os.path.realpath('./psf_files')
crd_dir = os.path.realpath('./crd_files')
# the path of sim/observable/potential calculation templates
template_dir = os.path.realpath('./templates')


# --------------- Example Targets ----------------
properties = []
special_properties = []

# DPPC areas
for temperature, surface_tension, wt in zip(
        [323.15, 333.15], [0, 0], [20, 10]
):
    prop = TargetProperty(
        prop_type="area",
        # property naming should follow this example!
        # keep the surface_tension as 0 if you run a NPT simulation,
        # this is needed by the code
        name="area_dppc_bilayer_{}_{}".format(surface_tension, temperature),
        system_type='bilayer',
        lipname="dppc",
        lipid=pc_demo,  # imported from def_lipids
        num_lipids=36,
        weight_factor=wt,
        temperature=temperature,
        surface_tension=surface_tension,
        psf_file=os.path.join(psf_dir, "bi_72_dppc_c36.psf"),
        crd_file=os.path.join(crd_dir, "bi_72_dppc_c36.crd"),
        root_dir=root,
        pot_template=os.path.join(template_dir, "potential"),
        obs_template=os.path.join(template_dir, "area"),
        sim_template=os.path.join(template_dir, "sim"),
        obs_file_format="area_{}.dat"
    )
    properties.append(prop)


# DPPC thicknesses
for temperature in [323.15, 333.15]:
    prop = TargetProperty(
        prop_type="db",
        name="db_dppc_bilayer_0_{}".format(temperature),
        system_type='bilayer',
        lipname="dppc",
        lipid=pc_demo,  # imported from def_lipids
        num_lipids=36,
        weight_factor=5,
        temperature=temperature,
        surface_tension=0,
        psf_file=os.path.join(psf_dir, "bi_72_dppc_c36.psf"),
        crd_file=os.path.join(crd_dir, "bi_72_dppc_c36.crd"),
        root_dir=root,
        pot_template=os.path.join(template_dir, "potential"),
        obs_template=os.path.join(template_dir, "db1"),
        obs_file_format="db_{}.dat"
    )
    properties.append(prop)


if __name__ == "__main__":
    for i, p in enumerate(properties):
        print("{}: {}".format(i, p.name))
