#!/usr/bin/env python

from fflip.ljpme.targetprop import *
from fflip.ljpme.optimizers import *
from def_lipids import *

# the path for outputs (except simulation)
root = os.path.realpath('./core_files/')
# the path of psf/crd files
psf_dir = os.path.realpath('./psf_files')
crd_dir = os.path.realpath('./crd_files')
# the path of sim/observable/potential calculation templates
template_dir = os.path.realpath('./templates')


# --------------- area_dppc_bilayer ----------------
properties = []
special_properties = []

# DPPC areas
for temperature, surface_tension, wt in zip(
        [323.15, 333.15], [0, 0], [10, 10]
):
    prop = TargetProperty(
        prop_type="area",
        name="area_dppc_bilayer_{}_{}".format(surface_tension, temperature),
        system_type='bilayer',
        lipid_name="dppc",
        lipid=match_lipid["dppc"],  # imported from def_lipids
        num_lipids=36,
        weight_factor=wt,
        temperature=temperature,
        surface_tension=surface_tension,
        psf_file=os.path.join(psf_dir, "dppc72.psf"),
        crd_file=os.path.join(crd_dir, "dppc72bi_centered.crd"),
        root_dir=root,
        pot_template=os.path.join(template_dir, "potential"),
        obs_template=os.path.join(template_dir, "area"),
        sim_template=os.path.join(template_dir, "sim"),
        # TODO: pass this to the calc_observable function of TargetProperty
        obs_file_format="area_{}.dat"
    )
    properties.append(prop)

if __name__ == "__main__":
    for i, p in enumerate(properties):
        print("{}: {}".format(i, p.name))
