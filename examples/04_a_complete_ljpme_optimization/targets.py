#!/usr/bin/env python

from fflip.ljpme.targetprop import *
from fflip.ljpme.optimizers import *


root = './scratch/'
psf_dir = './psf_files/'
crd_dir = './crd_files/'
template_dir = './templates/'


# --------------- area_dppc_bilayer ----------------
properties = []
special_properties = []


# DPPC areas
for temperature, surface_tension, wt in zip(
        [323.15, 333.15], [0, 0], [20, 10]
):
    prop = TargetProperty(
        prop_type="area",
        name="area_dppc_bilayer_{}_{}".format(surface_tension, temperature),
        system_type='bilayer',
        lipname="dppc",
        num_lipids=36,
        weight_factor=wt,
        temperature=temperature,
        surface_tension=surface_tension,
        psf_file=psf_dir + "bi_72_dppc_c36.psf",
        crd_file=crd_dir + "bi_72_dppc_c36.crd",
        root_dir=root,
        pot_template=template_dir + "potential",
        obs_template=template_dir + "area",
        sim_template=template_dir + "sim",
        # TODO: pass this to the calc_observable function of TargetProperty
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
        num_lipids=36,
        weight_factor=5,
        temperature=temperature,
        surface_tension=0,
        psf_file=psf_dir + "bi_72_dppc_c36.psf",
        crd_file=crd_dir + "bi_72_dppc_c36.crd",
        root_dir=root,
        pot_template=template_dir + "potential",
        obs_template=template_dir + "db1",
        # TODO: pass this to the calc_observable function of TargetProperty
        obs_file_format="db_{}.dat"
    )
    properties.append(prop)


if __name__ == "__main__":
    for i, p in enumerate(properties):
        print("{}: {}".format(i, p.name))
