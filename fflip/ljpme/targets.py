# -*- coding: utf-8 -*-


from fflip.ljpme.targetprop import *
# from fflip.ljpme.optimizers import *


root = '/u/alanyu/c36ljpme/fflow/'
psf_dir = '/u/alanyu/c36ljpme/structure_bank/psf_files/'
crd_dir = '/u/alanyu/c36ljpme/structure_bank/crd_files/'
template_dir = '/u/alanyu/c36ljpme/new_templates/'

useful_props = get_avail_exp_prop_names()

dppc_scd_names = get_sim_scd_names(
    '/u/alanyu/c36ljpme/fflow/.old-calc/runner/scd_dppc/block_data/*'
)
rdf_names = get_rdf_names_as_properties(
    '/u/alanyu/c36ljpme/fflow/.old-calc/runner/rdf/block_data/sparse*'
)

correlated_scd = []
for c in range(5, 16):
    correlated_scd.append('scd-c2{}'.format(c))
    correlated_scd.append('scd-c3{}'.format(c))

# --------------- area_dppc_bilayer ----------------
properties = []
special_properties = []


# DPPC bilayers
for temperature, surface_tension, wt in zip(
        [323.15, 323.15, 323.15, 333.15], [0, 5, -5, 0], [10, 0, 0, 7.5]
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
        # TODO: pass this to the calc_observable function of TargetProperty
        obs_file_format="area_{}.dat"
    )
    properties.append(prop)

# DLPC bilayer at 303.15K
prop = TargetProperty(
    prop_type='area',
    name='area_dlpc_bilayer_0_303.15',
    system_type='bilayer',
    lipname='dlpc',
    num_lipids=36,
    weight_factor=5,
    temperature=303.15,
    surface_tension=0,
    psf_file=psf_dir + 'bi_72_dlpc_c36.psf',
    crd_file=crd_dir + 'bi_72_dlpc_c36.crd',
    root_dir=root,
    pot_template=template_dir + 'potential',
    obs_template=template_dir + 'area',
    obs_file_format='area_{}.dat'
)
properties.append(prop)

# DPPC monolayers
for temperature, surface_tension in zip([321.15, 321.15, 321.15], [18, 40, 55]):
    prop = TargetProperty(
        prop_type="area",
        name="area_dppc_monolayer_{}_{}".format(surface_tension, temperature),
        system_type='monolayer',
        lipname="dppc",
        num_lipids=36,
        weight_factor=3,
        temperature=temperature,
        surface_tension=surface_tension,
        psf_file=psf_dir + "mono_72_dppc_c36.psf",
        crd_file=crd_dir + "mono_72_dppc_c36.crd",
        root_dir=root,
        pot_template=template_dir + "potential",
        obs_template=template_dir + "area",
        obs_file_format="area_{}.dat"
    )
    properties.append(prop)

# SCDs
for name in dppc_scd_names:
    if name in useful_props:
        if name in correlated_scd and name != 'scd-c25':
            wf = 0.5
        elif name in ['c11']:
            wf = 0.5
        elif name in ['scd-c23', 'scd-c24', 'scd-c25', 'scd-c3-hx',
                      'scd-c22-h2r', 'scd-c22-h2s']:
            wf = 0.5
        else:
            wf = 0.5
        prop = TargetProperty(
            prop_type='scd',
            name=name,
            system_type='bilayer',
            lipname='dppc',
            num_lipids=36,
            weight_factor=wf,
            temperature=323.15,
            surface_tension=0,
            psf_file=psf_dir + 'bi_72_dlpc_c36.psf',
            crd_file=crd_dir + 'bi_72_dlpc_c36.crd',
            root_dir=root,
            pot_template=template_dir + "potential",
            obs_template=template_dir + "dppc_scd",  # TODO
            obs_file_format=name + '_{}.dat'
        )
        properties.append(prop)

# RDFs
for name in rdf_names:
    short_name = name.split('_')[0]
    prop = TargetProperty(
        prop_type='rdf',
        name=name,
        system_type='bulk',
        lipname='prpc',
        num_lipids=9,
        weight_factor=1,
        temperature=298.15,
        surface_tension=0,
        psf_file=psf_dir + 'bulk_9_prpc_c36.psf',
        crd_file=crd_dir + 'bulk_9_prpc_c36.crd',
        root_dir=root,
        pot_template=1,  # TODO
        obs_template=1,  # TODO
        obs_file_format='sparse-' + short_name + '-{}.txt'
    )
    properties.append(prop)

special_properties = []

prop = SpecialProperty('dppc_bilayer_ka',
                       323.15,
                       1.5,
                       root,
                       parent_properties=properties[0:3],
                       exp_rel_dir='exp'
                       )
special_properties.append(prop)