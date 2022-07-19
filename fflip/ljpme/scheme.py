# -*- coding: utf-8 -*-

import os
from fflip.chm import *
from fflip.drude import *


class FolderNamingScheme(object):
    def __init__(self, target_property):
        self.tp = target_property

    def exp_folder(self, rel_dir='exp'):
        return os.path.join(self.tp.root_dir, rel_dir)

    def robustness_dir(self):
        return os.path.join(self.tp.root_dir, "robustness")

    def reweight_dir(self):
        return os.path.join(self.tp.root_dir, "reweightings")

    def trajectory_folder(self, iteration, traj_root):
        if self.tp.traj_rep < 0:
            subdir = "iter{}/{}_{}_{}_{}".format(
                iteration, self.tp.lipid_name, self.tp.system_type,
                self.tp.surface_tension, self.tp.temperature
                )
        else:
            subdir = "iter{}/{}_{}_{}_{}_{}".format(
                iteration, self.tp.lipid_name, self.tp.system_type,
                self.tp.surface_tension, self.tp.temperature, self.tp.traj_rep
                )
        return os.path.join(traj_root, subdir)

    def property_data_folder(self, iteration):
        if self.tp.traj_rep < 0:
            subdir = "observables/iter{}/{}_{}_{}_{}_{}".format(
                iteration, self.tp.prop_type, self.tp.lipid_name,
                self.tp.system_type, self.tp.surface_tension,
                self.tp.temperature
                )
        else:
            subdir = "observables/iter{}/{}_{}_{}_{}_{}_{}".format(
                iteration, self.tp.prop_type, self.tp.lipid_name,
                self.tp.system_type, self.tp.surface_tension,
                self.tp.temperature, self.tp.traj_rep
                )
        return os.path.join(self.tp.root_dir, subdir)

    def potential_data_folder(self, iteration):
        if self.tp.traj_rep < 0:
            subdir = "potentials/iter{}/potential_{}_{}_{}_{}_{}".format(
                iteration, self.tp.lipid_name, self.tp.system_type,
                self.tp.surface_tension, self.tp.temperature,
                self.tp.perturbation
                )
        else:
            subdir = "potentials/iter{}/potential_{}_{}_{}_{}_{}_{}".format(
                iteration, self.tp.lipid_name, self.tp.system_type,
                self.tp.surface_tension, self.tp.temperature,
                self.tp.perturbation, self.tp.traj_rep
                )
        return os.path.join(self.tp.root_dir, subdir)

    def reweighting_folder(self, iteration):
        return os.path.join(
            self.tp.reweight_dir, 'iter{}/{}'.format(
                iteration, self.tp.perturbation
            )
        )

    def reweighting_file_name(self):
        if self.tp.traj_rep < 0:
            return "{}_{}_{}_{}_{}".format(
                self.tp.prop_type, self.tp.lipid_name, self.tp.system_type,
                self.tp.surface_tension, self.tp.temperature
                )
        else:
            return "{}_{}_{}_{}_{}_{}".format(
                self.tp.prop_type, self.tp.lipid_name, self.tp.system_type,
                self.tp.surface_tension, self.tp.temperature, self.tp.traj_rep
                )
    def robustness_folder(self, iteration):
        return os.path.join(
            self.tp.robdir, 'iter{}/{}'.format(iteration, self.tp.perturbation)
        )

    def robustness_diff_file(self, iteration):
        if self.tp.traj_rep < 0:
            return os.path.join(
                self.tp.robdir, 'iter{}/{}/{}_diff.txt'.format(
                    iteration, self.tp.perturbation, self.tp.name
                )
            )
        else:
            return os.path.join(
                self.tp.robdir, 'iter{}/{}/{}_{}_diff.txt'.format(
                    iteration, self.tp.perturbation, self.tp.name,
                    self.tp.traj_rep
                )
            )

    def robustness_std_file(self, iteration):
        if self.tp.traj_rep < 0:
            return os.path.join(
                self.tp.robdir, 'iter{}/{}/{}_std.txt'.format(
                    iteration, self.tp.perturbation, self.tp.name
                )
            )
        else:
            return os.path.join(
                self.tp.robdir, 'iter{}/{}/{}_{}_std.txt'.format(
                    iteration, self.tp.perturbation, self.tp.name,
                    self.tp.traj_rep, self.tp.traj_rep
                )
            )


class SimOptScheme(object):
    def __init__(self, target_system):
        self.ts = target_system

    @property
    def critical_seqno(self):
        """Return:
            The first seqno to do the sensitivity analysis,
            the last seqno of the simulation
        """
        finder = {
            "drude_dppc_bilayer_0_323.15": (51, 200),
            "drude_dppc_bilayer_5_323.15": (51, 200),
            "drude_dppc_bilayer_-5_323.15": (51, 200),
            "drude_dppc_monolayer_18_321.15": (51, 200),
            "drude_dppc_monolayer_40_321.15": (51, 200),
            "drude_dppc_monolayer_55_321.15": (51, 200),
            "drude_dlpc_bilayer_0_303.15": (51, 200),
            "drude_dmpc_bilayer_0_303.15": (51, 200),
            "drude_popc_bilayer_0_303.15": (51, 200),
            "drude_dopc_bilayer_0_303.15": (51, 200),
            "drude_prpc_bulk_0_298.15": (11, 100),
            "dppc_bilayer_0_323.15": (51, 200),
            "dppc_bilayer_-5_323.15": (61, 300),
            "dppc_bilayer_5_323.15": (61, 300),
            "dppc_bilayer_0_333.15": (51, 200),
            "dppc_monolayer_18_321.15": (51, 200),
            "dppc_monolayer_40_321.15": (51, 200),
            "dppc_monolayer_55_321.15": (51, 200),
            "dlpc_bilayer_0_303.15": (51, 200),
            "dmpc_bilayer_0_303.15": (51, 200),
            "popc_bilayer_0_303.15": (51, 200),
            "prpc_bulk_0_298.15": (11, 100)
        }
        seqno_ = finder["{}_{}_{}_{}".format(
            self.ts.lipid_name, self.ts.system_type,
            self.ts.surface_tension, self.ts.temperature
        )]
        return seqno_

    @property
    def boxx(self):
        finder = {
            "drude_dppc_bilayer_0": 48,
            "drude_dppc_bilayer_-5": 48,
            "drude_dppc_bilayer_5": 48,
            "drude_dppc_monolayer_18": 45,
            "drude_dppc_monolayer_40": 49,
            "drude_dppc_monolayer_55": 52,
            "drude_dlpc_bilayer_0": 48,
            "drude_dmpc_bilayer_0": 48,
            "drude_popc_bilayer_0": 49,
            "drude_dopc_bilayer_0": 50,
            "drude_prpc_bulk_0": 42,
            "dppc_bilayer_0": 45.9,
            "dppc_bilayer_-5": 46,
            "dppc_bilayer_5": 46,
            "dppc_monolayer_18": 45,
            "dppc_monolayer_40": 49,
            "dppc_monolayer_55": 52,
            "dlpc_bilayer_0": 44.95,
            "dmpc_bilayer_0": 49.14,  # 45.6,
            "dspc_bilayer_0": 46.75,
            "popc_bilayer_0": 49.63,  # 46.4,
            "sopc_bilayer_0": 46,
            "prpc_bulk_0": 43
        }
        return finder["{}_{}_{}".format(
            self.ts.lipid_name, self.ts.system_type, self.ts.surface_tension
        )]

    @property
    def boxz(self):
        finder = {
            "drude_dppc_bilayer_0": 73,
            "drude_dppc_bilayer_-5": 73,
            "drude_dppc_bilayer_5": 73,
            "drude_dppc_monolayer_18": 223,
            "drude_dppc_monolayer_40": 223,
            "drude_dppc_monolayer_55": 223,
            "drude_dlpc_bilayer_0": 62,
            "drude_dmpc_bilayer_0": 63,
            "drude_popc_bilayer_0": 65,
            "drude_dopc_bilayer_0": 65,
            "drude_prpc_bulk_0": None,
            "dppc_bilayer_0": 72.7,
            "dppc_bilayer_-5": 64,
            "dppc_bilayer_5": 64,
            "dppc_monolayer_18": 223,
            "dppc_monolayer_40": 223,
            "dppc_monolayer_55": 223,
            "dlpc_bilayer_0": 67.25,
            "dmpc_bilayer_0": 56.20,  # 68.5,
            "dspc_bilayer_0": 78.85,
            "popc_bilayer_0": 64.32,  # 72.5,
            "sopc_bilayer_0": 78.8,
            "prpc_bulk_0": None
        }
        return finder["{}_{}_{}".format(
            self.ts.lipid_name, self.ts.system_type, self.ts.surface_tension
        )]

    @property
    def integrator(self):
        return 'L'

    @property
    def barostat(self):
        if self.ts.system_type != 'bulk':
            return "MCM"
        else:
            return "MC"

    @property
    def zmode(self):
        if self.ts.system_type == "monolayer":
            return 1
        elif self.ts.system_type == "bilayer":
            return 0
        else:
            assert self.ts.system_type == "bulk"
            return None


class LipidScheme(object):
    def __init__(self, lipid):
        self.lipid = lipid

    def parse(self, **kwargs):
        if type(self.lipid) == DrudeLipid:
            return self.lipid.parse_groups(**kwargs)
        elif type(self.lipid) == Lipid:
            return self.lipid.parse_groups(**kwargs)

    def param_ids(self, **kwargs):
        if type(self.lipid) == DrudeLipid:
            groups = self.lipid.parse_groups(**kwargs)
            ids = []
            for g in groups:
                ids.append(str(g.cgid) + '_' + str(g.internal_id))
            return ids
        elif type(self.lipid) == Lipid:
            return list(
                range(1, len(self.lipid.parse_groups(**kwargs)) + 1)
            )


# -------------------------------- Defaults ---------------------------------

def make_guess_of_block_size(btype, prop_type):
    if btype == 0:
        if prop_type == 'rdf':
            return 1
        else:
            return 10
    elif btype == 1:
        return 10
    else:
        raise Exception("type {} not supported, use 0 or 1".format(btype))


def make_guess_of_trajectory_range(name):
    # TODO: could add a check for last trajectory (through energy_rel_dir)
    if 'scd' in name:
        return 0, 0
    elif 'area' in name and 'bilayer' in name:
        return 0, 0
    elif 'area' in name and 'mono' in name:
        return 0, 0
    else:
        return 0, 0


def make_guess_of_intervals(name):
    if 'scd' in name:
        return 10, 10
    elif 'area' in name:
        return 10, 10
    elif 'db' in name or 'dhh' in name:
        return 10, 10
    else:
        return 10, 1


def make_guess_of_scaling(name):
    if 'area' in name:
        # (nm^2)
        return 1 / 0.6
    elif 'peak' in name:
        return 1 / 1
    elif 'foot' in name:
        return 1 / 0.5
    elif 'rmsd' in name:
        return 1 / 1
    elif 'scd' in name:
        return 1 / 0.15
    elif 'db' in name or 'dhh' in name:
        return 1 / 40
    elif 'ka' in name:
        return 1 / 200
    elif 'delta_area' in name:
        return 1 / 0.02


def make_guess_of_layertype(name):
    if 'bilayer' in name:
        return 'bilayer'
    elif 'mono' in name:
        return 'monolayer'
    elif 'scd' in name:
        return 'bilayer'
    else:
        return 'bulk'
