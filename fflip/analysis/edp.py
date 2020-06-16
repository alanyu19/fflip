# -*- coding: utf-8 -*-


from rflow.edp import *
from fflip.analysis.edp_util import *
from time import sleep
from random import randint


class ElectronDensityFactory:
    """
    Do production calculation and grouping
    """
    def __init__(self, psf_file, traj_template,
                 com_selection=None, box_size=10,
                 save_to_folder='.', sep=None):
        """
        Args:
            psf_file: str, the psf file
            traj_template: str, the trajectory template
            com_selection: str, the center of mass selection words (ref to
            mdtraj for this)
            save_to_folder: str, the folder to save all the atom data,
            usually don't need to change it.
        """
        self.psf_file = psf_file
        self.traj_template = traj_template
        self.com_selection = com_selection
        self.top = make_topology(psf_file)
        self.residues = find_resnames_from_psf(psf_file)
        if not os.path.isdir(save_to_folder):
            os.mkdir(save_to_folder)
        self.save_to_folder = save_to_folder
        self.box_size = box_size
        # this is currently a toy code, only accept
        self.sep = sep

    def __call__(self, first, last, skip=[]):
        """
        Args:
            first: int, the first trajectory index
            last: int, the last trajectory index
        Returns:
            None
        """
        sleep(randint(0, 5))
        trajs = TrajectoryIterator(
            first_sequence=first, last_sequence=last,
            filename_template=self.traj_template,
            topology_file=self.psf_file,
            atom_selection="all", load_function=md.load_dcd
        )
        for res in self.residues:
            if self.sep is None:
                subdir = os.path.join(
                    self.save_to_folder,
                    './atoms_{}'.format(res.lower())
                )
                if not os.path.isdir(subdir):
                    os.mkdir(subdir)
                blocks_dir = os.path.join(subdir, 'blocks')
                if not os.path.isdir(blocks_dir):
                    os.mkdir(blocks_dir)
            else:
                if res.lower() == 'tip3':
                    continue  # don't calculate water since (no interdig)
                subdir_upper = os.path.join(
                    self.save_to_folder,
                    './atoms_{}_upper'.format(res.lower())
                )
                subdir_lower = os.path.join(
                    self.save_to_folder,
                    './atoms_{}_lower'.format(res.lower())
                )
                if not os.path.isdir(subdir_upper):
                    os.mkdir(subdir_upper)
                if not os.path.isdir(subdir_lower):
                    os.mkdir(subdir_lower)
                blocks_dir_upper = os.path.join(subdir_upper, 'blocks')
                blocks_dir_lower = os.path.join(subdir_lower, 'blocks')
                if not os.path.isdir(blocks_dir_upper):
                    os.mkdir(blocks_dir_upper)
                if not os.path.isdir(blocks_dir_lower):
                    os.mkdir(blocks_dir_lower)

            atoms = find_atoms_from_psf(self.psf_file, res)
            for atom in atoms:
                if self.sep is None:
                    edc = ElectronDensityCalculator(
                        atom_selection="resname {} and name {}".format(
                            res.upper(), atom.upper()
                        ),
                        nbins=int(500 * (self.box_size / 10)),
                        com_selection=self.com_selection,
                        box_length_fixed=self.box_size,
                        topology_file=self.psf_file
                    )
                    for traj in trajs:
                        edc(traj)
                    np.savetxt(
                        os.path.join(
                            blocks_dir, "{}_{}.dat".format(atom, first)
                        ),
                        edc.adp_angstrom
                    )
                    np.savetxt(
                        os.path.join(
                            blocks_dir, "{}_{}_edp.dat".format(atom, first)
                        ),
                        edc.edp_angstrom
                    )
                else:
                    edc_upper = ElectronDensityCalculator(
                        atom_selection=
                        "resname {} and name {} and resid 0 to {}".format(
                            res.upper(), atom.upper(), self.sep
                        ),
                        nbins=int(500 * (self.box_size / 10)),
                        com_selection=self.com_selection,
                        box_length_fixed=self.box_size,
                        topology_file=self.psf_file
                    )
                    edc_lower = ElectronDensityCalculator(
                        atom_selection=
                        "resname {} and name {} and resid > {}".format(
                            res.upper(), atom.upper(), self.sep
                        ),
                        nbins=int(500 * (self.box_size / 10)),
                        com_selection=self.com_selection,
                        box_length_fixed=self.box_size,
                        topology_file=self.psf_file
                    )
                    for traj in trajs:
                        edc_upper(traj)
                        edc_lower(traj)
                    np.savetxt(
                        os.path.join(
                            blocks_dir_upper, "{}_{}.dat".format(atom, first)),
                        edc_upper.adp_angstrom
                    )
                    np.savetxt(
                        os.path.join(
                            blocks_dir_upper, "{}_{}_edp.dat".format(
                                atom, first
                            )
                        ),
                        edc_upper.edp_angstrom
                    )
                    np.savetxt(
                        os.path.join(
                            blocks_dir_lower, "{}_{}.dat".format(atom, first)),
                        edc_lower.adp_angstrom
                    )
                    np.savetxt(
                        os.path.join(
                            blocks_dir_lower, "{}_{}_edp.dat".format(
                                atom, first
                            )
                        ),
                        edc_lower.edp_angstrom
                    )

            if self.sep is None:
                z_bin_file = os.path.join(subdir, 'z.txt')
                if not os.path.isfile(z_bin_file):
                    np.savetxt(

                        z_bin_file, np.linspace(
                        -edc.box_length_fixed / 2 + edc.box_length_fixed / edc.num_bins / 2,
                        edc.box_length_fixed / 2 - edc.box_length_fixed / edc.num_bins / 2,
                        edc.num_bins
                    )
                               )
            else:
                z_bin_file1 = os.path.join(subdir_upper, 'z.txt')
                z_bin_file2 = os.path.join(subdir_lower, 'z.txt')
                edc = edc_upper
                np.savetxt(z_bin_file1, np.linspace(
                    -edc.box_length_fixed / 2 + edc.box_length_fixed / edc.num_bins / 2,
                    edc.box_length_fixed / 2 - edc.box_length_fixed / edc.num_bins / 2,
                    edc.num_bins
                )
                           )
                edc = edc_lower
                np.savetxt(z_bin_file2, np.linspace(
                    -edc.box_length_fixed / 2 + edc.box_length_fixed / edc.num_bins / 2,
                    edc.box_length_fixed / 2 - edc.box_length_fixed / edc.num_bins / 2,
                    edc.num_bins
                )
                           )

