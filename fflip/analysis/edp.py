# -*- coding: utf-8 -*-


import numpy as np
from rflow.observables import BinEdgeUpdater
from fflip.analysis.util import select_atoms, manually_select_atoms
from fflip.analysis.edp_util import *
from time import sleep
from random import randint


def smooth_edp(data, degree):
    triangle = np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1]))
    # up then down
    smoothed = []
    for i in range(degree, len(data) - degree * 2):
        point = data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    # Handle boundaries
    smoothed = [smoothed[0]] * int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed


class ElectronDensityCalculator(BinEdgeUpdater):
    """
    The electrom density calculator for only one atom name / type
    """
    def __init__(self, atom_selection, coordinate=2, nbins=500,
                 com_selection=None, use_fixed_box_length=True,
                 box_length_fixed=10, topology_file=None):
        """
        Make Sure that the atom_selection is one element only!!!
        Args:
            atom_selection:
            coordinate:
            nbins:
            com_selection: List of atom ids to calculate the com of the
            membrane, to make the distribution relative to
            the center of mass.
        """
        super().__init__(num_bins=nbins, coordinate=coordinate)
        self.atom_selection = atom_selection
        self.density = 0.0
        self.edensity = 0.0
        self.com_selection = com_selection
        self.other_coordinates = set((0, 1, 2))
        self.other_coordinates.remove(coordinate)
        self.other_coordinates = list(self.other_coordinates)
        self.average_cross_section = 0.0
        self.use_fixed_box_length = use_fixed_box_length
        self.box_length_fixed = box_length_fixed
        if "\'" in self.atom_selection:
            assert topology_file is not None and \
                isinstance(topology_file, str), \
                "Please provide a topology_file, it's possible that some " \
                "atoms can't be selected if you don't have that!"
            self.topology_file = topology_file

    def __call__(self, trajectory):
        # get the average box edge and total number of frames
        super().__call__(trajectory)
        # get the average cross section
        cross_section = \
            trajectory.unitcell_lengths[:, self.other_coordinates[0]] * \
            trajectory.unitcell_lengths[:, self.other_coordinates[1]]
        self.average_cross_section = self.previous_n_frames * \
                                     self.average_cross_section + \
                                     trajectory.n_frames * cross_section.mean()
        self.average_cross_section /= self.n_frames
        # handling bad charactor "'"
        # this a loose condition, but enough for now (yalun)
        if "\'" in self.atom_selection:
            atom_names = []
            select_words = self.atom_selection.split()
            for iw, w in enumerate(select_words):
                if 'name' in w and "\'" in select_words[iw+1]:
                    atom_names.append(select_words[iw+1])
            atom_ids = []
            for atom_name in atom_names:
                selection = manually_select_atoms(self.topology_file, atom_name)
                atom_ids += selection
            atom_ids = np.array(atom_ids)
        else:
            atom_ids = select_atoms(trajectory, self.atom_selection)
        # get the element, you can see the reason why we only allow one atom
        # type here
        atom_for_element = trajectory.topology.atom(atom_ids[0])
        com_ids = select_atoms(trajectory, self.com_selection)
        # normalized with respect to the box edge
        # NOTE: the normalized may exceed/not fill 1 when using box_length_fixed
        normalized = normalizecoor(
            trajectory, self.coordinate, subselect=atom_ids,
            com_selection=com_ids,
            use_fixed_box_length=self.use_fixed_box_length,
            box_length_fixed=self.box_length_fixed
        )
        # weights for histograming, inverse of the volume of each bin
        weights = 1 / cross_section.transpose()
        weights = np.tile(weights.reshape(-1, 1), (1, normalized.shape[1]))
        if self.use_fixed_box_length:
            weights = weights * self.num_bins / self.box_length_fixed
        weights = weights / trajectory.n_frames
        histogram = np.histogram(
            normalized, bins=self.num_bins, range=(0, 1),
            weights=weights
        )  # this is !much! faster than manual bins
        density = histogram[0]
        self.density = trajectory.n_frames * density + self.previous_n_frames \
            * self.density
        self.density /= self.n_frames
        edensity = histogram[0] * atom_for_element.element.number
        # atomic number is also the eletron number
        self.edensity = trajectory.n_frames * edensity + \
            self.previous_n_frames * self.edensity
        self.edensity /= self.n_frames
        if self.use_fixed_box_length:
            self.edp = self.edensity  # electron density
            self.adp = self.density  # just the atom density
        else:
            self.edp = self.edensity * self.num_bins / self.average_box_size
            self.adp = self.density * self.num_bins / self.average_box_size
        self.edp_angstrom = self.edp / 1000
        self.adp_angstrom = self.adp / 1000
        return np.array([histogram[0]])


class MembraneThickness(BinEdgeUpdater):
    def __init__(self, atom_selection='water',
                 coordinate=2, nbins=500, edge_bins=100,
                 com_selection='not water', use_fixed_box_length=True,
                 box_length_fixed=10, smooth_points=5):
        super().__init__(num_bins=nbins, coordinate=coordinate)
        self.atom_selection = atom_selection
        self.edge_bins = edge_bins
        self.density = 0.0
        self.edensity = 0.0
        self.com_selection = com_selection
        self.other_coordinates = set((0, 1, 2))
        self.other_coordinates.remove(coordinate)
        self.other_coordinates = list(self.other_coordinates)
        self.average_cross_section = 0.0
        self.use_fixed_box_length = use_fixed_box_length
        self.box_length_fixed = box_length_fixed
        self.smooth_points = smooth_points

    def __call__(self, trajectory):
        super().__call__(trajectory)
        # this a loose condition, but enough for now (yalun)
        atom_ids = select_atoms(trajectory, self.atom_selection)
        com_ids = select_atoms(trajectory, self.com_selection)
        # normalized with respect to the box edge
        normalized = normalizecoor(
            trajectory, self.coordinate, subselect=atom_ids,
            com_selection=com_ids,
            use_fixed_box_length=self.use_fixed_box_length,
            box_length_fixed=self.box_length_fixed
        )
        db_per_frame = []
        for fr in range(trajectory.n_frames):
            histogram = np.histogram(
                normalized[fr], bins=self.num_bins, range=(0, 1)
            )  # this is !much! faster than manual bins
            distrib = np.array([histogram[0]])[0]
            data = smooth_edp(distrib, 12)
            rho_max_w_n = np.amax(data[:int(self.num_bins/2)])
            rho_mid_w_n = rho_max_w_n / 2.0
            rho_max_w_p = np.amax(data[int(self.num_bins/2):])
            rho_mid_w_p = rho_max_w_p / 2.0
            db_n = np.interp(
                rho_mid_w_n, np.flipud(
                    data[self.edge_bins:int(self.num_bins/2)]
                ),
                np.flipud(
                    np.arange(
                        -5 * self.box_length_fixed * (
                                1 - (2 * self.edge_bins / self.num_bins)
                        ), 0, 0.2
                    )
                )
            )
            db_p = np.interp(
                rho_mid_w_p,
                data[int(self.num_bins/2):-self.edge_bins],
                np.arange(
                    0,
                    5 * self.box_length_fixed * (
                            1 - (2 * self.edge_bins / self.num_bins)
                    ), 0.2
                )
            )
            db_per_frame.append(db_p - db_n)
        return np.array(db_per_frame)


class MembraneDhh(BinEdgeUpdater):
    def __init__(self, atom_selection='resname ',
                 coordinate=2, nbins=500, edge_bins=100,
                 com_selection='not water', use_fixed_box_length=True,
                 box_length_fixed=10, smooth_points=5):
        super().__init__(num_bins=nbins, coordinate=coordinate)
        self.atom_selection = atom_selection
        self.edge_bins = edge_bins
        self.density = 0.0
        self.edensity = 0.0
        self.com_selection = com_selection
        self.other_coordinates = set((0, 1, 2))
        self.other_coordinates.remove(coordinate)
        self.other_coordinates = list(self.other_coordinates)
        self.average_cross_section = 0.0
        self.use_fixed_box_length = use_fixed_box_length
        self.box_length_fixed = box_length_fixed
        self.smooth_points = smooth_points
        self.masses = None

    def __call__(self, trajectory):
        super().__call__(trajectory)
        # this a loose condition, but enough for now (yalun)
        atom_ids = select_atoms(trajectory, self.atom_selection)
        com_ids = select_atoms(trajectory, self.com_selection)
        # normalized with respect to the box edge
        normalized = normalizecoor(
            trajectory, self.coordinate, subselect=atom_ids,
            com_selection=com_ids,
            use_fixed_box_length=self.use_fixed_box_length,
            box_length_fixed=self.box_length_fixed
        )
        if self.masses is None:
            for i, a in enumerate(trajectory.topology.atoms):
                assert i == a.index
            self.masses = np.array([atom.element.mass for atom in trajectory.topology.atoms])
        dhh_per_frame = []
        for fr in range(trajectory.n_frames):
            histogram = np.histogram(
                normalized[fr], bins=self.num_bins, range=(0, 1),
                weights=self.masses
            )  # this is !much! faster than manual bins
            distrib = np.array([histogram[0]])[0]
            data = smooth_edp(distrib, 12)
            rho_max_n = np.amax(data[:int(self.num_bins/2)])
            n_position = np.where(data==rho_max_n)
            rho_max_p = np.amax(data[int(self.num_bins/2):])
            p_position = np.where(data==rho_max_p)
            dhh_per_frame.append(
                self.box_length_fixed * (p_position[-1] - n_position[0])[0] * 10 / self.num_bins
            )
        return np.array(dhh_per_frame)


class ElectronDensityFactory:
    """
    Do production calculation and grouping
    """
    def __init__(self, psf_file, traj_template,
                 com_selection=None, box_size=10,
                 save_to_folder='.', sep=None, pop_drude=False):
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
        self.pop_drude = pop_drude

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
                if res.lower() == 'tip3' or res.lower() in skip:
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

            atoms_temp = find_atoms_from_psf(self.psf_file, res)
            if self.pop_drude:
                atoms = []
                for _atom in atoms_temp:
                    if _atom[0].upper() != 'D':
                        atoms.append(_atom)
            else:
                atoms = atoms_temp
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

