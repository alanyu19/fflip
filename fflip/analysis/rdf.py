# -*- coding: utf-8 -*-

import os
import numpy as np
import mdtraj as md
from simtk.openmm.app import CharmmPsfFile


def manually_select_res_atom(topology_file, atom_selection):
    if "\'" in atom_selection:
        # Do our own selection to avoid error in mdtraj for H3' and O3' in
        # sterols and possible future bad naming of atoms
        select_words = atom_selection.split()
        atom_names = []
        res_names = []
        for iw, w in enumerate(select_words):
            if 'resname' in w:
                res_names.append(select_words[iw + 1])
            if 'name' in w and "\'" in select_words[iw + 1]:
                atom_names.append(select_words[iw + 1])
        atom_ids = []
        with open(topology_file, 'r') as topfile:
            lines = topfile.readlines()
        for line in lines:
            for res in res_names:
                if res.upper() in line:
                    for atom in atom_names:
                        if atom.upper() in line:
                            atom_ids.append(int(line.strip().split()[0]))
        return np.array(atom_ids)


class RDF(object):
    def __init__(self, psf, name, atom_selections,
                 r_range=(0.0, 2.0), bin_width=0.005,
                 dimension=3, separate_leaflet=False,
                 method='yalun', save_blocks=True):
        """
        Args:
            psf: the psf file (str) or CharmmPsfFile
            name: str, the name users give to the RDF
            atom_selections:
            r_range: list of two floats, in nanometer
            bin_width: float, the bin width
            dimension: int, the dimension, can be 2 or 3
            separate_leaflet: bool, if want to calculate two leaflets separately
            method: str, 'mdtraj' or 'yalun' (default)
            save_blocks: bool, if saving the block data
        """
        if isinstance(psf, CharmmPsfFile):
            psf = psf
        else:
            assert isinstance(psf, str)
            self.psf_file = psf
            psf = CharmmPsfFile(psf)
        self.topology = md.Topology.from_openmm(psf.topology)
        self.r_range = list(r_range)
        self.bin_width = bin_width
        self.bins = np.arange(self.r_range[0], self.r_range[1], self.bin_width)
        self.name = name
        self.sele_words = atom_selections
        if "\'" in atom_selections[0]:
           self.atom1 = manually_select_res_atom(
               self.psf_file, atom_selections[0]
           )
        else:
            self.atom1 = self.topology.select(atom_selections[0])
        if "\'" in atom_selections[1]:
            self.atom2 = manually_select_res_atom(
                self.psf_file, atom_selections[1]
            )
        else:
            self.atom2 = self.topology.select(atom_selections[1])
        self.natom1 = self.atom1.shape[0]
        self.natom2 = self.atom2.shape[0]
        self.dimension = dimension
        # the leaflet feature is for lipid bilayer RDF
        self.separate_leaflet = separate_leaflet
        self.method = method
        self.save_blocks = save_blocks
        self.count_traj = 0

    @property
    def v_bin_size(self):
        if self.dimension == 3:
            return 4 * np.pi * self.bin_width 
        elif self.dimension == 2:
            return 2 * np.pi * self.bin_width
        else:
            raise Exception(
                "Dimension ({})not supported!".format(self.dimension)
            )
    
    @property
    def pairs(self):
        """the atom pairs generator"""
        pairs = []
        grid1, grid2 = np.meshgrid(self.atom1, self.atom2)
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                if grid1[i][j] != grid2[i][j]:
                    pairs.append([grid1[i][j], grid2[i][j]])
        return np.array(pairs)
    
    @property
    def comb_vectors(self):
        """for PBC images"""
        comb_vectors = []
        if self.dimension == 3:
            for i in [0, 1, -1]:
                for j in [0, 1, -1]:
                    for k in [0, 1, -1]:
                        comb_vectors.append([i, j, k])
        elif self.dimension == 2:
            for i in [0, 1, -1]:
                for j in [0, 1, -1]:
                    comb_vectors.append([i, j])
        return comb_vectors

    @staticmethod
    def make_pairs(atoms1, atoms2):
        pairs = []
        grid1, grid2 = np.meshgrid(atoms1, atoms2)
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                if not grid1[i][j] == grid2[i][j]:
                    pairs.append([grid1[i][j], grid2[i][j]])
        return pairs

    def up_low_atoms_from_recentered(self, traj):
        # attention:
        # 1) can not be used for water and those have a high diffusion constant
        # 2) assumption: the bilayer is lying in the xy plane
        # 3) must use centered trajectory!
        n_frames = traj.n_frames
        zmean_1 = np.mean(traj.xyz[int(n_frames / 2), self.atom1, 2])
        zmean_2 = np.mean(traj.xyz[int(n_frames / 2), self.atom2, 2])
        upper_atom1 = self.atom1[
            np.where(traj.xyz[int(n_frames / 2), self.atom1, 2] > zmean_1)
        ]
        upper_atom2 = self.atom2[
            np.where(traj.xyz[int(n_frames / 2), self.atom2, 2] > zmean_2)
        ]
        lower_atom1 = self.atom1[
            np.where(traj.xyz[int(n_frames / 2), self.atom1, 2] <= zmean_1)
        ]
        lower_atom2 = self.atom2[
            np.where(traj.xyz[int(n_frames / 2), self.atom2, 2] <= zmean_2)
        ]
        return upper_atom1, upper_atom2, lower_atom1, lower_atom2

    def up_low_pairs_from_recentered(self, traj):
        # attention:
        # 1) can not be used for water and those have a high diffusion constant
        # 2) assumption: the bilayer is lying in the xy plane
        # 3) must use centered trajectory!
        upper_atom1, upper_atom2, lower_atom1, lower_atom2 = \
            self.up_low_atoms_from_recentered(traj)
        upper_pairs = self.make_pairs(upper_atom1, upper_atom2)
        lower_pairs = self.make_pairs(lower_atom1, lower_atom2)
        return np.array(upper_pairs), np.array(lower_pairs)

    def natom_up_low_from_recentered(self, traj):
        """
        For getting the bulk density, assume that the diffusion is not very
        fast, so that the bulk density won't change (statistically),
        this (and the upper_lower_pairs method) can be improved by per frame
        calculation, but will take significantly longer time.
        """
        upper_atom1, upper_atom2, lower_atom1, lower_atom2 = \
            self.up_low_atoms_from_recentered(traj)
        natom1_upper = np.shape(upper_atom1)[0]
        natom1_lower = np.shape(lower_atom1)[0]
        natom2_upper = np.shape(upper_atom2)[0]
        natom2_lower = np.shape(lower_atom2)[0]
        return natom1_upper, natom2_upper, natom1_lower, natom2_lower

    def hard_up_low_atoms(self, recipe):
        if 'explicit_grouping' in recipe:
            groups = recipe['explicit_grouping']
            assert isinstance(groups, tuple)
            assert len(groups) == 2
            assert isinstance(groups[0], list)
            assert isinstance(groups[1], list)
            upper_atom1 = []
            lower_atom1 = []
            for atom in self.atom1:
                if atom in groups[0]:
                    upper_atom1.append(atom)
                else:
                    assert atom in groups[1]
                    lower_atom1.append(atom)
            upper_atom2 = []
            lower_atom2 = []
            for atom in self.atom2:
                if atom in groups[0]:
                    upper_atom2.append(atom)
                else:
                    assert atom in groups[1]
                    lower_atom2.append(atom)
            return upper_atom1, upper_atom2, lower_atom1, lower_atom2
        elif 'hard_border' in recipe:
            assert isinstance(recipe['hard_border'], int)
            upper_atom1 = []
            lower_atom1 = []
            for atom in self.atom1:
                if atom <= recipe['hard_border']:
                    upper_atom1.append(atom)
                else:
                    lower_atom1.append(atom)
            upper_atom2 = []
            lower_atom2 = []
            for atom in self.atom2:
                if atom <= recipe['hard_border']:
                    upper_atom2.append(atom)
                else:
                    lower_atom2.append(atom)
            return upper_atom1, upper_atom2, lower_atom1, lower_atom2
        else:
            pass  # being stupid for now

    def hard_up_low_pairs(self, recipe=dict()):
        upper_atom1, upper_atom2, lower_atom1, lower_atom2 = \
            self.hard_up_low_atoms(recipe)
        upper_pairs = self.make_pairs(upper_atom1, upper_atom2)
        lower_pairs = self.make_pairs(lower_atom1, lower_atom2)
        return np.array(upper_pairs), np.array(lower_pairs)

    def hard_natom_up_low(self, recipe=dict()):
        upper_atom1, upper_atom2, lower_atom1, lower_atom2 = \
            self.hard_up_low_atoms(recipe)
        natom1_upper = np.shape(upper_atom1)[0]
        natom1_lower = np.shape(lower_atom1)[0]
        natom2_upper = np.shape(upper_atom2)[0]
        natom2_lower = np.shape(lower_atom2)[0]
        return natom1_upper, natom2_upper, natom1_lower, natom2_lower

    @staticmethod
    def box_edges(traj):
        return traj.unitcell_lengths[:, :]
    
    def get_distances(self, abc, xyz, pairs):
        r2_all_images = []
        for comb in self.comb_vectors:
            if self.dimension == 3:
                num_pairs = pairs.shape[0]
                r2_image = (
                    xyz[:, pairs.transpose()[0], 0] -
                    xyz[:, pairs.transpose()[1], 0] +
                    comb[0] * np.tile(abc[:, 0], (num_pairs, 1)).transpose()
                ) ** 2 + (
                    xyz[:, pairs.transpose()[0], 1] -
                    xyz[:, pairs.transpose()[1], 1] +
                    comb[1] * np.tile(abc[:, 1], (num_pairs, 1)).transpose()
                ) ** 2 + (
                    xyz[:, pairs.transpose()[0], 2] -
                    xyz[:, pairs.transpose()[1], 2] +
                    comb[1] * np.tile(abc[:, 2], (num_pairs, 1)).transpose()
                ) ** 2

            elif self.dimension == 2:
                num_pairs = pairs.shape[0]
                r2_image = (
                    xyz[:, pairs.transpose()[0], 0] -
                    xyz[:, pairs.transpose()[1], 0] +
                    comb[0] * np.tile(abc[:, 0], (num_pairs, 1)).transpose()
                ) ** 2 + (
                    xyz[:, pairs.transpose()[0], 1] -
                    xyz[:, pairs.transpose()[1], 1] +
                    comb[1] * np.tile(abc[:, 1], (num_pairs, 1)).transpose()
                ) ** 2
            else:
                raise Exception(
                    "You are NOT studying relativity, select dimension under 4!"
                )
            r2_all_images.append(r2_image)
        r2_min = np.min(np.array(r2_all_images), axis=0).transpose()
        return np.sqrt(r2_min)
               
    def calc(self, traj, recentered, up_low_recipe, save_sparse):
        abc = self.box_edges(traj)
        # get bulk density, assume that the 2d density is on x-y plane
        volumns = (abc[:, 0] * abc[:, 1] * abc[:, 2]) if self.dimension == 3 \
            else (abc[:, 0] * abc[:, 1])
        if not self.separate_leaflet:
            dens_bulk = self.natom2 / volumns if not np.array_equal(
                self.atom1, self.atom2
            ) else (self.natom2 - 1) / volumns
            r2_min = self.get_distances(abc, traj.xyz, self.pairs)
            rdf = np.zeros(int(self.bins.shape[0] - 1))
            # create this even if there is no need to save the sparse RDF
            sparse = []
            for i in range(traj.n_frames):
                hist = np.histogram(r2_min[:, i], self.bins)
                # save sparse RDF each frame for reweighting
                if save_sparse:
                    sparse.append(hist[0])
                    # this is a very very entry level code,
                    # should include peak/valley index later
                rdf = rdf + hist[0] / (
                        self.v_bin_size * (
                            self.bins[:-1] + self.bin_width / 2
                        ) ** (int(self.dimension - 1))) / dens_bulk[i] / \
                    self.natom1
            rdf = rdf / traj.n_frames
            if save_sparse:
                np.savetxt(
                    'sparse-{}-{}.txt'.format(self.name, self.count_traj),
                    np.array(sparse)
                )
            return self.bins[:-1] + self.bin_width, rdf
        else:
            # upper/lower pairs
            if recentered:
                up_pairs, low_pairs = self.up_low_pairs_from_recentered(traj)
                natom1_upper, natom2_upper, natom1_lower, natom2_lower = \
                    self.natom_up_low_from_recentered(traj)
            else:
                up_pairs, low_pairs = self.hard_up_low_pairs(up_low_recipe)
                natom1_upper, natom2_upper, natom1_lower, natom2_lower = \
                    self.hard_natom_up_low(up_low_recipe)
            # bulk density
            dens_bulk_upper = natom2_upper / volumns if not \
                np.array_equal(self.atom1, self.atom2) else \
                (natom2_upper - 1) / volumns
            dens_bulk_lower = natom2_lower / volumns if not \
                np.array_equal(self.atom1, self.atom2) else \
                (natom2_lower - 1) / volumns
            r2_min_upper = self.get_distances(abc, traj.xyz, up_pairs)
            r2_min_lower = self.get_distances(abc, traj.xyz, low_pairs)
            rdf_upper = np.zeros(int(self.bins.shape[0] - 1))
            rdf_lower = np.zeros(int(self.bins.shape[0] - 1))
            for i in range(traj.n_frames):
                hist_upper = np.histogram(r2_min_upper[:, i], self.bins)
                hist_lower = np.histogram(r2_min_lower[:, i], self.bins)
                rdf_upper = rdf_upper + hist_upper[0] / (
                        self.v_bin_size * (
                            self.bins[:-1] + self.bin_width/ 2
                        ) ** (int(self.dimension - 1))
                ) / dens_bulk_upper[i] / natom1_upper
                rdf_lower = rdf_lower + hist_lower[0] / (
                        self.v_bin_size * (
                            self.bins[:-1] + self.bin_width / 2
                        )**(int(self.dimension - 1))
                ) / dens_bulk_lower[i] / natom1_lower
            rdf_upper = rdf_upper / traj.n_frames
            rdf_lower = rdf_lower / traj.n_frames
            rdf_avg = (rdf_upper * natom1_upper + rdf_lower * natom1_lower) / \
                      (natom1_upper + natom1_lower)
            return self.bins[:-1] + self.bin_width, \
                np.array([rdf_avg, rdf_upper, rdf_lower])

    def save_sparse(self, traj, begin, bin_width, subfolder):
        rdf_all_frames = []
        for frm in range(traj.n_frames):
            radius, frame_rdf = md.compute_rdf(
                traj[frm], self.pairs, self.r_range, bin_width=bin_width
            )
            rdf_all_frames.append(frame_rdf)
        np.savetxt('{}/sparse-{}-{}.txt'.format(
            subfolder, self.name, self.count_traj + begin - 1),
            np.array(rdf_all_frames)
        )
        np.savetxt('{}/r.txt'.format(subfolder), np.array(radius))

    def __call__(
        self, traj, first_trj_index=1,
        recentered=False, up_low_recipe=None,
        verbose=1, print_interval=10, save_blocks_interval=5, save_to=".",
        save_sparse=False, sparse_bin_width=0.02, sparse_subfolder="."
    ):
        if self.dimension == 2 and not recentered:
            assert up_low_recipe is not None,\
                "Please provide a recipe for telling upper leaflet from lower!"

        self.count_traj += 1
        if verbose >= 2:
            print('Calculating <{}> rdf for trajectory {} ...'.format(
                self.name, self.count_traj))
        elif verbose == 1:
            if (self.count_traj - 1) % print_interval == 0:
                print('Calculating <{}> rdf for trajectory {} ...'.format(
                    self.name, self.count_traj))
        else:
            pass

        if self.method == 'mdtraj':
            # only support 3d distances/rdf
            if save_sparse:
                self.save_sparse(
                    traj, first_trj_index, sparse_bin_width, sparse_subfolder
                )
            self.radius, rdf_tmp = md.compute_rdf(
                traj, self.pairs, self.r_range
            )
            # confirm that the calculation is correct here (future work ...)  
            if self.count_traj == 1:
                self.rdf = rdf_tmp
            else:
                self.rdf = (self.rdf * (self.count_traj - 1) + rdf_tmp) / \
                           self.count_traj

        elif self.method == 'yalun':
            self.radius, rdf_tmp = self.calc(
                traj, recentered, up_low_recipe, save_sparse
            )
            if self.count_traj == 1:
                self.rdf = rdf_tmp
            else:
                self.rdf = (self.rdf * (self.count_traj - 1) + rdf_tmp) / \
                           self.count_traj
            # Save block average RDF to text file
            if self.save_blocks:
                if (self.count_traj - 1) % save_blocks_interval == 0:
                    self.last_block_rdf = rdf_tmp
                else:
                    self.last_block_rdf = \
                        (self.last_block_rdf *
                         ((self.count_traj - 1) % save_blocks_interval) +
                         rdf_tmp) / \
                        ((self.count_traj - 1) % save_blocks_interval + 1)
                if self.count_traj % save_blocks_interval == 0:
                    if self.dimension == 2:
                        file_name = os.path.join(
                            save_to, 'rdf-{}-{}-{}.txt'.format(
                                self.name,
                                self.count_traj - save_blocks_interval + first_trj_index,
                                self.count_traj + first_trj_index - 1
                            )
                        )
                        np.savetxt(
                            file_name, np.swapaxes(
                                np.array(
                                    [self.radius,
                                     self.last_block_rdf[0],
                                     self.last_block_rdf[1],
                                     self.last_block_rdf[2]]
                                ), 0, 1
                            )
                        )
                    elif self.dimension == 3:
                        file_name = os.path.join(
                            save_to,'rdf-{}-{}-{}.txt'.format(
                                self.name,
                                self.count_traj - save_blocks_interval + first_trj_index,
                                self.count_traj + first_trj_index - 1
                            )
                        )
                        np.savetxt(
                            file_name, np.swapaxes(
                                np.array([self.radius, self.last_block_rdf]),
                                0, 1
                            )
                        )
        else:
            raise Exception('Method not accepted!')
    
    def plot(self, cut=1.2, color='blue'):
        from matplotlib import pyplot as plt
        plt.plot(self.radius[:int(cut / self.bin_width)],
                 self.rdf[:int(cut / self.bin_width)] if self.dimension == 3
                 else self.rdf[0, :int(cut / self.bin_width)], color=color)
        plt.title(self.name)
        plt.show()
    
    def save_to_file(self):
        file_name = 'rdf-{}.txt'.format(self.name)
        if self.dimension == 3:
            np.savetxt(file_name, np.swapaxes(
                np.array([self.radius, self.rdf]), 0, 1) 
                )
        elif self.dimension == 2:
            np.savetxt(file_name, np.swapaxes(
                np.array([self.radius, self.rdf[0], self.rdf[1], self.rdf[2]]),
                0, 1
            )
                      )


class Coordination(object):
    
    def __init__(self, psf, name, atom_selections, r_range=[0, 0.4],
                 dimension=3, separate_leaflet=False, method='count'):
        topology = md.Topology.from_openmm(psf.topology)
        self.topology = topology
        self.r_range = r_range
        # self.bin_width = bin_width
        # self.bins = np.arange(self.r_range[0], self.r_range[1],
        # self.bin_width)
        self.name = name
        self.sele_words = atom_selections
        self.atom1 = self.topology.select(atom_selections[0])
        self.atom2 = self.topology.select(atom_selections[1])
        self.natom1 = self.atom1.shape[0]
        self.natom2 = self.atom2.shape[0]
        self.dimension = dimension
        # the leaflet feature is for some bilayers
        self.separate_leaflet = separate_leaflet
        self.method = method
        self.count_traj = 0
    
    @property
    def pairs(self):
        pairs = []
        grid1, grid2 = np.meshgrid(self.atom1, self.atom2)
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                pairs.append([grid1[i][j], grid2[i][j]])
        return np.array(pairs)
    
    @property
    def comb_vectors(self):
        comb_vectors = []
        if self.dimension == 3:
            for i in [0, 1, -1]:
                for j in [0, 1, -1]:
                    for k in [0, 1, -1]:
                        comb_vectors.append([i,j,k])
        elif self.dimension == 2:
            for i in [0, 1, -1]:
                for j in [0, 1, -1]:
                    comb_vectors.append([i,j])
        return comb_vectors
    
    # attention: currently not very correct for water and those have a high diffusion constant
    # assumption, the bilayer is lying in the xy plane
    def upper_lower_pairs(self, traj):
        n_frames = traj.n_frames
        zmean_1 = np.mean(traj.xyz[int(n_frames/2), self.atom1, 2])
        zmean_2 = np.mean(traj.xyz[int(n_frames/2), self.atom2, 2])
        self.upper_atom1 = self.atom1[np.where(traj.xyz[int(n_frames/2), self.atom1, 2] > zmean_1)]
        self.upper_atom2 = self.atom2[np.where(traj.xyz[int(n_frames/2), self.atom2, 2] > zmean_2)]
        self.lower_atom1 = self.atom1[np.where(traj.xyz[int(n_frames/2), self.atom1, 2] <= zmean_1)]
        self.lower_atom2 = self.atom2[np.where(traj.xyz[int(n_frames/2), self.atom2, 2] <= zmean_2)]
        upper_pairs = []
        grid1, grid2 = np.meshgrid(self.upper_atom1, self.upper_atom2)
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                upper_pairs.append([grid1[i][j], grid2[i][j]])
        lower_pairs = []
        grid1, grid2 = np.meshgrid(self.lower_atom1, self.lower_atom2)
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                lower_pairs.append([grid1[i][j], grid2[i][j]])
        return np.array(upper_pairs), np.array(lower_pairs)

    def natom_upper_lower(self, traj):
        # for getting the bulk density, suppose the diffusion is not very fast,
        # so that the bulk density won't change (statistically), this (and upper_lower_pairs)
        # can be improved by per frame calculation, but will take significantly more time
        n_frames = traj.n_frames
        zmean_1 = np.mean(traj.xyz[int(n_frames/2), self.atom1, 2])
        zmean_2 = np.mean(traj.xyz[int(n_frames/2), self.atom2, 2])
        natom1_upper = np.shape(
            self.atom1[np.where(traj.xyz[int(n_frames/2), self.atom1, 2] > zmean_1)]
        )[0]
        natom1_lower = np.shape(
            self.atom1[np.where(traj.xyz[int(n_frames/2), self.atom1, 2] <= zmean_1)]
        )[0]
        natom2_upper = np.shape(
            self.atom2[np.where(traj.xyz[int(n_frames/2), self.atom2, 2] > zmean_2)]
        )[0]
        natom2_lower = np.shape(
            self.atom2[np.where(traj.xyz[int(n_frames/2), self.atom2, 2] <= zmean_2)]
        )[0]
        return natom1_upper, natom2_upper, natom1_lower, natom2_lower

    def box_edges(self, traj):
        return traj.unitcell_lengths[:,:]
    
    def get_distances(self, abc, xyz, pairs):
        # particle_a and particle_b should be the atom indexes
        r2_all_images = []
        for comb in self.comb_vectors:
            r2_image = []
            if self.dimension == 3:
                for a, b in pairs:                    
                    # xyz[frame, atom, axis]
                    if a != b:
                        r2_image.append((xyz[:, a, 0] - xyz[:, b, 0] + comb[0] * abc[:, 0])**2 + \
                                        (xyz[:, a, 1] - xyz[:, b, 1] + comb[1] * abc[:, 1])**2 + \
                                        (xyz[:, a, 2] - xyz[:, b, 2] + comb[2] * abc[:, 2])**2
                                       )
            elif self.dimension == 2:
                for a,  b in pairs:
                    if a != b:
                        r2_image.append((xyz[:, a, 0] - xyz[:, b, 0] + comb[0] * abc[:, 0])**2 + \
                                        (xyz[:, a, 1] - xyz[:, b, 1] + comb[1] * abc[:, 1])**2
                                       )
            else:
                print("Warning: you are on the planet earth, select dimension under 3")
            r2_all_images.append(r2_image)
        r2_min = np.min(np.array(r2_all_images), axis = 0)
        # print(r2_min.shape) r2_min[pairs, frames]
        return np.sqrt(r2_min)
               
    def calc(self, traj, calc_avg):
        abc = self.box_edges(traj) # T
        if self.separate_leaflet == False :
            r2_min = self.get_distances(abc, traj.xyz, self.pairs)
            count_all_pairs = ((r2_min > self.r_range[0]) & (r2_min <= self.r_range[1])).sum(axis=0)
            if calc_avg:
                count_sum = ((r2_min > self.r_range[0]) & (r2_min <= self.r_range[1])).sum()
            count_avg = count_sum/ self.natom1/ traj.n_frames
            count_in_range = count_all_pairs/ self.natom1
            return count_in_range, count_avg
        else:
            upper_pairs, lower_pairs = self.upper_lower_pairs(traj)
            natom1_upper, natom2_upper, natom1_lower, natom2_lower = self.natom_upper_lower(traj)
            r2_min_upper = self.get_distances(abc, traj.xyz, upper_pairs)
            r2_min_lower = self.get_distances(abc, traj.xyz, lower_pairs)
            count_all_upper = ((r2_min_upper > self.r_range[0]) & (r2_min_upper <= self.r_range[1])).sum(axis=0)
            count_all_lower = ((r2_min_lower > self.r_range[0]) & (r2_min_lower <= self.r_range[1])).sum(axis=0)
            count_sum_upper = ((r2_min_upper > self.r_range[0]) & (r2_min_upper <= self.r_range[1])).sum()
            count_sum_lower = ((r2_min_lower > self.r_range[0]) & (r2_min_lower <= self.r_range[1])).sum()
            count_in_range_upper = count_all_upper/ natom1_upper if natom1_upper != 0 else 0
            count_in_range_lower = count_all_lower/ natom1_lower if natom1_lower != 0 else 0
            count_in_range = (count_in_range_lower + count_in_range_upper)/ 2
            count_avg = (count_sum_lower + count_sum_upper)/ (natom1_lower + natom1_upper)/ traj.n_frames
            return count_in_range, count_avg
            
    def __call__(self, traj, verbose = 1, print_interval = 10, calc_avg = True):
        self.count_traj += 1
        if verbose >= 2:
            print('Calculating <{}> coordination number for trajectory {} ...'.format(
                self.name, self.count_traj))
        elif verbose == 1:
            if (self.count_traj - 1) % print_interval == 0:
                print('Calculating <{}> coordination number for trajectory {} ...'.format(
                    self.name, self.count_traj))
        else:
            pass
        if self.method == 'rdf':
            assert self.dimension == 3
            abc = self.box_edges(traj)
            # get bulk density, assume that the 2d density is on x-y plane
            volumns = (abc[:, 0] * abc[:, 1] * abc[:, 2])
            volumn_avg = np.mean(volumns)
            self.radius, rdf = md.compute_rdf(traj, self.pairs, self.r_range)
            dens_bulk = self.natom2/ volumn_avg if not np.array_equal(self.atom1, self.atom2) else \
            (self.natom2 - 1)/ volumn_avg
            coordination_number_avg = (4 * np.pi * dens_bulk * (self.radius[1] - self.radius[0]) * (
                self.radius[int(self.r_range[0]/0.005): int(self.r_range[1]/0.005)]
            )**2 * rdf[int(self.r_range[0]/0.005): int(self.r_range[1]/0.005)]).sum()
            if self.count_traj == 1:
                self.coordination_number_avg = coordination_number_avg
            else:
                self.coordination_number_avg = (self.coordination_number_avg * (self.count_traj - 1
                                                                               ) + coordination_number_avg
                                           )/ self.count_traj
            return None # Currently don't support per frame calculation
        elif self.method == 'count':
            coordination_number_per_frame, coordination_number_avg = self.calc(traj, calc_avg)
            if self.count_traj == 1:
                self.coordination_number_avg = coordination_number_avg
            else:
                self.coordination_number_avg = (self.coordination_number_avg * (self.count_traj - 1
                                                                               ) + coordination_number_avg
                                           )/ self.count_traj
            # this is for the Timeseries, so no need to do appending
            return coordination_number_per_frame
        else:
            print('Method not accepted!')
