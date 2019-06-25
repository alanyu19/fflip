# -*- coding: utf-8 -*-

import numpy as np
import mdtraj as md

class rdf(object):
    def __init__(self, psf, name, atom_selections, r_range = [0,2], bin_width = 0.005,
                 dimension = 3, separate_leaflet = False, method = 'yalun', save_blocks = True):
        topology = md.Topology.from_openmm(psf.topology)
        self.topology = topology
        self.r_range = r_range
        self.bin_width = bin_width
        self.bins = np.arange(self.r_range[0], self.r_range[1], self.bin_width)
        self.name = name
        self.sele_words = atom_selections
        self.atom1 = self.topology.select(atom_selections[0])
        self.atom2 = self.topology.select(atom_selections[1])
        self.natom1 = self.atom1.shape[0]
        self.natom2 = self.atom2.shape[0]
        self.dimension = dimension
        # the leatlet feature is for some bilayer rdfs
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
            print('f**k, dimension not supported!')
    
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
     
    def get_xyz(self, traj):
        # xyz[frame, atom, axis]
        return traj.xyz[:,self.atom1,:], traj.xyz[:,self.atom2,:]
    
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
        return np.sqrt(r2_min)
               
    def calc(self, traj):
        abc = self.box_edges(traj)
        # get bulk density, assume that the 2d density is on x-y plane
        volumns = (abc[:, 0] * abc[:, 1] * abc[:, 2]) if self.dimension == 3 else \
                  (abc[:, 0] * abc[:, 1])
        if self.separate_leaflet == False :
            dens_bulk = self.natom2/ volumns if not np.array_equal(self.atom1, self.atom2) else \
            (self.natom2 - 1)/ volumns
            r2_min = self.get_distances(abc, traj.xyz, self.pairs)
            rdf = np.zeros(int(self.bins.shape[0] - 1))
            for i in range(traj.n_frames):
                hist = np.histogram(r2_min[:, i], self.bins)
                rdf = rdf + hist[0]/ (self.v_bin_size * (self.bins[:-1] + self.bin_width/ 2
                                                        )**(int(self.dimension - 1))
                                    )/ dens_bulk[i]/ self.natom1
            rdf = rdf/ traj.n_frames
            return self.bins[:-1] + self.bin_width, rdf
        else:
            upper_pairs, lower_pairs = self.upper_lower_pairs(traj)
            natom1_upper, natom2_upper, natom1_lower, natom2_lower = self.natom_upper_lower(traj)
            dens_bulk_upper = natom2_upper/ volumns if not np.array_equal(self.atom1, self.atom2) else \
            (natom2_upper - 1)/ volumns
            dens_bulk_lower = natom2_lower/ volumns if not np.array_equal(self.atom1, self.atom2) else \
            (natom2_lower - 1)/ volumns
            r2_min_upper = self.get_distances(abc, traj.xyz, upper_pairs)
            r2_min_lower = self.get_distances(abc, traj.xyz, lower_pairs)
            rdf_upper = np.zeros(int(self.bins.shape[0] - 1))
            rdf_lower = np.zeros(int(self.bins.shape[0] - 1))
            for i in range(traj.n_frames):
                hist_upper = np.histogram(r2_min_upper[:, i], self.bins)
                hist_lower = np.histogram(r2_min_lower[:, i], self.bins)
                rdf_upper = rdf_upper + hist_upper[0]/ (self.v_bin_size * (self.bins[:-1] + self.bin_width/ 2
                                                                       )**(int(self.dimension - 1))
                                                    )/ dens_bulk_upper[i]/ natom1_upper
                rdf_lower = rdf_lower + hist_lower[0]/ (self.v_bin_size * (self.bins[:-1] + self.bin_width/ 2
                                                                       )**(int(self.dimension - 1))
                                                    )/ dens_bulk_lower[i]/ natom1_lower
            rdf_upper = rdf_upper/ traj.n_frames
            rdf_lower = rdf_lower/ traj.n_frames
            rdf_avg = (rdf_upper * natom1_upper + rdf_lower * natom1_lower
                      )/ (natom1_upper + natom1_lower)
            return self.bins[:-1] + self.bin_width, np.array([rdf_avg, rdf_upper, rdf_lower])
            
    def __call__(self, traj, begin, verbose = 1, print_interval = 10, save_blocks_interval = 5):
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
            self.radius, rdf_tmp = md.compute_rdf(traj, self.pairs, self.r_range)
            # confirm that the calculation is correct here (future work ...)  
            if self.count_traj == 1:
                self.rdf = rdf_tmp
            else:
                self.rdf = (self.rdf * (self.count_traj - 1) + rdf_tmp)/ self.count_traj
        elif self.method == 'yalun':
            self.radius, rdf_tmp = self.calc(traj)
            if self.count_traj == 1:
                self.rdf = rdf_tmp
            else:
                self.rdf = (self.rdf * (self.count_traj - 1) + rdf_tmp)/ self.count_traj
            if self.save_blocks == True:
                if (self.count_traj - 1) % save_blocks_interval == 0:
                    self.last_block_rdf = rdf_tmp
                else:
                    self.last_block_rdf = (self.last_block_rdf * ((self.count_traj - 1) % save_blocks_interval) + rdf_tmp
                                           )/ ((self.count_traj - 1) % save_blocks_interval + 1)
                if self.count_traj % save_blocks_interval == 0:
                    if self.dimension == 2:
                        np.savetxt('rdf-{}-{}-{}.txt'.format(
                            self.name, self.count_traj - save_blocks_interval + begin, self.count_traj + begin - 1
                        ), np.swapaxes(
                            np.array([self.radius, self.last_block_rdf[0], self.last_block_rdf[1], self.last_block_rdf[2]]),
                            0, 1
                        )
                        )
                    elif self.dimension == 3:
                        np.savetxt('rdf-{}-{}-{}.txt'.format(
                            self.name, self.count_traj - save_blocks_interval + begin, self.count_traj + begin - 1
                        ), np.swapaxes(np.array([self.radius, self.last_block_rdf]), 0, 1))
        else:
            print('Method not accepted!')
    
    def plot(self, cut = 1.2, color = 'blue'):
        from matplotlib import pyplot as plt
        plt.plot(self.radius[:int(cut/ self.bin_width)], 
                 self.rdf[:int(cut/ self.bin_width)] if self.dimension == 3 else \
                 self.rdf[0,:int(cut/ self.bin_width)], 
                 color = color)
        plt.title(self.name)
        plt.show()
    
    def save_to_file(self):
        file_name = 'rdf-{}.txt'.format(self.name)
        if self.dimension == 3:
            np.savetxt(file_name, np.array(self.radius, self.rdf))
        elif self.dimension == 2:
            np.savetxt(file_name, np.swapaxes(
                np.array([self.radius, self.rdf[0], self.rdf[1], self.rdf[2]]), 0, 1
            )
                      )