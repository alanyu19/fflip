# -*- coding: utf-8 -*-

import numpy as np
import mdtraj as md
from sklearn.cluster import DBSCAN
from simtk.openmm.app import CharmmPsfFile
from fflip.analysis.util import manually_select_res_atom
from matplotlib import pyplot as plt


# Just like LJ radius
def radius(head):
    if head == 'pe':
        return 0.29
    elif head == 'pc':
        return 0.29
    elif head == 'pg':
        return 0.29
    elif head == 'pi':
        return 0.29
    elif head == 'ps':
        return 0.29
    elif head == 'sito':
        return 0.29


def rpatom(head):
    phos = ['P', 'O11', 'O12', 'O13', 'O14']
    if head.lower() == 'cer':
        return ['O1', 'HO1']
    elif head.lower() == 'sito':
        return ["O3", "H3\'"]
    elif head.lower() == 'pe':
        return phos + ['C11', 'C12', 'N']
    elif head.lower() == 'pc':
        return phos + ['C11', 'C12', 'N']
    elif head.lower() == 'pg':
        return phos + ['C11', 'C12', 'C13', 'OC2', 'OC3']
    elif head.lower() == 'pi':
        return phos + ['C11', 'C12', 'C13', 'C14', 'O4', 'C15', 'C16']
    elif head.lower() == 'ps':
        return phos + ['C11', 'C12', 'C13', 'N']
    else:
        raise Exception('Not anticipated head type: {}'.format(head))


def select_identifier(psf_, topology_, atom_selection):
    if "\'" in atom_selection:
        atoms = manually_select_res_atom(
            psf_, atom_selection
        )
    else:
        atoms = topology_.select(atom_selection)
    return atoms


def select_representative(psf_, topology_, res, atom_names):
    part1 = "resname {} and ".format(res)
    part2 = ""
    for i, a in enumerate(atom_names):
        if i == 0:
            part2 = "name {}".format(a)
        else:
            part2 = part2 + " or name {}".format(a)
    sele_words = part1 + "({})".format(part2)
    if "\'" in sele_words:
        atoms = manually_select_res_atom(
            psf_, sele_words
        )
    else:
        atoms = topology_.select(sele_words)
    return atoms


def filter_top_bot_atoms(candidates, recipe):
    if 'explicit_grouping' in recipe:
        groups = recipe['explicit_grouping']
        assert isinstance(groups, tuple)
        assert len(groups) == 2
        assert isinstance(groups[0], list)
        assert isinstance(groups[1], list)
        top_atoms = []
        bot_atoms = []
        for atom in candidates:
            if atom in groups[0]:
                top_atoms.append(atom)
            else:
                assert atom in groups[1]
                bot_atoms.append(atom)
        return top_atoms, bot_atoms
    elif 'hard_border' in recipe:
        assert isinstance(recipe['hard_border'], int)
        top_atoms = []
        bot_atoms = []
        for atom in candidates:
            if atom <= recipe['hard_border']:
                top_atoms.append(atom)
            else:
                bot_atoms.append(atom)
        return top_atoms, bot_atoms
    else:
        raise Exception


def idatom(liptype):
    scheme = {
        'g': 'C2',
        'c': 'O3',
        's': 'C2S',
    }
    assert liptype in scheme
    return scheme[liptype]


class ClusterLip(object):
    def __init__(self, psf_file, top_res_info, bot_res_info, leaflet_recipe,
                 min_lipids=3, skip_every_n_frames=10,
                 do_top=True, do_bot=True, clfix=dict()):
        self.min_lipids = min_lipids
        self.skip = skip_every_n_frames
        self.do_top = do_top
        self.do_bot = do_bot
        self.clfix = clfix
        self.psf_file = psf_file
        psf = CharmmPsfFile(psf_file)
        self.topology = md.Topology.from_openmm(psf.topology)
        top_atoms = dict()
        bot_atoms = dict()
        top_resnames = dict()
        bot_resnames = dict()
        res_count_top = 0
        res_count_bot = 0
        for res in top_res_info:
            atoms = select_identifier(
                self.psf_file, self.topology,
                atom_selection="resname {} and name {}".format(
                    res[0], idatom(res[1].split('.')[0])
                )
            )
            # ids/id_atoms are used to find the specific residue used in the
            # following steps
            top_id_atoms, dummy = filter_top_bot_atoms(atoms, leaflet_recipe)
            # top_ids[res[0]] = top_id_atoms
            rp_atom_names = rpatom(res[1].split('.')[1])
            rp_atoms = select_representative(
                self.psf_file, self.topology, res[0], rp_atom_names
            )
            # then group rp_atoms to their mother residues
            for index, tida in enumerate(top_id_atoms):
                res_count_top += 1
                # store this info to retrieve later
                top_resnames[res_count_top] = res
                top_atoms[res_count_top] = []
                ref_residue = self.topology.atom(tida).residue
                for rpa in rp_atoms:
                    if self.topology.atom(rpa).residue == ref_residue:
                        top_atoms[res_count_top].append(int(rpa))

        for res in bot_res_info:
            atoms = select_identifier(
                self.psf_file, self.topology,
                atom_selection='resname {} and name {}'.format(
                    res[0], idatom(res[1].split('.')[0])
                )
            )
            dummy, bot_id_atoms = filter_top_bot_atoms(atoms, leaflet_recipe)
            # bot_ids[res[0]] = bot_id_atoms
            rp_atom_names = rpatom(res[1].split('.')[1])
            rp_atoms = select_representative(
                self.psf_file, self.topology, res[0], rp_atom_names
            )
            # then group rp_atoms to their mother residues
            for index, tida in enumerate(bot_id_atoms):
                res_count_bot += 1
                bot_resnames[res_count_bot] = res
                bot_atoms[res_count_bot] = []
                ref_residue = self.topology.atom(tida).residue
                for rpa in rp_atoms:
                    if self.topology.atom(rpa).residue == ref_residue:
                        bot_atoms[res_count_bot].append(int(rpa))
        self.top_atoms = top_atoms  # a dictionary
        self.bot_atoms = bot_atoms  # a dictionary
        self.top_res = top_resnames
        self.bot_res = bot_resnames
        self.top_count = res_count_top
        self.bot_count = res_count_bot
        self.count_traj = 0
        self.count_frame = 0

    def _top_scales(self):
        d_list = []
        for p in self.top_pairs:
            res1_head = self.top_res[p[0]][1].split('.')[1].lower()
            res2_head = self.top_res[p[1]][1].split('.')[1].lower()
            if res1_head + '-' + res2_head in self.clfix:
                d_list.append(self.clfix[res1_head + '-' + res2_head])
            elif res2_head + '-' + res1_head in self.clfix:
                d_list.append(self.clfix[res2_head + '-' + res1_head])
            else:
                d_list.append(radius(res1_head) + radius(res2_head))
        return np.array(d_list)

    @property
    def top_scales(self):
        return self._top_scales()

    def _bot_scales(self):
        d_list = []
        for p in self.bot_pairs:
            res1_head = self.bot_res[p[0]][1].split('.')[1].lower()
            res2_head = self.bot_res[p[1]][1].split('.')[1].lower()
            if res1_head + '-' + res2_head in self.clfix:
                d_list.append(self.clfix[res1_head + '-' + res2_head])
            elif res2_head + '-' + res1_head in self.clfix:
                d_list.append(self.clfix[res2_head + '-' + res1_head])
            else:
                d_list.append(radius(res1_head) + radius(res2_head))
        return np.array(d_list)

    @property
    def bot_scales(self):
        return self._bot_scales()

    @property
    def top_pairs(self):
        """the pair generator using sudo-resid generated in __init__"""
        pairs = []
        grid1, grid2 = np.meshgrid(
            np.arange(1, self.top_count + 1),
            np.arange(1, self.top_count + 1)
        )
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                pairs.append([grid1[i][j], grid2[i][j]])
        return np.array(pairs)

    @property
    def bot_pairs(self):
        """the pair generator using sudo-resid generated in __init__"""
        pairs = []
        grid1, grid2 = np.meshgrid(
            np.arange(1, self.bot_count + 1),
            np.arange(1, self.bot_count + 1)
        )
        for i in range(grid1.shape[0]):
            for j in range(grid1.shape[1]):
                pairs.append([grid1[i][j], grid2[i][j]])
        return np.array(pairs)

    @property
    def comb_vectors(self):
        """for PBC images"""
        comb_vectors = []
        for i in [0, 1, -1]:
            for j in [0, 1, -1]:
                comb_vectors.append([i, j])
        return comb_vectors

    @staticmethod
    def box_edges(traj):
        return traj.unitcell_lengths[:, :]
    
    def get_distances(self, abc, xyz, pairs, atoms):
        r2_all_images = []
        group1 = []
        group2 = []
        # first (0) elements in a pair
        for e1 in pairs.transpose()[0]:
            group1.append(atoms[e1])
        # second (1) elements in a pair
        for e2 in pairs.transpose()[1]:
            group2.append(atoms[e2])
        # loop over periodic images
        for comb in self.comb_vectors:
            r2_image = []
            for g1, g2 in zip(group1, group2):
                r2_temp = (
                    xyz[::self.skip, g1, 0].mean(axis=1) -
                    xyz[::self.skip, g2, 0].mean(axis=1) +
                    comb[0] * abc[::self.skip, 0]
                ) ** 2 + (
                    xyz[::self.skip, g1, 1].mean(axis=1) -
                    xyz[::self.skip, g2, 1].mean(axis=1) +
                    comb[1] * abc[::self.skip, 1]
                ) ** 2
                r2_image.append(r2_temp)
            r2_all_images.append(r2_image)
        # transpose to make frames row
        r2_min = np.min(np.array(r2_all_images), axis=0).transpose()
        return np.sqrt(r2_min)

    def calc(self, traj):
        abc = self.box_edges(traj)
        top_cluster_sizes = []
        top_num_clustering_lipids = []
        top_labels = []
        bot_cluster_sizes = []
        bot_num_clustering_lipids = []
        bot_labels = []
        # top leaflet
        if self.do_top:
            dbs_top = DBSCAN(
                eps=1.0, min_samples=self.min_lipids, metric='precomputed'
            )
            num_top_pairs = np.shape(self.top_pairs)[0]
            matrix_edge = int(np.sqrt(num_top_pairs))

            r2_all = self.get_distances(
                abc, traj.xyz, self.top_pairs, self.top_atoms
            )
            top_scales = self._top_scales()
            for r2_frame in r2_all:
                matrix = (r2_frame / top_scales).reshape(
                    (matrix_edge, matrix_edge)
                )
                dbs_top.fit(matrix)
                top_num_clustering_lipids.append(
                    (dbs_top.labels_ != -1).sum()
                )
                num_clusters = dbs_top.labels_.max()
                for lb in range(num_clusters):
                    top_cluster_sizes.append((dbs_top.labels_ == lb).sum())
                top_labels.append(dbs_top.labels_)
        # bottom leaflet
        if self.do_bot:
            dbs_bot = DBSCAN(
                eps=1.0, min_samples=self.min_lipids, metric='precomputed'
            )
            num_bot_pairs = np.shape(self.bot_pairs)[0]
            matrix_edge = int(np.sqrt(num_bot_pairs))

            r2_all = self.get_distances(
                abc, traj.xyz, self.bot_pairs, self.bot_atoms
            )
            bot_scales = self._bot_scales()
            for r2_frame in r2_all:
                matrix = (r2_frame / bot_scales).reshape(
                    (matrix_edge, matrix_edge)
                )
                dbs_bot.fit(matrix)
                bot_num_clustering_lipids.append(
                    (dbs_bot.labels_ != -1).sum()
                )
                num_clusters = dbs_bot.labels_.max()
                for lb in range(num_clusters):
                    bot_cluster_sizes.append((dbs_bot.labels_ == lb).sum())
                bot_labels.append(dbs_bot.labels_)
        return top_cluster_sizes, top_num_clustering_lipids, top_labels,\
            bot_cluster_sizes, bot_num_clustering_lipids, bot_labels

    @property
    def top_res_indexes(self):
        res_indexes = {}
        for r in self.top_res:
            if self.top_res[r][0] not in res_indexes:
                res_indexes[self.top_res[r][0]] = [r]
            else:
                res_indexes[self.top_res[r][0]].append(r)
        return res_indexes

    @property
    def bot_res_indexes(self):
        res_indexes = {}
        for r in self.bot_res:
            if self.bot_res[r][0] not in res_indexes:
                res_indexes[self.bot_res[r][0]] = [r]
            else:
                res_indexes[self.bot_res[r][0]].append(r)
        return res_indexes

    def analyze_res(self):
        # First assert some calculation is done
        assert self.count_traj >= 1
        # Start analysis
        top_res_info = {}
        bot_res_info = {}
        if self.do_top:
            top_res_info = {}
            labels_all_frames = np.array(self.tlb)
            in_cluster_all_frames = labels_all_frames != -1
            for resname in self.top_res_indexes:
                # convert to zero based
                indexes = np.array(self.top_res_indexes[resname]) - 1
                top_res_info[resname] = in_cluster_all_frames[:, indexes].sum(
                    axis=1
                )
        if self.do_bot:
            labels_all_frames = np.array(self.blb)
            in_cluster_all_frames = labels_all_frames != -1
            for resname in self.bot_res_indexes:
                # convert to zero based
                indexes = np.array(self.bot_res_indexes[resname]) - 1
                bot_res_info[resname] = in_cluster_all_frames[:, indexes].sum(
                    axis=1
                )
        return top_res_info, bot_res_info

    def __call__(self, traj):
        self.count_traj += 1
        self.count_frame += int(traj.xyz.shape[0] / self.skip)
        if self.count_traj == 1:
            self.tcs, self.tncl, self.tlb, self.bcs, self.bncl, self.blb = \
                self.calc(traj)
        else:
            tcs, tncl, tlb, bcs, bncl, blb = \
                self.calc(traj)
            # append
            self.tcs += tcs
            self.tncl += tncl
            self.tlb += tlb
            self.bcs += bcs
            self.bncl += bncl
            self.blb += blb

    def plot_and_save(self):
        top_res_info, bot_res_info = self.analyze_res()
        plt.rc('font', size=9)
        if self.do_top:
            # First plot the distribution of the cluster sizes
            h, b = np.histogram(
                self.tcs, bins=10, range=(0, 10), normed=True
            )
            plt.figure(figsize=(3, 2))
            plt.bar(b[:-1], h, color='b')
            plt.title('top size')
            plt.ylabel('Probability', rotation=90)
            plt.savefig('topsize.png', dpi=330, bbox_inches='tight')
            plt.close()
            # Then plot the average number of lipids involving clustering
            plt.figure(figsize=(3, 2))
            plt.plot(
                np.arange(0, len(self.tncl)) * self.skip, self.tncl, color='b'
            )
            plt.xlabel("Frame #")
            plt.ylabel("# of lipids in clusters", rotation=90)
            plt.title("Average (top): {}".format(np.array(self.tncl).mean()))
            plt.savefig('topnum.png', dpi=330, bbox_inches='tight')
            plt.close()
            # Each residue type
            perc_comp_list_top = []
            perc_cls_list_top = []
            txtheader = "        "
            compline = "Comp    "
            clusline = "Incl    "
            for rcount, resname in enumerate(top_res_info):
                # save time series first
                np.savetxt(
                    '{}_vs_time_top.txt'.format(
                        resname.lower()), top_res_info[resname]
                )
                txtheader += "{0:<10s}".format(resname)
                compline += "{0:<10.4f}".format(
                    len(self.top_res_indexes[resname]) / self.top_count
                )
                perc_comp_list_top.append(
                    len(self.top_res_indexes[resname]) / self.top_count
                )
                clusline += "{0:<10.4f}".format(
                    top_res_info[resname].sum() / np.array(self.tncl).sum()
                )
                perc_cls_list_top.append(
                    top_res_info[resname].sum() / np.array(self.tncl).sum()
                )
            # -- plot the percentage and compare with the composition
            plt.figure(figsize=(3, 2))
            plt.bar(
                np.arange(0, rcount + 1) - 0.16, perc_comp_list_top,
                label='Composition', width=0.3, color='gray'
            )
            plt.bar(
                np.arange(0, rcount + 1) + 0.16, perc_cls_list_top,
                label='In cluster', width=0.3, color='red'
            )
            plt.legend(fontsize=6)
            plt.xticks(
                range(0, len(self.top_res_indexes.keys())),
                self.top_res_indexes.keys(), fontsize=8
            )
            plt.title('Top %')
            plt.savefig('topperc.png', dpi=330, bbox_inches='tight')
            plt.close()
            # write out the same info to txt file
            with open('topperc.txt', 'w') as f:
                f.write(txtheader + '\n' + compline + '\n' + clusline)

        if self.do_bot:
            # First plot the distribution of the cluster sizes
            h, b = np.histogram(
                self.bcs, bins=10, range=(0, 10), normed=True
            )
            plt.figure(figsize=(3, 2))
            plt.bar(b[:-1], h, color='g')
            plt.title('bot size')
            plt.ylabel('Probability', rotation=90)
            plt.savefig('botsize.png', dpi=330, bbox_inches='tight')
            plt.close()
            # Then plot the average number of lipids involving clustering
            plt.figure(figsize=(3, 2))
            plt.plot(
                np.arange(0, len(self.bncl)) * self.skip, self.bncl, color='g'
            )
            plt.xlabel("Frame #")
            plt.ylabel("# of lipids in clusters", rotation=90)
            plt.title("Average (bot): {}".format(np.array(self.bncl).mean()))
            plt.savefig('botnum.png', dpi=330, bbox_inches='tight')
            plt.close()
            # Each residue type
            perc_comp_list_bot = []
            perc_cls_list_bot = []
            txtheader = "        "
            compline = "Comp    "
            clusline = "Incl    "
            for rcount, resname in enumerate(bot_res_info):
                # save time series first
                np.savetxt(
                    '{}_vs_time_bot.txt'.format(
                        resname.lower()), bot_res_info[resname]
                )
                txtheader += "{0:<10s}".format(resname)
                compline += "{0:<10.4f}".format(
                    len(self.bot_res_indexes[resname]) / self.bot_count
                )
                perc_comp_list_bot.append(
                    len(self.bot_res_indexes[resname]) / self.bot_count
                )
                clusline += "{0:<10.4f}".format(
                    bot_res_info[resname].sum() / np.array(self.bncl).sum()
                )
                perc_cls_list_bot.append(
                    top_res_info[resname].sum() / np.array(self.bncl).sum()
                )
            # -- plot the percentage and compare with the composition
            plt.figure(figsize=(3, 2))
            plt.bar(
                np.arange(0, rcount + 1) - 0.16, perc_comp_list_bot,
                label='Composition', width=0.3, color='gray'
            )
            plt.bar(
                np.arange(0, rcount + 1) + 0.16, perc_cls_list_bot,
                label='In cluster', width=0.3, color='darkorange'
            )
            plt.legend(fontsize=6)
            plt.xticks(
                range(0, len(self.bot_res_indexes.keys())),
                self.bot_res_indexes.keys(), fontsize=8
            )
            plt.title('Bottom %')
            plt.savefig('botperc.png', dpi=330, bbox_inches='tight')
            plt.close()
            # write out the same info to txt file
            with open('botperc.txt', 'w') as f:
                f.write(txtheader + '\n' + compline + '\n' + clusline)
