# -*- coding: utf-8 -*-

import numpy as np
import mdtraj as md
from sklearn.cluster import DBSCAN
from simtk.openmm.app import CharmmPsfFile
from fflip.analysis.util import manually_select_res_atom
from matplotlib import pyplot as plt
import matplotlib.patches as patches


def color_(head):
    if head == 'pe':
        return 'yellow'
    elif head == 'pc':
        return 'yellow'
    elif head == 'pg':
        return 'yellow'
    elif head == 'pi':
        return 'yellow'
    elif head == 'ps':
        return 'yellow'
    elif head == 'sito':
        return 'b'
    elif head == 'cer':
        return 'r'


def rpatom(head):
    phos = ['P', 'O11', 'O12', 'O13', 'O14']
    if head.lower() == 'cer':
        return ['C2S']
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
                 radius, min_lipids=3, skip_every_n_frames=10,
                 do_top=True, do_bot=True, clfix=dict()):
        self.min_lipids = min_lipids
        self.skip = skip_every_n_frames
        self.do_top = do_top
        self.do_bot = do_bot
        self.clfix = clfix
        self.psf_file = psf_file
        self.radius = radius
        psf = CharmmPsfFile(psf_file)
        self.topology = md.Topology.from_openmm(psf.topology)
        self.top_reslib = []
        self.bot_reslib = []
        top_atoms = dict()
        bot_atoms = dict()
        top_resnames = dict()
        bot_resnames = dict()
        res_count_top = 0
        res_count_bot = 0
        for res in top_res_info:
            self.top_reslib.append(res[0])
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
            self.bot_reslib.append(res[0])
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
                d_list.append(self.radius[res1_head] + self.radius[res2_head])
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
                d_list.append(self.radius[res1_head] + self.radius[res2_head])
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

    def write_xy(self, traj, itraj, frame, xy_file='{}_{}.xy'):
        # abc, xyz here are for a single frame
        frame = frame - 1
        topstring = ""
        botstring = ""
        xyz = traj.xyz[frame]
        abc = self.box_edges(traj)[frame]
        for i, comb in enumerate(self.comb_vectors):
            for index in self.top_atoms:
                x = xyz[self.top_atoms[index], 0].mean() + comb[0] * abc[0]
                y = xyz[self.top_atoms[index], 1].mean() + comb[1] * abc[1]
                resname = self.top_res[index]
                topstring += "{0:>5s}{1:>4d}{2:>3d}{3:>10.5f}{4:>10.5f}\n".\
                    format(resname[0], index, i + 1, x, y)
            for index in self.bot_atoms:
                x = xyz[self.bot_atoms[index], 0].mean() + comb[0] * abc[0]
                y = xyz[self.bot_atoms[index], 1].mean() + comb[1] * abc[1]
                resname = self.bot_res[index]
                botstring += "{0:>5s}{1:>4d}{2:>3d}{3:>10.5f}{4:>10.5f}\n".\
                    format(resname[0], index, i + 1, x, y)
        with open('top_' + xy_file.format(itraj, frame + 1), 'w') as f:
            f.write(str(abc[0]) + '\t' + str(abc[1]) + '\n')
            f.write(topstring)
        with open('bot_' + xy_file.format(itraj, frame + 1), 'w') as f:
            f.write(str(abc[0]) + '\t' + str(abc[1]) + '\n')
            f.write(botstring)

    def plot_use_xy(
            self, itraj, frame, leaflet, xy_file='{}_{}_{}.xy',
            plot=True, edge=1.5, addlabel=True, highcontrast=False
    ):
        with open(xy_file.format(leaflet, itraj, frame), 'r') as fr:
            abc_ = fr.readline()
        xcell = float(abc_.split()[0])
        ycell = float(abc_.split()[1])
        data_ = np.loadtxt(
            xy_file.format(leaflet, itraj, frame), dtype=str, skiprows=1
        )
        xy_ = data_[:, 3:].astype(np.float)
        assert data_.shape[0] % 9 == 0
        x_ = xy_[:, 0]
        xmax = x_.max()
        xmin = x_.min()
        x_ = x_ - (xmax + xmin) / 2
        y_ = xy_[:, 1]
        ymax = y_.max()
        ymin = y_.min()
        y_ = y_ - (ymax + ymin) / 2
        r2 = []
        cutoff = []
        for a_ in range(len(x_)):
            r2.append([])
            cutoff.append([])
            index1 = a_ % int(data_.shape[0] / 9) + 1
            res1_head = self.top_res[index1][1].split('.')[1].lower()
            for b_ in range(len(x_)):
                r2[a_].append(
                    (x_[a_] - x_[b_])**2 + (y_[a_] - y_[b_])**2
                )
                index2 = b_ % int(data_.shape[0] / 9) + 1
                res2_head = self.top_res[index2][1].split('.')[1].lower()
                if res1_head + '-' + res2_head in self.clfix:
                    cutoff[a_].append(self.clfix[res1_head + '-' + res2_head])
                elif res2_head + '-' + res1_head in self.clfix:
                    cutoff[a_].append(self.clfix[res2_head + '-' + res1_head])
                else:
                    cutoff[a_].append(
                        self.radius[res1_head] + self.radius[res2_head]
                    )
        dbstop = DBSCAN(
            eps=1.0, min_samples=self.min_lipids, metric='precomputed'
        )
        matrix = np.sqrt(np.array(r2)) / np.array(cutoff)
        dbstop.fit(matrix)
        if not plot:
            return dbstop  # add bot latter
        plt.figure(figsize=(2.25, 2.25))
        plt.xlim(-edge * xcell / 2, edge * xcell / 2)
        plt.ylim(-edge * ycell / 2, edge * ycell / 2)
        ax = plt.gca()
        for c_ in range(len(x_)):
            indexc = c_ % int(data_.shape[0] / 9) + 1
            res_head = self.top_res[indexc][1].split('.')[1].lower()
            if dbstop.labels_[c_] != -1:
                circle = plt.Circle(
                    (x_[c_], y_[c_]), self.radius[res_head],
                    color=color_(res_head)
                )
            else:
                if highcontrast:
                    circle = plt.Circle(
                        (x_[c_], y_[c_]), self.radius[res_head],
                        color='grey', fill=None
                    )
                else:
                    circle = plt.Circle(
                        (x_[c_], y_[c_]), self.radius[res_head],
                        color=color_(res_head), linewidth=0.8, fill=None
                    )
            ax.add_artist(circle)
        rect = patches.Rectangle(
            (-xcell/2, -ycell/2), xcell, ycell, linewidth=2, linestyle='--',
            edgecolor='k', facecolor='none', zorder=9
        )
        ax = plt.gca()
        ax.add_patch(rect)
        if addlabel:
            plt.xlabel('x [nm]', fontsize=11, labelpad=2)
            plt.ylabel('y [nm]', fontsize=11, labelpad=2)
        plt.savefig(
            'visualize_{}_{}_{}.png'.format(leaflet, itraj, frame),
            dpi=300, bbox_inches='tight'
        )

    def plot_cluster(
        self, itraj, frame, leaflets=['top', 'bot'], edge=1.5,
        addlabel=True, highcontrast=False
    ):
        for leaflet in leaflets:
            self.plot_use_xy(
                itraj, frame, leaflet, plot=True, edge=edge,
                addlabel=addlabel, highcontrast=highcontrast
            )

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
        """
        Returns: top/bot dicts in cluster per frame for resnames
        """
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

    def __call__(self, traj, save_labels=True, first_traj_index=1):
        self.count_traj += 1
        self.count_frame += int(traj.xyz.shape[0] / self.skip)
        if self.count_traj == 1:
            self.tcs, self.tncl, self.tlb, self.bcs, self.bncl, self.blb = \
                self.calc(traj)
            tlb = self.tlb
            blb = self.blb
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
        if save_labels and self.do_top:
            np.savetxt(
                './labels/top_labels_{}.txt'.format(
                    self.count_traj + first_traj_index - 1
                ),
                np.array(tlb),
                fmt='%3d'
            )
            self.top_label_saved = True
        if save_labels and self.do_bot:
            np.savetxt(
                './labels/bot_labels_{}.txt'.format(
                    self.count_traj + first_traj_index - 1
                ),
                np.array(blb),
                fmt='%3d'
            )
            self.bot_label_saved = True

    def plot_and_save(self, max_size=10):
        top_res_info, bot_res_info = self.analyze_res()
        plt.rc('font', size=9)
        if self.do_top:
            # First plot the distribution of the cluster sizes
            h, b = np.histogram(
                self.tcs, bins=max_size, range=(0, max_size), normed=True
            )
            np.savetxt(
                "topsize.txt", np.array([b[:-1], h]).swapaxes(0, 1), fmt='%.5f',
                header="Size    Perc"
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
            plt.title("Average (top): {0:.2f}".format(np.array(self.tncl).mean()))
            plt.savefig('topnum.png', dpi=330, bbox_inches='tight')
            plt.close()
            # Each residue type
            perc_comp_list_top = []
            perc_cls_list_top = []
            txtheader = "        "
            compline = "Comp    "
            clusline = "Incl    "
            diffline = "Gain    "
            for rcount, resname in enumerate(top_res_info):
                # save time series first
                np.savetxt(
                    '{}_vs_time_top.txt'.format(
                        resname.lower()), top_res_info[resname], fmt='%3d'
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
                diffline += "{0:<10.4f}".format(
                    top_res_info[resname].sum() / np.array(self.tncl).sum() -
                    len(self.top_res_indexes[resname]) / self.top_count
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
                f.write(
                    txtheader + '\n' + compline + '\n' + clusline + '\n'
                    + diffline + '\n'
                )
        if self.do_bot:
            # First plot the distribution of the cluster sizes
            h, b = np.histogram(
                self.bcs, bins=max_size, range=(0, max_size), normed=True
            )
            np.savetxt(
                "botsize.txt", np.array([b[:-1], h]).swapaxes(0, 1), fmt='%.5f',
                header="Size    Perc"
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
            plt.title("Average (bot): {0:.2f}".format(np.array(self.bncl).mean()))
            plt.savefig('botnum.png', dpi=330, bbox_inches='tight')
            plt.close()
            # Each residue type
            perc_comp_list_bot = []
            perc_cls_list_bot = []
            txtheader = "        "
            compline = "Comp    "
            clusline = "Incl    "
            diffline = "Gain    "
            for rcount, resname in enumerate(bot_res_info):
                # save time series first
                np.savetxt(
                    '{}_vs_time_bot.txt'.format(
                        resname.lower()), bot_res_info[resname], fmt='%3d'
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
                    bot_res_info[resname].sum() / np.array(self.bncl).sum()
                )
                diffline += "{0:<10.4f}".format(
                    bot_res_info[resname].sum() / np.array(self.bncl).sum() -
                    len(self.bot_res_indexes[resname]) / self.bot_count
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
                f.write(
                    txtheader + '\n' + compline + '\n' + clusline + '\n'
                    + diffline + '\n'
                )

    def analyze_more(self, first_traj_index, last_traj_index):
        # Top leaflet
        if self.do_top:
            assert self.top_label_saved
            topfile = open('top_detail.txt', 'w')
            header = "     Trj     Frm   Label    Size"
            for resname_ in self.top_reslib:
                header += "{0:>8s}".format(resname_)
            topfile.write(header + '\n')
            for i in range(first_traj_index, last_traj_index + 1):
                top_res = np.array(
                    [self.top_res[i][0] for i in range(1, self.top_count + 1)],
                    dtype='str'
                )
                labels_traj = np.loadtxt(
                    'labels/top_labels_{}.txt'.format(i)
                )
                for fi, labels in enumerate(labels_traj):
                    num_clusters = labels.max()
                    for lb in range(int(num_clusters)):
                        size = (labels == lb).sum()
                        res_count = ""
                        for resname_ in self.top_reslib:
                            res_count += "{0:>8d}".format(
                                ((labels == lb) * (top_res == resname_)).sum()
                            )
                        topfile.write(
                            "{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4}\n".format(
                                i, self.skip * fi + 1, lb + 1, size, res_count
                            )
                        )
            topfile.close()
        # Bottom leaflet
        if self.do_bot:
            assert self.bot_label_saved
            botfile = open('bot_detail.txt', 'w')
            header = "     Trj     Frm   Label    Size"
            for resname_ in self.bot_reslib:
                header += "{0:>8s}".format(resname_)
            botfile.write(header + '\n')
            for i in range(first_traj_index, last_traj_index + 1):
                bot_res = np.array(
                    [self.bot_res[i][0] for i in range(1, self.bot_count + 1)],
                    dtype='str'
                )
                labels_traj = np.loadtxt(
                    'labels/bot_labels_{}.txt'.format(i)
                )
                for fi, labels in enumerate(labels_traj):
                    num_clusters = labels.max()
                    for lb in range(int(num_clusters)):
                        size = (labels == lb).sum()
                        res_count = ""
                        for resname_ in self.bot_reslib:
                            res_count += "{0:>8d}".format(
                                ((labels == lb) * (bot_res == resname_)).sum()
                            )
                        botfile.write(
                            "{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4}\n".format(
                                i, self.skip * fi + 1, lb + 1, size, res_count
                            )
                        )
            botfile.close()
