#!/usr/bin/env python

# Authors:Xiaohong Zhuang, Yalun Yu, Dr.Jeffery B. Klauda
# Date:9/27/2016
# Calculate the number of clusters and number of lipids in each cluster


import numpy as np
from sklearn.cluster import DBSCAN

Nframe = 1


# Find line number for (x_pt, y_pt)
def res_info(x_pt, y_pt, block_data):
    for entry in block_data:
        if "{:.5f}".format(float(x_pt)) == "{:.5f}".format(float(entry[3])) \
            and \
                "{:.5f}".format(float(y_pt)) == "{:.5f}".format(
            float(entry[4])
        ):
            resn = entry[0]
            resid = entry[1]
            nim = entry[2]
            return resn, resid, nim


def cluster_calculation(
        Dcut, index_trj, leaflet, nlip, nframe, minsp=3,
        box_file_template="{}/box-{}.txt", xy_file_temlate="{}/xy-{}.dat",
        num_file_template="{}/num-{}.txt", tot_file_template="{}/tot-{}.txt"
):
    print(
        'line\tIdyn\tIf\tIc\tIlip\tResn\tResid\tNimg\tX_coord\t    Y_coord'
    )

    # Box data
    filetoread1 = box_file_template.format(leaflet, index_trj)
    data1 = np.genfromtxt(filetoread1)
    left_all_frame = data1[:, 0]
    right_all_frame = data1[:, 1]
    # box_all_frame = data1[:, 2]

    # Read the data file
    filetoread2 = xy_file_temlate.format(leaflet, index_trj)
    Data = np.genfromtxt(filetoread2, dtype=str)
    info1s = []
    info2s = []
    for fr in range(nframe):
        info1, info2 = cluster_bridge(
            Data, left_all_frame[fr], right_all_frame[fr],
            Dcut, index_trj, fr, nlip, minsp=minsp
        )
        if info1 is not None:
            info1s.append(info1)
            info2s += info2
    file_to_write = num_file_template.format(leaflet, index_trj)
    with open(file_to_write, 'w') as f:
        # Write title, index of frame, Index of cluster,
        # number of lipids in each cluster
        f.write('%s\t%s\t%s\t%s\n' % ('Idyn', 'If2', 'Ic', 'Nlc'))
        for string in info2s:
            f.write(string)
    file_to_write = tot_file_template.format(leaflet, index_trj)
    with open(file_to_write, 'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\n' % ('Idyn', 'If1', 'Nlf', 'Ncf', 'Nlca'))
        for string in info1s:
            f.write(string)


def cluster_bridge(Data, left, right, Dcut, index_trj, frame, nlip, minsp=3):
    Idyn = index_trj
    lpf = nlip * 9  # lines per frame
    Fr = int(frame)  # starts from 0
    Data = Data[Fr * lpf:(Fr + 1) * lpf]
    data_ = Data[:, 3:]
    data = Data[:, 3:].astype(float)
    # Find number of rows and columns
    # NR = np.shape(data)[0]
    # NC = np.shape(data)[1]

    # calculte clustering of single block
    info1, info2 = run_dbscan(data, data_, Data, Idyn, Dcut, left, right, Fr, minsp=minsp)
    return info1, info2


def run_dbscan(data, data_, Data, Idyn, Dcut, left, right, frame, minsp=3):
    line = 0
    # Compute DBSCAN
    db = DBSCAN(eps=Dcut, min_samples=minsp).fit(data)
    # print db.get_params(deep=True)
    labels = db.labels_
    # components = db.components_
    # Number of clusters in labels, ignoring noise (-1) if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    clusters_1 = []
    clusters_2 = []

    # Find size of each cluster
    if n_clusters_ > 0:
        # Find clusters have (any) coordinate inside main image
        clusters_1 = []
        # loop over clusters
        for i in range(0, n_clusters_, 1):
            comp_1 = data[labels == i]
            comp_1_ = data_[labels == i]
            check1 = any(
                left < xy[0] < right and right > xy[1] > left for xy in comp_1
            )
            if check1:
                clusters_1.append(comp_1)

    # Find clusters not dulipicate due to the boudary ,i.e.
    # not share any residue, use the last cluster
    return_info_1 = None
    if len(clusters_1) > 0:
        clusters_2 = []
        Lresid_ac = []
        # Check if resid dulplicate inside each cluster
        for i in range(0, len(clusters_1), 1):
            comp_1c = clusters_1[i]
            # Find resid of all lipids in the clusters
            Lresid_1c = []
            for p in range(0, len(comp_1c), 1):
                resn, resid, nim = res_info(comp_1c[p][0], comp_1c[p][1], Data)
                Lresid_1c.append(resid)
            for q in range(0, len(comp_1c), 1):
                Occ = Lresid_1c.count(Lresid_1c[q])
                if Occ > 1:
                    for j in range(0, len(comp_1c), 1):
                        # Find line in the coordinate file
                        resn, resid, nim = res_info(
                            comp_1c[j][0], comp_1c[j][1], Data
                        )
                        # print line, frame, cluster,
                        # component and x y coordinate in the cluster
                        print(
                            line, '\t', Idyn - 1, '\t', frame, '\t', i, '\t', j,
                            '\t', resn, '\t', resid, '\t', nim, '\t',
                            comp_1c[j][0], '\t', comp_1c[j][1]
                        )
                        line += 1
                    print(
                        "Resid", Lresid_1c[q], "is duplicated", Occ,
                        "times in a cluster for Idyn(dyn index)=", Idyn - 1,
                        ", If(frame index)=", frame, ", Ic(cluster index)=", i
                    )
                    raise Exception(
                        "Error! Please use lower Dcut! \n"
                    )
            Lresid_ac.append(Lresid_1c)

        # Check if resid dulplicate inside among clusters
        for i in range(0, len(clusters_1), 1):
            comp_2 = clusters_1[i]
            check2list = []
            for k in range(i + 1, len(clusters_1), 1):
                check2 = any(r in Lresid_ac[i] for r in Lresid_ac[k])
                check2list.append(check2)
            if True not in check2list:
                clusters_2.append(comp_2)

        # save frame number and tot lipids in cluster, and number of clusters,
        # and frame averaged cluster size
        Sc = sum(len(x) for x in clusters_2)
        Nc = len(clusters_2)
        return_info_1 = '%i\t%i\t%i\t%i\t%2.4f\n' % (
            Idyn - 1, frame, Sc, Nc, Sc / float(Nc)
        )
    return_info_2 = []
    if len(clusters_2) > 0 and frame >= 0:
        for i in range(0, len(clusters_2), 1):
            comp_f = clusters_2[i]
            # save frame number and number of clusters
            return_info_2.append('%s\t%s\t%s\t%s\n' % (Idyn-1, frame, i, len(comp_f)))
            for j in range(0, len(comp_f), 1):
                # Find line in the coordinate file
                resn, resid, nim = res_info(comp_f[j][0], comp_f[j][1], Data)
                # print line, frame, cluster, component and x y coordinate
                # in the cluster
                print(
                    "{0:<8d}{1:<8d}{2:<8d}{3:<8d}{4:<8d}{5:<8s}{6:<8s}{7:<8s}{8:<12f}{9:<12f}".format(
                        line, Idyn-1, frame, i, j, resn, resid, nim, comp_f[j][0], comp_f[j][1]
                    )
                )
                line += 1
    return return_info_1, return_info_2
