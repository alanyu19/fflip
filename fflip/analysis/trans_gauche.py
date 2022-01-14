# -*- coding: utf-8 -*-


from rflow.trajectory import *
from fflip.ljpme.torsionfuncs import DihedralAngle


def get_all_dihedrals_in_a_tail(
    trajectories, topology, sn, first_carbon, last_carbon,
    verbose=False
):
    # Trajectories should be rickflow CharmmTrajectoryIterator, sn can be 2 and 3
    all_carbon = []
    for carbon in range(first_carbon, last_carbon - 2):
        if verbose:
            print('Carbon {} ...'.format(carbon), end = " ")
        data = []
        for trj in trajectories:
            dihe_ang = DihedralAngle(topology, 'C{}{}'.format(sn, carbon), 
                                      'C{}{}'.format(sn, carbon+1), 
                                      'C{}{}'.format(sn, carbon+2),
                                      'C{}{}'.format(sn, carbon+3))
            data.append(dihe_ang(trj))
        data = np.array(data) * 180/ np.pi
        all_carbon.append(data)
    if verbose:
        print("Finished current trajectory set.")
    return all_carbon


def tg_ratio(
    psf_file, traj_template, first_traj, last_traj,
    starting_carbons, ending_carbons, save_file=True,
    save_template={'sn2': 'sn2_tg_ratio_{}_{}.dat', 'sn1': 'sn1_tg_ratio_{}_{}.dat'}
):
    psf = CharmmPsfFile(psf_file)
    topology = md.Topology.from_openmm(psf.topology)
    
    trajs = CharmmTrajectoryIterator(
        first_sequence=first_traj, last_sequence=last_traj,
        filename_template=traj_template, 
        topology_file=psf_file,
        atom_selection="all", load_function=md.load_dcd
    )
    sn2 = get_all_dihedrals_in_a_tail(
        trajs, topology, 2, starting_carbons['sn2'], ending_carbons['sn2']
    )
    sn1 = get_all_dihedrals_in_a_tail(
        trajs, topology, 3, starting_carbons['sn1'], ending_carbons['sn1']
    )

    sn2 = np.array(sn2)
    sn1 = np.array(sn1)
    # sn1.shape == (# of carbons, # of trajectory files, # of frames per traj, # of lipids)

    gauche_sn2 = np.sum((120 >= sn2) & (sn2 > -120), axis = (0, 3))
    trans_sn2 = np.sum(120 < sn2, axis = (0, 3)) + np.sum((sn2 <= -120), axis = (0, 3))
    gauche_sn1 = np.sum((120 >= sn1) & (sn1 > -120), axis = (0, 3))
    trans_sn1 = np.sum(120 < sn1, axis = (0, 3)) + np.sum((sn1 <= -120), axis = (0, 3))
    # print(trans_sn1.shape)

    gauche_sn2 = gauche_sn2.flatten()
    trans_sn2 = trans_sn2.flatten()
    gauche_sn1 = gauche_sn1.flatten()
    trans_sn1 = trans_sn1.flatten()
    
    if save_file:
        np.savetxt(save_template['sn2'].format(first_traj, last_traj), trans_sn2/gauche_sn2)
        np.savetxt(save_template['sn1'].format(first_traj, last_traj), trans_sn1/gauche_sn1)

    return {'sn2': (trans_sn2/gauche_sn2).mean(), 'sn1': (trans_sn1/gauche_sn1).mean()}

