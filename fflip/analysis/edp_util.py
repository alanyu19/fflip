#!/usr/bin/python

from rflow.trajectory import *


def normalizecoor(trajectory, coordinates=2, com_selection=None, subselect="all",
              use_fixed_box_length=False, box_length_fixed=None):
    """
    Normalize the trajectory so that all coordinates are in [0,1] and the center of
    mass of the membrane is at 0.5. [copied and modified from rickflow]

    Args:
        trajectory:     An mdtraj trajectory object.
        coordinates:    0,1, or 2 (for x,y,z); or a list
        com_selection:  Selection of the membrane (to get the center of mass).
                        Can be a list of ints or a selection string.
        subselect:      Atom selection (usually the permeant). Can be a list of ints or a selection string.
                        The normalized array will only contain the atoms in this selection.

    Returns:
        np.array: The normalized coordinates.

    """
    if use_fixed_box_length:
        assert box_length_fixed!=None and (
                isinstance(box_length_fixed, float) or \
            isinstance(box_length_fixed, int)
        )

    membrane_center = center_of_mass_of_selection(trajectory, com_selection, coordinates)

    selected = select_atoms(trajectory, subselect)
    # normalize z coordinates: scale to [0,1] and shift membrane center to 0.5
    z_normalized = trajectory.xyz[:, selected,
                   coordinates].transpose() - membrane_center.transpose()

    if not use_fixed_box_length:
        z_normalized /= trajectory.unitcell_lengths[:, coordinates].transpose()
    else:
        assert (z_normalized < box_length_fixed).all()
        z_normalized /= box_length_fixed

    if com_selection is not None and len(com_selection) > 0:
        z_normalized += 0.5  # shift so that com is at 0.5

    z_normalized = np.mod(z_normalized, 1.0).transpose()

    return z_normalized


def find_resnames_from_psf(psf_file):
    res_list = []
    reading = False
    with open(psf_file, 'r') as psf:
        lines = psf.readlines()
        for l_number, l in enumerate(lines):
            if reading == True:
                _segid = l.strip().split()[1]
                _resid = l.strip().split()[2]
                _resname = l.strip().split()[3]
                _atmname = l.strip().split()[4]
                if _resname not in res_list: 
                    res_list.append(_resname)
            if '!NBOND' in lines[l_number+3] or '!NBOND' in lines[l_number+2] \
                    or '!NBOND' in lines[l_number+1]:
                break
            if '!NATOM' in l:
                reading = True
    return res_list


def find_atoms_from_psf(psf_file, resname):
    # constructing atom list from atlml file
    atom_list = []
    reading = False
    with open(psf_file, 'r') as psf:
        lines = psf.readlines()
        for l_number, l in enumerate(lines):
            if reading == True:
                _segid = l.strip().split()[1]
                _resid = l.strip().split()[2]
                _resname = l.strip().split()[3]
                _atmname = l.strip().split()[4]
                if resname.lower() == _resname.lower() and \
                        _atmname.lower() not in atom_list:
                    atom_list.append(_atmname.lower())
            if '!NBOND' in lines[l_number+3] or '!NBOND' in lines[l_number+2]\
                    or '!NBOND' in lines[l_number+1]:
                break
            if '!NATOM' in l:
                reading = True
    return atom_list


def gddradd(topology_file, traj_template, first, last, resn, atmn, bounds,
            axis=2, bins=500, verbose='v'):
    # Get Density Distribution by Resname and Atomname
    # units are in nanometer (could write a wrapper to
    # make full use of the simtk.unit)
    psf = CharmmPsfFile(topology_file)
    trajs = TrajectoryIterator(
        first_sequence=first, last_sequence=last,
        filename_template=traj_template, 
        topology_file=topology_file,
        atom_selection="all", load_function=md.load_dcd)
    """
    use vector for lb & up when your axis is list, e.g. calculating 2D/3D density
    """
    if not type(axis) == int:
        axis = list(axis)
    else:
        axis = list([axis])
    # find the bin size
    bin_size = 1
    for bindex, bds in enumerate(bounds):
        if type(bins) == int:
            bin_size = bin_size * (bds[1] - bds[0]) / bins
        else:
            bin_size = bin_size * (bds[1] - bds[0]) / list(bins)[bindex]
    omm_topology = md.Topology.from_openmm(psf.topology)
    selection = omm_topology.select(
        'name {} and resname {}'.format(atmn.upper(), resn.upper())
    )
    if len(selection) == 0:
        print("Error: No atom selected, please check your selection words")
    num_selection = len(list(selection))
    
    # find out the axis that are not used in histograming
    # and use it (them) for weights
    other_axis = []
    for i in [0, 1, 2]:
        if i not in axis:
            other_axis.append(i)
    total_frames = 0
    for i, traj in enumerate(trajs):
        if verbose == 'v':
            print('Calculating trajectory {}'.format(i + first))
        weights = np.ones(traj.n_frames)
        for ax in other_axis:
            weights = weights * (1 / traj.unitcell_lengths[:, ax])
        weights = weights / bin_size
        weights = weights.repeat(num_selection)
        onaxis = traj.xyz[:, :, axis]
        sample = onaxis[:, selection, :]
        sample = sample.reshape(-1, len(axis))
        total_frames += traj.n_frames
        h, edges = np.histogramdd(
            sample, bins=bins, range=bounds, weights=weights, density=None
        )
        if i == 0:
            H = h
        else:
            H = H + h
    return H / total_frames, edges


def gddra(topology_file, traj_template, first, last, resn, atmn, bounds,
          axis=2, bins=500, verbose='v'):
    # Get Density Distribution by Resname and Atomname
    # units are in nanometer (could write a wrapper
    # to make full use of the simtk.unit)
    psf = CharmmPsfFile(topology_file)
    trajs = TrajectoryIterator(
        first_sequence=first, last_sequence=last,
        filename_template=traj_template, 
        topology_file=topology_file,
        atom_selection="all", load_function=md.load_dcd)
    """
    use vector for lb & up when your axis is list, e.g. calculating 2/3D density
    """
    if not type(axis) == int:
        axis = list(axis)
    else:
        axis = list([axis])
    # find the bin size
    bin_size = 1
    for bindex, bds in enumerate(bounds):
        if type(bins) == int:
            bin_size = bin_size * (bds[1] - bds[0]) / bins
        else:
            bin_size = bin_size * (bds[1] - bds[0]) / list(bins)[bindex]
    omm_topology = md.Topology.from_openmm(psf.topology)
    if "\'" in atmn:
        # Do our own selection to avoid error in mdtraj for H3' and O3' in
        # sterols and possible future bad naming of atoms
        atmindxs = []
        with open(topology_file, 'r') as topfile:
            lines = topfile.readlines()
        for line in lines:
            if atmn.upper() in line:
                atmindxs.append(int(line.strip().split()[0]))
        selection = np.array(atmindxs)
        print(selection)
    else:
        selection = omm_topology.select(
            'name {} and resname {}'.format(atmn.upper(), resn.upper())
        )
    if len(selection) == 0:
        print("Error: No atom selected, please check your selection words")
    num_selection = len(list(selection))
    
    # find out the axis that are not used in histograming
    # and use it (them) for weights
    other_axis = []
    for i in [0, 1, 2]:
        if i not in axis:
            other_axis.append(i)
    total_frames = 0
    for i, traj in enumerate(trajs):
        if verbose == 'v':
            print('Calculating trajectory {}'.format(i + first))
        weights = np.ones(traj.n_frames)
        for ax in other_axis:
            weights = weights * (1 / traj.unitcell_lengths[:, ax])
        weights = weights / bin_size
        weights = weights.repeat(num_selection)
        onaxis = traj.xyz[:, :, axis]
        sample = onaxis[:, selection, :]
        sample = sample.reshape(-1, len(axis))
        total_frames += traj.n_frames
        weights = weights.reshape(-1, 1)
        print(sample.shape, weights.shape)
        h, edges = np.histogram(
            sample, bins=bins, range=bounds[0], weights=weights, density=None
        )
        if i == 0:
            H = h
        else:
            H = H + h
    return H / total_frames, edges


def get_atom_density_along_z_axis_for_residue(
        psf_file, traj_template, first, last, resname, bounds, axis, bins,
        save_to_where='.', verbose='v'
):
    # ONLY for jbk group EDP purpose
    binsize = 0.2
    resdir = os.path.join(save_to_where, 'atoms_{}_mdtraj'.format(resname.lower()))
    os.system('mkdir {}'.format(resdir))
    atoms = find_atoms_from_psf(psf_file, resname)
    for atm in atoms:
        H, e = gddra(
            psf_file, traj_template, first, last, resname, atm,
            bounds=[(-5, 5)], axis=2, bins=500, verbose=verbose
        )
        # unit conversion from nanometer to Angstrom
        data = np.swapaxes(np.array([10 * e[1:] - 0.1, H / 1000]), 0, 1)
        np.savetxt(resdir + '/{}_{}_{}.dat'.format(atm, first, last), data)


def combine_to_average(psf_file, path_to_data='.', z_unit_in_A=True, pop_drude=False):
    atom_folders = glob.glob(os.path.join(path_to_data, 'atoms_*'))
    residues = find_resnames_from_psf(psf_file)
    # compare the residue names got from psf_file
    # and those appeared in atoms_*_mdtraj folder
    resf = []
    for fn in atom_folders:
        resf.append(fn.split('_')[1])
    for r in residues:
        if r.lower() not in resf:
            print('Warning: density data for residue {} not found'.format(r))
    # averaging data for each res + atom
    for fn, resn in zip(atom_folders, resf):
        atoms_raw = find_atoms_from_psf(psf_file, resn)
        if pop_drude:
            atoms = []
            for _atom in atoms_raw:
                if _atom[0].upper() != 'D':
                    atoms.append(_atom)
                else:
                   pass
        else:
            atoms = atoms_raw
        zdata = np.loadtxt(os.path.join(fn, 'z.txt'))
        if z_unit_in_A:
            zdata = 10 * zdata
        for atom in atoms:
            data = []
            afiles = glob.glob(os.path.join(fn, 'blocks/{}_*[!_edp].dat'.format(
                atom)))
            for af in afiles:
                data.append(np.loadtxt(af))
            atom_average = np.mean(np.array(data), axis=0)
            to_save = np.array([zdata, atom_average]).transpose()
            np.savetxt(os.path.join(fn, '{}.dat'.format(atom)), to_save)
