import pytest
import glob
import copy
import mdtraj as md
from simtk.openmm.app import CharmmPsfFile, CharmmParameterSet
from simtk.openmm.app import PME, LJPME, HBonds
import simtk.unit as u
from fflip.omm.util import filter_solution
from fflip.omm.playpara import brutal_nonbonded_parameter, change_nb_exceptions, find_nb_parameter
from fflip.ljpme.param import gen_param_offset


@pytest.fixture(scope="module")
def traj():
    return 0
    # return md.load(abspath('data/psm_dopc_ld.hdf5'))


@pytest.fixture(scope="function")
def system_72dppc_c36ljpme(psf_72dppc_c36):
    psf = psf_72dppc_c36
    parameter_files = glob.glob("data/toppar/*")
    parameters = CharmmParameterSet(*parameter_files)
    topology = md.Topology.from_openmm(psf.topology)
    psf.setBox(50, 50, 50)
    return psf.createSystem(
        parameters,
        nonbondedMethod=LJPME,
        constraints=HBonds,
        nonbondedCutoff=1.0 * u.nanometer,
        ewaldErrorTolerance=0.0005
    )


def test_parameterenergies(pc_c36, psf_72dppc_c36, system_72dppc_c36ljpme):
    parameter_sets = pc_c36.parse_groups()[:4]
    all_offsets = [gen_param_offset(ps, amount=0.01) for i, ps in enumerate(parameter_sets)]
    topology = md.Topology.from_openmm(psf_72dppc_c36.topology)
    old_parameters = []
    for g in parameter_sets:
        old_parameters.append(find_nb_parameter(system_72dppc_c36ljpme, topology, [g.center_names[0]], [g.par_type]))
    for g, offset in zip(parameter_sets, all_offsets):
        system_new, oldp, newp = brutal_nonbonded_parameter(
            system_72dppc_c36ljpme, topology, g, offset
        )
    change_nb_exceptions(psf_72dppc_c36, system_new)
    reference_parameters = [0.175, 0.2969, 0.4226, -0.485]
    for g, old_parameter, reference in zip(parameter_sets, old_parameters, reference_parameters):
        new_parameter = find_nb_parameter(system_new, topology, [g.center_names[0]], [g.par_type])
        # print(g, old_parameter, new_parameter)
        assert round(new_parameter[0], 4) == reference


def test_brutal_nonbonded_parameters_from_solution_file(pc_c36, psf_72dppc_c36, system_72dppc_c36ljpme):
    solution = filter_solution("data/solution_c36.txt")
    parameter_sets = pc_c36.parse_groups()
    all_offsets = [gen_param_offset(ps, amount=solution[i]) for i, ps in enumerate(parameter_sets)]
    topology = md.Topology.from_openmm(psf_72dppc_c36.topology)
    old_parameters = []
    for g in parameter_sets:
        old_parameters.append(find_nb_parameter(system_72dppc_c36ljpme, topology, [g.center_names[0]], [g.par_type]))
    for g, offset in zip(parameter_sets, all_offsets):
        system_new, oldp, newp = brutal_nonbonded_parameter(
            system_72dppc_c36ljpme, topology, g, offset
        )
    change_nb_exceptions(psf_72dppc_c36, system_new)
    reference_parameters = [0.1933, 0.2894, 0.4134, -0.4932, 0.9393, 0.2977, 0.4876, -0.66, -0.2181, 0.09, 0.0721]
    for g, old_parameter, reference in zip(parameter_sets, old_parameters, reference_parameters):
        new_parameter = find_nb_parameter(system_new, topology, [g.center_names[0]], [g.par_type])
        # print(g, old_parameter, new_parameter)
        assert round(new_parameter[0], 4) == reference
    

