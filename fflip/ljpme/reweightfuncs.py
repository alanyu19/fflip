# -*- coding: utf-8 -*-

import numpy as np

from rflow.observables import TimeSeries
from fflip.omm.util import beta_kjmol


class DataLoader(object):
    def __init__(self, dimension=1, load_func=np.loadtxt):
        self.dimension = dimension
        self.load_func = load_func

    def __call__(self, *args, **kwargs):
        """
        - kwargs: must include file_to_load
        """
        data = self.load_func(kwargs['file_to_load'])
        # reshape to match other data
        return data


class CoreReweighter(object):
    def __init__(self, **kwargs):
        """
        The reweighter to do weightings
        Args:
            **kwargs:
                - temperature: defaul=323.15, the original temperature,
                if temperature_perturbed is not provided, it's also used as
                the temperature of the reweighted state.
                - temperature_perturbed: the temperature of the reweighted state
        """
        if "temperature" not in kwargs:
            self.beta = beta_kjmol(323.15)
        else:
            self.beta = beta_kjmol(float(kwargs["temperature"]))
        if "temperature_perturbed" not in kwargs:
            self.beta2 = self.beta
        else:
            self.beta2 = beta_kjmol(float(kwargs["temperature_perturbed"]))

    def __call__(self, property_data, o_energy_data, p_energy_data):
        """
        One of the property_data and perturbed_energy data can be 2D array.
        Args:
            property_data: the property data array, can be 1D (area, scd ...)
            or 2D (RDF, each row is a snapshot, columns for different distance.
            o_energy_data: the original energy comes with the property.
            p_energy_data: the perturbed energy, can be 1D or 2D array.
            If a 2D array is provided, it calculates the reweighted property
            for different perturbations of parameter.
        Returns:
            The original property, the perturbed property
        """
        assert \
            (len(property_data) == len(o_energy_data) == len(p_energy_data)), \
            "Data are of different size!"
        o_energy_data = np.array(o_energy_data)
        p_energy_data = np.array(p_energy_data)
        property_data = np.array(property_data)
        # if all one dimension
        if len(p_energy_data.shape) == 1 and len(property_data.shape) == 1:
            tune = self.beta2 * np.mean(p_energy_data) - \
                   self.beta * np.mean(o_energy_data)
            obs = np.sum(
                property_data * np.exp(
                    -self.beta2 * p_energy_data +
                    self.beta * o_energy_data + tune
                )
            )
            ptf = np.sum(
                np.exp(
                    -self.beta2 * p_energy_data +
                    self.beta * o_energy_data + tune
                )
            )
            original = np.mean(property_data)
            perturbed = round(obs/ptf, 6)
            return original, perturbed
        elif len(p_energy_data.shape) == 1 and len(property_data.shape) == 2:
            # TODO: This part is not checked / tested
            ncol = property_data.shape[1]
            # nrow = property_data.shape[0]
            tune = self.beta2 * np.mean(p_energy_data) - \
                self.beta * np.mean(o_energy_data)
            p_energy_data = np.tile(
                np.reshape(p_energy_data, (-1, 1)), (1, ncol)
            )
            o_energy_data = np.tile(
                np.reshape(o_energy_data, (-1, 1)), (1, ncol)
            )
            obs = np.sum(property_data * np.exp(
                -self.beta2 * p_energy_data + self.beta * o_energy_data + tune
            ), axis=0)
            ptf = np.sum(np.exp(
                -self.beta2 * p_energy_data + self.beta * o_energy_data + tune
            ), axis=0)
            original = np.mean(property_data, axis=0)
            perturbed = obs/ptf
            return original, perturbed
        elif len(p_energy_data.shape) == 2 and len(property_data.shape) == 1:
            original = 0
            perturbed = []
            for column in range(p_energy_data.shape[1]):
                p_energy_data_column = p_energy_data[:, column]
                tune = self.beta2 * np.mean(p_energy_data_column) - \
                    self.beta * np.mean(o_energy_data)
                obs = np.sum(property_data * np.exp(
                    -self.beta2 * p_energy_data_column +
                    self.beta * o_energy_data + tune
                ))
                ptf = np.sum(np.exp(
                    -self.beta2 * p_energy_data_column +
                    self.beta * o_energy_data + tune
                ))
                original = np.mean(property_data)
                perturbed_column = round(obs/ptf, 5)
                perturbed.append(perturbed_column)
            perturbed = np.array(perturbed)
            original = round(float(original), 6)
            return original, perturbed


class SimpleReweighter(object):
    def __init__(self, **kwargs):
        import warnings
        warnings.warn("This object will be depracated in the future version, "
                      "please adapt your code to the new CoreReweighter!")
        if "temperature" not in kwargs:
            self.beta = beta_kjmol(323.15)
        else:
            self.beta = beta_kjmol(float(kwargs["temperature"]))
        if "temperature_perturbed" not in kwargs:
            self.beta2 = self.beta
        else:
            self.beta2 = beta_kjmol(float(kwargs["temperature_perturbed"]))

    def __call__(self, property_data, o_energy_data, p_energy_data):
        """the p_energy_data can be both 1D or 2D"""
        o_energy_data = np.array(o_energy_data)
        p_energy_data = np.array(p_energy_data)
        property_data = np.array(property_data)
        assert len(property_data) == len(o_energy_data) == len(p_energy_data), \
            "Data are of different size!"
        original = 0
        if p_energy_data.shape[1] == 1:
            tune = self.beta2 * np.mean(p_energy_data) - self.beta * np.mean(
                o_energy_data
            )
            obs = np.sum(
                property_data * np.exp(
                    -self.beta2 * p_energy_data[:, 0] +
                    self.beta * o_energy_data + tune
                )
            )
            ptf = np.sum(
                np.exp(
                    -self.beta2 * p_energy_data[:, 0] +
                    self.beta * o_energy_data + tune
                )
            )
            original = np.mean(property_data)
            perturbed = round(obs/ptf, 5)
        else:
            perturbed = []
            for column in range(p_energy_data.shape[1]):
                p_energy_data_column = p_energy_data[:, column]
                tune = self.beta2 * np.mean(p_energy_data_column) - \
                    self.beta * np.mean(o_energy_data)
                obs = np.sum(
                    property_data * np.exp(
                        -self.beta2 * p_energy_data_column +
                        self.beta * o_energy_data + tune
                    )
                )
                ptf = np.sum(
                    np.exp(
                        -self.beta2 * p_energy_data_column +
                        self.beta * o_energy_data + tune
                    )
                )
                original = np.mean(property_data)
                perturbed_column = round(obs/ptf, 5)
                perturbed.append(perturbed_column)
            perturbed = np.array(perturbed)
        ''' Return both the original and perturbed '''
        original = round(float(original), 5)
        return original, perturbed


def reweight_many_params(param_ids, temperature,
                         path_to_prp, path_to_org, path_to_ptb,
                         starting_traj, ending_traj,
                         traj_interval_energy, traj_interval_prop=None):
    """
    A function to reweight property for series of parameters,
    the input energy should be in kj/mole
    Args:
        param_ids: int, indexers in the observable/potential files
        temperature: float, the temperature, currently support only one
        path_to_prp: str
        path_to_org: str
        path_to_ptb: str
        starting_traj: int, the first trajectory index
        ending_traj: int, the last trajectory index
        traj_interval_energy: int, the block size in terms of trajectory count
        traj_interval_prop: int, the block size in terms of trajectory count

    Returns:
        the original and a list (for different parameters) of the perturbed
        observable.
    """
    # TODO: add a file checking function here to make sure that all the files
    #  needed are available
    if traj_interval_prop is None:
        traj_interval_prop = traj_interval_energy
    old = 0
    new_list = []
    for par in param_ids:
        rt = CoreReweighter(temperature=temperature)
        property_series = TimeSeries(evaluator=DataLoader())
        original_series = TimeSeries(evaluator=DataLoader())
        perturbed_series = TimeSeries(evaluator=DataLoader())
        for trj in range(starting_traj, ending_traj + 1, traj_interval_energy):
            original_series(file_to_load=path_to_org.format(trj))
            perturbed_series(file_to_load=path_to_ptb.format(trj, str(par)))
        for trj in range(starting_traj, ending_traj + 1, traj_interval_prop):
            property_series(file_to_load=path_to_prp.format(trj))
        old, new = rt(
            property_series.data, original_series.data, perturbed_series.data
        )
        new_list.append(new)
    return old, new_list


def reweight_one_param(prop, o_energies, p_energies, obeta, pbeta):
    """
    Get reweighted properties based on original and perturbed energies.
    Input data should be np array of same size.
    """
    assert (len(prop) == len(o_energies) == len(p_energies)),\
        "Data are of different size!"
    tune = pbeta * np.mean(p_energies) - obeta * np.mean(o_energies)
    obs = np.sum(prop * np.exp(-pbeta * p_energies + obeta * o_energies + tune))
    ptf = np.sum(np.exp(-pbeta * p_energies + obeta * o_energies + tune))
    original = np.mean(prop)
    perturbed = obs/ptf
    
    """Return both the original and perturbed, plus the average energy shift"""
    return original, perturbed, tune 
