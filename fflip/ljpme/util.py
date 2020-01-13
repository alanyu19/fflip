#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import glob
import pandas as pd


class ReweightingError(Exception):
    pass


class NoSimDataError(ReweightingError):
    def __init__(self):
        print('There is no simulated data, please run the reweight method first!')


class NoRewDataError(ReweightingError):
    def __init__(self):
        print('There is no  reweighted data, please run the reweight method first!')


class OtherReweightError(ReweightingError):
    def __int__(self, message):
        print(message)


class FutureResult(object):
    def get_result(self, path_to_file, file_name, time_to_wait, loading_method, delete = False):
        while not os.path.isfile(path_to_file + file_name):
            time.sleep(time_to_wait)
        if loading_method == 'numpy':
            result = np.genfromtxt(path_to_file + file_name)
        elif loading_method == 'pandas':
            result = pd.read_csv(path_to_file + file_name)
        else:
            with open(path_to_file + file_name, 'r') as file:
                result = file.readlines()
        if delete == True:
            os.system("rm -rf {}".format(path_to_file))
        return result



def move_traj(fromdir, todir):
    # Use with care
    os.system("mv {}/dyn?.dcd {}".format(fromdir, todir))
    os.system("mv {}/dyn??.dcd {}".format(fromdir, todir))
    os.system("mv {}/dyn???.dcd {}".format(fromdir, todir))


def construct_indexes(inp, properties):
    if inp=='all' or inp=='All':
        return list(range(len(properties)))
    elif ',' in inp:
        indexes = []
        values = inp.split(',')
        for v in values:
            assert 0 <= int(v) < len(properties)
            indexes.append(int(v))
        return indexes
    elif '-' in inp:
        first_last = inp.split('-')
        assert len(first_last) == 2
        first = int(first_last[0])
        last = int(first_last[1])
        assert len(properties) > last > first >= 0
        return list(range(first, last + 1))
    else:
        try:
            index = int(inp)
            assert 0 <= index < len(properties)
            return [index]
        except:
            raise Exception('Invalid property index(es)!')


def parse_first_last(inp):
    if ',' in inp:
        indexes = []
        values = inp.split(',')
        for v in values:
            indexes.append(int(v))
        return indexes
    else:
        try:
            index = int(inp)
            return [index]
        except:
            raise Exception('Invalid trajectory index(es)!')


def on_cluster(executable, executable_args_list, *args, **kwargs):
    out_dir = kwargs['out_dir']
    if os.path.isdir(out_dir):
        os.system("rm -rf {}".format(out_dir))
    os.system("mkdir {}".format(out_dir))
    if 'job_time' in kwargs:
        run_time = kwargs['job_time']
    else:
        run_time = '01:00:00'
    if 'partition' in kwargs:
        partition = kwargs['partition']
    else:
        partition = "ivy,sbr,hpodall,spodall"
    if 'ntasks' in kwargs:
        ntasks = kwargs['ntasks']
    else:
        ntasks = 1
    if 'conda_env' in kwargs:
        conda = kwargs['conda_env']
    else:
        conda = 'drude'
    with open(kwargs["submit_script"], 'w+') as f:
        f.write(
            "#!/bin/bash\n" + \
            "#SBATCH --output=./{}/{}.out\n".format(out_dir, kwargs["slurm_name"]) + \
            "#SBATCH --error=./{}/{}.err\n".format(out_dir, kwargs["slurm_name"]) + \
            "#SBATCH --time={}\n".format(run_time) + \
            "#SBATCH --partition={}\n".format(partition) + \
            "#SBATCH --ntasks={}\n".format(ntasks) + '\n'
        )
        f.write(
            "source ~/.bashrc\n" + \
            "conda activate {}\n".format(conda)
        )
        argstring = ''
        for arg in executable_args_list:
            argstring += ' {}'.format(str(arg))
        argstring += ' {}'.format(kwargs['out_dir'])
        f.write("\n" + "python " + executable + " " + argstring + " >& ./{}/{}.out".format(out_dir, kwargs["exec_name"]))
    os.system('sbatch {}'.format(kwargs['submit_script']))
    time.sleep(1)
    os.system("rm -f {}".format(kwargs['submit_script']))




class sensitivity_evaluator(object):
    def __init__(self, ngroups, exp_x, exp, sim_x, sim, rew, sens_type = 1, n_peaks = 2, n_foots = 1):
        """
        Args:
            -- exp: the experimental value(s)
            -- sim: the simulated value(s)
            -- rew: the reweighted value(s)
            -- sens_type: the catagory of the property (1: area/scd, 2: rdf)
        """
        self.ngroups = ngroups
        self.sim = sim
        self.exp = exp
        if sens_type == 2:
            self.exp_x = exp_x
            self.sim_x = sim_x
        self.rew = rew
        self.sens_type = sens_type
        self.n_peaks = n_peaks
        self.n_foots = n_foots

    @property
    def diff_sim_exp(self):
        if self.sens_type == 1:
            """
            Area / Order parameter
            """
            return self.sim - self.exp
        elif self.sens_type == 2:
            """
            Things like rdf which contain both positions and magnitudes
            """
            x_exp, y_exp = find_rdf_peaks_and_foots(
                self.exp_x, self.exp, first_n_peaks = self.n_peaks, first_n_foots = self.n_foots, smooth_window_size = 3
            )
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim_x, self.sim, first_n_peaks = self.n_peaks, first_n_foots = self.n_foots, smooth_window_size = 1
            )
            return x_sim - x_exp, y_sim - y_exp

    @property
    def rel_diff_sim_exp(self):
        if self.sens_type == 1:
            """
            Area / Order parameter
            """
            return (self.sim - self.exp) /self.exp
        elif self.sens_type == 2:
            """
            Things like rdf which contain both positions and magnitudes
            """
            x_exp, y_exp = find_rdf_peaks_and_foots(
                self.exp_x, self.exp, first_n_peaks = self.n_peaks, first_n_foots = self.n_foots, smooth_window_size = 3
            )
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim_x, self.sim, first_n_peaks = self.n_peaks, first_n_foots = self.n_foots, smooth_window_size = 1
            )
            return (x_sim - x_exp)/x_exp, (y_sim - y_exp)/y_exp

    @property
    def diff_rew_sim(self):
        if self.sens_type == 1:
            rew = self.rew
            sim_tiled = np.tile(self.sim, self.ngroups)
            return rew - sim_tiled
        elif self.sens_type == 2:
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim_x, self.sim, first_n_peaks=self.n_peaks, first_n_foots=self.n_foots, smooth_window_size=1
            )
            r_list, peak_foot_value_list = find_rdf_peaks_and_foots(
                self.sim_x, self.rew, first_n_peaks = self.n_peaks, first_n_foots = self.n_foots, smooth_window_size = 1
            )
            diff_list = []
            for i, (r, pfv) in enumerate(zip(r_list, peak_foot_value_list)):
                diff_list.append(np.array([r - x_sim, pfv - y_sim]))
            return np.array(diff_list)

    @property
    def rel_diff_rew_sim(self):
        if self.sens_type == 1:
            rew = self.rew
            sim_tiled = np.tile(self.sim, self.ngroups)
            return (rew - sim_tiled) /self.exp
        elif self.sens_type == 2:
            x_sim, y_sim = find_rdf_peaks_and_foots(
                self.sim_x, self.sim, first_n_peaks=self.n_peaks, first_n_foots=self.n_foots, smooth_window_size=1
            )
            r_list, peak_foot_value_list = find_rdf_peaks_and_foots(
                self.sim_x, self.rew, first_n_peaks = self.n_peaks, first_n_foots = self.n_foots, smooth_window_size = 1
            )
            x_exp, y_exp = find_rdf_peaks_and_foots(
                self.exp_x, self.exp, first_n_peaks=self.n_peaks, first_n_foots=self.n_foots, smooth_window_size=3
            )
            rel_diff_list = []
            for i, (r, pfv) in enumerate(zip(r_list, peak_foot_value_list)):
                rel_diff_list.append(np.array([(r - x_sim) /x_exp, (pfv - y_sim) /y_exp]))
            return np.array(rel_diff_list)
