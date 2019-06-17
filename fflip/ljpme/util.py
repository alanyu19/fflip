# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

def reweight(property_name, folders, ginfo_file = "./gtcnp_order.txt", result_dir = "./rew", temperature=323.15, decimal_places=3):
    # write this into a class with a name
    if not os.path.isdir(result_dir):
        os.system("mkdir {}".format(result_dir))
    rt = reweightProperty(temperature=temperature)
    sensitivity = []
    value = []
    for f in folders:
        property_data = np.loadtxt("{}/{}.dat".format(f, property_name))
        if property_name == 'area' or property_name == 'areas':
            # maybe we can do this in rickflow
            property_data = property_data * 100
        o_energy_data = np.loadtxt("{}/original.dat".format(f))
        p_energy_data = np.loadtxt("{}/perturbed.dat".format(f))
        if len(p_energy_data.shape) == 1:
            p_energy_data = np.reshape(p_energy_data, (-1,1))
        old, new = rt(property_data, o_energy_data, p_energy_data)
        fl = f.split("/")[-1]
        np.savetxt(result_dir + '/org_{}_avg_{}.dat'.format(property_name, fl), np.round(np.reshape(old, (1)), decimal_places))
        np.savetxt(result_dir + '/rew_{}_avg_{}.dat'.format(property_name, fl), np.round(np.reshape(new, (-1,1)), decimal_places))
        np.savetxt(result_dir + '/sensitivity_{}_{}.dat'.format(property_name, fl), np.round(np.reshape((new-old), (-1,1)), decimal_places))
        sensitivity.append(new - old)
        value.append(new)
    sensitivity_avg = np.round(np.mean(np.array(sensitivity), axis=0), decimal_places)
    sensitivity_std = np.round(np.std(np.array(sensitivity), axis=0), decimal_places)
    value_avg = np.round(np.mean(np.array(value), axis=0), decimal_places)
    value_std = np.round(np.std(np.array(value), axis=0), decimal_places)
    format_string = "{0} {1:>8." + str(decimal_places) + "f} +- {2:<5." + str(decimal_places) + "f} {3:>8." + str(decimal_places) + "f} +- {4:<5." + str(decimal_places) + "f}\n"
    with open(result_dir + "/all-{}-sens.txt".format(property_name), 'w') as ftw:
        with open(ginfo_file) as ftr:
            glines = ftr.readlines()
            sensitivity_avg = np.reshape(sensitivity_avg, (-1))
            sensitivity_std = np.reshape(sensitivity_std, (-1))
            value_avg = np.reshape(value_avg, (-1))
            value_std = np.reshape(value_std, (-1))
            assert len(glines) == sensitivity_avg.shape[0] == sensitivity_std.shape[0]
            for j in range(len(glines)):
                ginfo = glines[j].rstrip()
                sens_avg = sensitivity_avg[j]
                sens_std = sensitivity_std[j]
                val_avg = value_avg[j]
                val_std = value_std[j]
                ftw.write(format_string.format(ginfo, sens_avg, sens_std, val_avg, val_std))
                #ftw.write("{0} {1:>8.3f} +- {2:<5.3f}\n".format(ginfo, sens_avg, sens_std))
    # return
    return sensitivity_avg, sensitivity_std, value_avg, value_std


def create_dict(area_weight = 1000, scd_weight = 1):
    # read atom names and types of parameters:
    para_names = [np.nan, np.nan] # two place holders, same for the nexts
    para_types = [np.nan, np.nan]
    para_changes = [np.nan, np.nan]

    with open('./rew/all-area-sens.txt', 'r') as f:
        for line in f.readlines():
            line = line.strip().split()
            nm = line[0]; para_names.append(nm)
            tp = line[1]; para_types.append(tp)
            cg = line[2]; para_changes.append(cg)


    # initialize the dictionary:
    prop_dict = {}

    # fill in first few columns and create key for area
    prop_dict['name'] = para_names
    prop_dict['type'] = para_types
    # weight is always the first and the second should be deviation from the experiment
    # sim area is 58.6556 and the exp is 63.0
    area_dev = 63.0 - 58.6556
    prop_dict['area'] = [area_weight, area_dev]

    # read scd_diffs between exp and sim:
    scd_info = pd.read_csv("scd_table.csv")
    scd_name = scd_info['name']
    scd_diff = scd_info['diff']
    for i, key in enumerate(scd_name):
        prop_dict[key] = [scd_weight, round(scd_diff[i], 3)]

    # create key for stderr of sensitivity and gaussian_charge - current(during optimization) charge
    # the second should be activated after finding the QM charges
    prop_dict['sens err'] = [1, None] # no deviation available, not going to use the standard for area and scd
    # (to be activated) prop_dict['qm-mm'] = [1, None] # no deviation available, not going to use the standard for area and scd
     

    # read in sensitivity information:
    # 1. area/lipid + the stderr of the robustness of it
    with open('./rew/all-area-sens.txt', 'r') as f:
        for name_cmp, line in zip(prop_dict['name'][2:], f.readlines()):
            line = line.strip().split()
            assert line[0] == name_cmp # I think this is strong enough without the para_type checking (YYL)
            prop_dict['area'].append(round(float(line[3])/1, 4)) # thefactor devided here should have no influce, but be careful anyway
            prop_dict['sens err'].append(np.abs(round(float(line[5])/float(line[3]), 4)))

    # 2. scd
    def fill_scd_info(name):
        with open('./rew/all-{}-sens.txt'.format(name), 'r') as f:
            # loop over parameters
            for name_cmp, line in zip(prop_dict['name'][2:], f.readlines()):
                line = line.strip().split()
                assert line[0] == name_cmp # I think this is strong enough without the para_type checking (YYL)
                # thefactor devided here should have no influce, but be careful anyway
                prop_dict[name].append(round(float(line[3])/1, 4))
    # 2.1 scd for special C22
    fill_scd_info('scd-c22-h2s')
    fill_scd_info('scd-c22-h2r')
    # 2.2 scd for other sn-2 carbons
    for carbon in range(2, 16):
        fill_scd_info('scd-c3{}'.format(carbon))
    # 2.3 scd for sn-1 carbons
    for carbon in range(3, 16):
        fill_scd_info('scd-c2{}'.format(carbon))

    return prop_dict

