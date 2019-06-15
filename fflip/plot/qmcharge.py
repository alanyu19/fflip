#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

class charge_compare(object):
    def __init__(self, atoms, mm_charges, qm_fitting_methods, qm_data_path='./'):
        self.atoms = atoms
        self.mm_charges = mm_charges
        self.qm_fitting_methods = qm_fitting_methods
        self.qm_data_path = qm_data_path
    
    @property
    def qm_dict(self, suffix = ['{}_avg.dat', '{}_err.dat']):
        dic = {}
        for i, method in enumerate(self.qm_fitting_methods):
            avg = np.loadtxt(self.qm_data_path + suffix[0].format(method))
            err = np.loadtxt(self.qm_data_path + suffix[1].format(method))
            dic[method] = [avg, err]
        return dic
        
    def gen_table(self, ff = 'C36'):
        tdic = {}
        tdic['atom'] = self.atoms
        tdic[ff] = self.mm_charges
        for method in self.qm_dict.keys():
            tdic[method] = self.qm_dict[method][0]
            tdic[method + '_err'] = self.qm_dict[method][1]
        df = pd.DataFrame(tdic)
        print(df)
        df.to_csv('qm-charges.csv')
        
    def plot(self, group_boundaries = [0, 16, 27, 35, 45], 
             fig_name = './DPPC-charge.png', ff = 'C36', 
             colors = ['deeppink', 'royalblue', 'gold',], dpi = 300):
        # the group_boundaries here is for dppc, same for fig_name
        font = {'fontname': 'sans-serif'}
          # can add more
        xl = []       
        f, a = plt.subplots(len(group_boundaries) - 1, 1, 
                            figsize=(12, len(group_boundaries) * 9))        
        for i in range(len(self.atoms)):
            xl.append(i)
        
        lowers = group_boundaries[:-1]
        uppers = group_boundaries[1:]
        
        for k in range(len(lowers)):
            
            ax = plt.subplot(len(lowers), 1, k+1)
            
            for j, method in enumerate(self.qm_dict.keys()):
                # avg
                avgs = self.qm_dict[method][0][lowers[k]:uppers[k]]
                # err
                errs = self.qm_dict[method][1][lowers[k]:uppers[k]]
                ax.errorbar(xl[lowers[k]:uppers[k]], avgs, yerr = errs, 
                             fmt = 'o-', markersize = 6, linewidth = 2, 
                             capsize = 4, capthick = 4, color = colors[j], 
                             label = self.qm_fitting_methods[j])
                
            ax.errorbar(xl[lowers[k]:uppers[k]], self.mm_charges[lowers[k]:uppers[k]], 
                         fmt = 'o', markersize = 6, linewidth = 2, 
                         capsize = 4, capthick = 4, color = 'black',
                         label = ff)
            
            ax.legend(fontsize = 18, loc = 1)
            
            ax.spines['bottom'].set_linewidth(1.2)
            ax.spines['top'].set_linewidth(1.2)
            ax.spines['left'].set_linewidth(1.2)
            ax.spines['right'].set_linewidth(1.2)
            
            plt.ylabel("Partial Charge [e]", fontsize = 22, **font, rotation = 90,)
            plt.yticks(**font, fontsize = 22)
            #plt.ylim((0.2, 0.4))
            plt.xticks(xl[lowers[k]:uppers[k]], self.atoms[lowers[k]:uppers[k]], 
                       fontsize = 22, rotation = 90)

        plt.subplots_adjust(hspace = 0.25)
        plt.savefig(fig_name, dpi = dpi)
   