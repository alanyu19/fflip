import numpy as np
from matplotlib import pyplot as plt


# one = [[0.34,  1,   180.0], [0.05,  2,     0.0], [1.94,  3,     0.0],
#       [0.16,  4,   180.0], [0.05,  5,   180.0], [0.03,  6,     0.0]]
#
# dict_1 = {'name': 'use_c36', 'para': one}
#
# two = [[0.10,  1,   180.0], [0.16,  2,     0.0], [1.93,  3,     0.0],
#       [0.06,  4,   180.0], [0.04,  5,   180.0], [0.05,  6,     0.0]]
#
# dict_2 = {'name': 'use_qm', 'para': two}
#
# dicts = [dict_1, dict_2]


def plot_charmm_dihedral_cosine_series(input_dictionaries, file_to_save=0):
    phi = np.arange(0, 2 * np.pi, 2 * np.pi/200)
    
    plot_list = []
    legend = []        
    for count, one_dict in enumerate(input_dictionaries):
        plot_list.append(0 * (1 + np.cos(0)))
        # use_c36 = 0 * (1 + np.cos(0))
        for m in one_dict['para']:
            plot_list[count] += m[0] * (1 + np.cos(m[1]*phi-np.pi*m[2]/180))
        legend.append(one_dict['name'])
        
    fig, ax = plt.subplots(figsize=(8, 5))       
    for count, one_plot in enumerate(plot_list):          
        ax.plot(phi/np.pi*180, one_plot)
        #ax.plot(phi/np.pi*180, use_qm)
    ax.legend(legend, fontsize=14, loc='upper right')
    plt.xlabel("phi ($^\circ$)", fontsize=24, )
    plt.ylabel("energy (kcal/mol)", fontsize=24, )
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    if file_to_save:
        plt.savefig(file_to_save, format='png', 
                    dpi=300, bbox_inches = 'tight')
    
#plot_charmm_dihedral_cosine_series(dicts)
