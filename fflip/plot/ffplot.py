import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import minimize


def find_distance_to_curve(x, f, xp, fp):
    """ a function that can be used to find the distance 
    of points to a 'curve', the curve should contain enough
    data so that data can be interpolated linearly. """
    fake = np.interp(x, xp, fp)
    return np.sum(np.abs(f - fake))


class distance(object):
    def __init__(self, exp_q, exp_f, sim_q, sim_f):
        self.exp_q = exp_q
        self.exp_f = exp_f
        self.sim_q = sim_q
        self.sim_f = sim_f
    def __call__(self, x):
        scaled_sim_f = x * self.sim_f
        return find_distance_to_curve(self.exp_q, self.exp_f, 
                                      self.sim_q, scaled_sim_f)


# pass the distance object to the optimizer
def scale_fq_data_on(distance):
    # starting point is set to be 1
    print("Starting Optimization ...")
    res = minimize(distance, 1, method='nelder-mead',
                   options={'xtol': 1e-6, 'disp': True})
    return res.x


def get_ffs(sim_dict_list, exp_dict_list, name_of_lipid, ff_type, 
            exp_to_scale_sim_on = 1, scale_lower_cutoff = 0, save_fig = True):
    """ The dictionary should contain the following information:
    1. the data (should be in shape n*2
    2. the legend
    3. color ? """
    font = {'fontname':'serif', 'weight':'bold'}
    sim_color = ['red']
    exp_color = ['wheat', 'lightblue', 'lightgray']
    # First thing is to get the reference experimental data
    ref_q = exp_dict_list[exp_to_scale_sim_on]['data'][:, 0]
    ref_f = exp_dict_list[exp_to_scale_sim_on]['data'][:, 1]
    ref_q = ref_q[np.where(ref_f > scale_lower_cutoff)]
    ref_f = ref_f[np.where(ref_f > scale_lower_cutoff)]
    assert ref_q.shape[0] == ref_f.shape[0]
    #sim_q = [], sim_f = [], exp_q = [], exp_f = []
    average_distances = []
    for sim_dict in sim_dict_list:
        #sim_q.append(sim_dict['data'][:0])
        #sim_f.append(sim_dict['data'][:1])
        dist = distance(ref_q, ref_f, sim_dict['data'][:, 0], sim_dict['data'][:, 1])
        scaling_factor = scale_fq_data_on(dist)
        print("Scaling factor for {} is {}".format(sim_dict['legend'], scaling_factor))
        sim_dict['data'][:, 1] = sim_dict['data'][:, 1] * scaling_factor
        average_distances.append(dist(scaling_factor)/np.shape(ref_q)[0])

    font = {'fontname':'serif', 'weight':'bold'}

    fig, ax = plt.subplots(figsize=(12, 10))
    plt.xlim((0, 1))
    plt.ylim((0, 2.5))
    plt.xlabel("q", fontsize=24, **font,)
    plt.ylabel("F(q)", fontsize=24, **font, rotation=90,)
    plt.xticks(np.arange(0, 1, 0.1), **font)
    plt.yticks(np.arange(0, 2.5, 0.5), **font)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2) #set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    plt.title("{}".format(name_of_lipid), fontsize=36, y=0.9, fontdict= {'verticalalignment': 'baseline'},)

    dx = -1/7
    dy = -1/7
    offset = matplotlib.transforms.ScaledTranslation(dx, 0, fig.dpi_scale_trans)
    plt.setp(ax.yaxis.get_majorticklabels())
    for label in ax.yaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    offset = matplotlib.transforms.ScaledTranslation(0, dy, fig.dpi_scale_trans)
    plt.setp(ax.xaxis.get_majorticklabels())
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

    for i, ff in enumerate(sim_dict_list):
        ax.plot(ff['data'][:, 0], ff['data'][:, 1], '-', linewidth=3,
                label=ff['legend'], markersize = 4, color = sim_color[i])
        ax.legend(loc='upper right', fontsize =20, handlelength=2)
    for i, ff in enumerate(exp_dict_list):
        ax.plot(ff['data'][:, 0], ff['data'][:, 1], 'o', linewidth=3,
                label=ff['legend'], markersize = 4, color = exp_color[i])
        ax.legend(loc='upper right', fontsize=20, handlelength=2)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(24)
        #tick.label.set_rotation('vertical')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(24)
        #tick.label.set_rotation('vertical')
    #ax.legend.bbox_transform.transform_angles(90)
    if save_fig:
        plt.savefig('{}-{}.png'.format(name_of_lipid, ff_type), dpi = 150, 
                    transparent = False, bbox_inches = 'tight')
    plt.show(transparent = True) #wangyihang zhenshuai
    print(average_distances)


#data = np.genfromtxt("./xopc/popc/sim_to_exp_output/popc.xff", skip_header=1)
#sim_dict1 = {'data': data, 'legend': 'fucked'}
#data = np.genfromtxt("./xopc/experimental_data/POPC_ORI@30Cin0D2006.xff", skip_header=9)
#exp_dict1 = {'data': data, 'legend': 'fucker_1'}
#data = np.genfromtxt("./xopc/experimental_data/POPC_ULV@30Cin0D.xff", skip_header=9)
#exp_dict2 = {'data': data, 'legend': 'fucker_2'}
#data = np.genfromtxt("./xopc/experimental_data/DOPCexpXulv.xff", skip_header=9)
#exp_dict3 = {'data': data, 'legend': 'fucker_3'}
#
#get_ffs([sim_dict1], [exp_dict1, exp_dict2, exp_dict3], 'POPC', 'xff',
#        exp_to_scale_sim_on = 1, scale_lower_cutoff = 1.5)