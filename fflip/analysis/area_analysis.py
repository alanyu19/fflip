# -*- coding: utf-8 -*-
from fflip.plot.area_plot import *
from fflip.analysis.util import *
import numpy as np

def get_area_and_ka(
    nlipid, temperature, block_size, step_size, force_calc=False, force_fraction=0.6, plot=True,
    file_to_read='area.dat', data_col=-1, time_to_skip=10, plot_title='', plot_only=False
):
    if plot_only:
        data0 = np.loadtxt(file_to_read)
        data0 = 100 * data0
        file_to_save = 'sa.png'
        y_max = int(data0.max() + (data0.max() - data0.min())*0.3)
        y_min = int(data0.min() - (data0.max() - data0.min())*0.3)
        plot_area(file_to_save, data0, skip_ns=0, interval=step_size, area_range=(y_min, y_max), title=plot_title)
        return
    data0 = np.loadtxt(file_to_read)
    data_dimension = len(list(data0.shape))
    if data_dimension == 1:
        pass
    elif data_dimension == 2:
        print('dimension is two')
        data0 = data0[:,data_col]
    if data0[0] > 50:
        # A^2
        multiply_factor = 1
    elif 1 > data0[0] > 0.3:
        # nm^2
        multiply_factor = 100
    data0 = data0 * multiply_factor
    kb = 1.3806
    sad, data_equil, eq_start_at, use_last_steps, block_size = \
    get_equil_data(data0, step_size=step_size, nskip=500, block_size=block_size)
    print("Equilibrium reached after {} ns".format(eq_start_at))
    print("Block size is {} ns".format(block_size))
    if use_last_steps < 3 * block_size * step_size or force_calc == True: 
        print("Ignoring the staring point of euiqlibrium above, using last {}% data ...".format(int(force_fraction*100)))
        print("\n=.=")
        total_length = data0.shape[0]
        laststeps = int(total_length * force_fraction)
        data_to_use = data0[-laststeps:]
        eq_start_at = int((1 - force_fraction) * (total_length/step_size)) # refresh for plotting
    else:
        data_to_use = data_equil
    mean, err, bn, skip = calc_block_avg(data_to_use, step_size, block_size, time_to_skip=time_to_skip)
    # Calculate Compressibility Modulus using largest possible data set
    N_steps= np.shape(data_to_use)[0]     # number of data points
    print("Using last {} ns data for compressibility calculation ...".format(round(N_steps*step_size,2)))
    mean2, err2, bn2, skip2 = calc_block_avg(data_to_use, step_size, block_size, time_to_skip=0, verbose = False)
    VarS=0      # VarS is variance after the system reaches equilibrium
    for i in range(0, N_steps, 1):
        VarS=VarS+(data_to_use[i] - mean2)**2
    Variance=VarS/N_steps
    Ka = kb * temperature * mean / (nlipid * Variance * 1000) # order is corrected/recovered using 1000
    print('Final SA_avg/lip and err: %3.2f +- %3.2f' % (mean, err))
    print('Compressibility modulus is: {0:.2f}'.format(Ka))
    # plot module is in fflip.plot
    if plot:
        # skip_frames = int(eq_start_at/step_size)
        file_to_save = 'sa.png'
        y_max = int(data0.max() + (data0.max() - data0.min())*0.3)
        y_min = int(data0.min() - (data0.max() - data0.min())*0.3)
        plot_area(file_to_save, data0, skip_ns=eq_start_at, interval=step_size, area_range=(y_min, y_max), title=plot_title)
