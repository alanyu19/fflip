# -*- coding: utf-8 -*-
from fflip.plot import area_plot

def get_area_and_ka(nlipid, temperature, block_size, step_size, force_calc = False, force_fraction = 0.6, plot = True,
        file_to_read = 'area.dat', data_dimension = 1, data_col = -1, multiply_factor = 100, time_to_skip = 10):
    data0 = np.loadtxt(file_to_read)
    if data_dimension == 1:
        pass
    elif data_dimension == 2:
        data0 = data[:,data_col]
    data0 = data0 * multiply_factor
    kb = 1.3806
    sad, data_equil, eq_start_at, use_last_steps, block_size = \
    get_equil_data(data0, step_size = step_size, nskip = 500, block_size = block_size)
    print("Equilibrium reached after {} ns".format(eq_start_at))
    print("Block size is {} ns".format(block_size))
    if use_last_steps < 3 * block_size * step_size or force_calc == True: 
        print("Ignoring the staring point of euiqlibrium above, using last {}% data ...".format(int(force_fraction*100)))
        print("\n=.=")
        total_length = data0.shape[0]
        laststeps = int(total_length * force_fraction)
        data_to_use = data0[-laststeps:]
    else:
        data_to_use = data_equil
    mean, err, bn, skip = calc_block_avg(data_to_use, step_size, block_size, time_to_skip = time_to_skip)
    # Calculate Compressibility Modulus using largest possible data set
    N_steps= np.shape(data_to_use)[0]     # number of data points
    print("Using last {} ns data for compressibility calculation ...".format(round(N_steps*step_size,2)))
    mean2, err2, bn2, skip2 = calc_block_avg(data_to_use, step_size, block_size, time_to_skip = 0, verbose = False)
    VarS=0      # VarS is variance after the system reaches equilibrium
    for i in range(0, N_steps, 1):
        VarS=VarS+(data_to_use[i] - mean2)**2
    Variance=VarS/N_steps
    Ka = kb * temperature * mean / (nlipid * Variance * 1000) # order is corrected/recovered using 1000
    print('Final SA_avg/lip and err: %3.2f +- %3.2f' % (mean, err))
    print('Compressibility modulus is: {0:.2f}'.format(Ka))
    # plot module is in fflip.plot
    if plot:
        file_to_save = 'sa.png'
        y_max = area.max() + (area.max() - area.min())*0.3
        y_min = area.min() - (area.max() - area.min())*0.3
        plot_area(file_to_save, file_to_read, skip_ns = 0, interval = 0.002, unit_of_area_data = 'nm', area_range = (y_min, y_max)):
