# -*- coding: utf-8 -*-


def plot_area(file_to_save, data, skip_ns=0, interval=0.002, area_range=(35, 80), title=''):
    import numpy as np
    from matplotlib import pyplot as plt
    skip_frames = int(skip_ns/interval)
    data_skipped = data[skip_frames:]
    cumulative = np.cumsum(data_skipped)/np.arange(1, data_skipped.shape[0]+1, 1)
    pattern = ['-','-']
    linewidth = [3, 3]
    markersize = [3, 3]
    lg = ['instantaneous','cumulative']
    colors = ['deepskyblue', 'deeppink']
    csfont = {'fontname':'Comic Sans MS'}
    font = {'fontname': 'serif', 'weight': 500}
    fig, ax = plt.subplots(figsize=(4, 3))
    
    x = np.arange(0, data.shape[0] * interval, interval)
    plt.xlim((0, data.shape[0] * interval))
    plt.ylim(area_range)
    plt.xlabel("Time [ns]", fontsize=14, **font,) #verticalalignment="center")
    plt.ylabel("Area / Lipid [$\mathrm{\AA}^2$]", fontsize=14, **font,)# rotation=90,)
    # plt.ylabel("Area / Lipid [" + "$\mathrm{\AA}^2$]" + "$^2$]", fontsize=14, **font,)# rotation=90,)

    plt.xticks(np.arange(0, data.shape[0] * interval, 
                         int(((data.shape[0] + 1) * interval)/10) + 5), **font
    )
    plt.yticks(np.arange(area_range[0], area_range[1], 4), **font)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2) #set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    
    x = np.arange(0, data.shape[0] * interval, interval)

    for i, y in enumerate([data, cumulative]):
        if i==0:
            ax.plot(x, y, pattern[i], linewidth=linewidth[i], markersize=markersize[i], 
                    label=lg[i], color=colors[i])
        else:
            ax.plot(x[skip_frames:], y, pattern[i], linewidth=linewidth[i], markersize=markersize[i],
                    label = lg[i], color = colors[i])
        ax.legend(loc='lower left', fontsize=8, handlelength=2, edgecolor='k')
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
    ax.text(.5,.9, title, fontsize=14, horizontalalignment='center', 
            transform=ax.transAxes)
    plt.savefig(file_to_save, dpi=300, bbox_inches = 'tight')
