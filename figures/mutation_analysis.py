__author__ = 'davidmurphy'

import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, \
    cst_from_fldr
from figures.other_code.sub_rate_figs import neutcons_uratio_single
from figures.common_functions import format_panels

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'
figdir = final_dir + '/mainfigs'


#%%
# pct = [95, 94, 93, 92]
# s_ape = [1 - neutcons_uratio_single('ape', pc=p)[0].mean() for p in pct]
# s_cad = [1 - neutcons_uratio_single('cadd', pc=p)[0].mean() for p in pct]
# s_fis = [1 - neutcons_uratio_single('fish', pc=p)[0].mean() for p in pct]


#%% FIGURE 4 DATA
def figure_4_data():
    pct = [95, 94, 93]
    s_ape = [1 - neutcons_uratio_single('ape', pc=p)[0].mean() for p in pct]
    s_cad = [1 - neutcons_uratio_single('cadd', pc=p)[0].mean() for p in pct]
    s_fis = [1 - neutcons_uratio_single('fish', pc=p)[0].mean() for p in pct]

    # get inferred udel for ape and CADD
    u_ape = []
    for p in pct:
        if p != 94:
            fldr = 'ape_cons{}_gmask'.format(p)
        else:
            fldr = 'ape_cons94_gmask_mnb_378'
        # fldr = 'ape_cons{}_gmask'.format(p)
        cst = cst_from_fldr(fldr)
        udel = cst.stat.utot[0] * 1e8
        u_ape.append(udel)
    u_cad = []
    for p in pct:
        fldr = 'cadd{}_gmask'.format(p)
        cst = cst_from_fldr(fldr)
        udel = cst.stat.utot[0] * 1e8
        u_cad.append(udel)
    u_fis = []
    for p in pct:
        if p != 94:
            fldr = 'fish_cons{}_gmask'.format(p)
        else:
            fldr = 'fish_cons94_gmask_mnb_378'
        cst = cst_from_fldr(fldr)
        udel = cst.stat.utot[0] * 1e8
        u_fis.append(udel)

    # return s_ape, s_cad, u_ape, u_cad
    return s_ape, s_cad, s_fis, u_ape, u_cad, u_fis


s_ape, s_cad, s_fis, u_ape, u_cad, u_fis = figure_4_data()


#%% FIGURE 4
def figure_4(s_ape, s_cad, s_fis, u_ape, u_cad, u_fis):
    #  upper and lower bounds of u total
    ulo = 1.29
    uhi = 1.51
    pct = [95, 94, 93]

    # plt.figure(figsize=(5, 3))
    fig = plt.figure(figsize=(4.875, 1.95))
    plt.subplots_adjust(wspace=0.05, left=0.13, right=0.995, top=0.995, bottom=0.21)
    # panel A: ape rates
    ax1 = plt.subplot(131)
    format_panels(ax1)

    x1 = np.arange(len(pct))
    utop = np.array([u / ulo for u in u_ape])
    ubot = np.array([u / uhi for u in u_ape])
    uheight = utop-ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_ape, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('ape', labelpad=2)
    plt.ylabel('proportion of\ndeleterious mutations', labelpad=2)
    plt.yticks(x=0.06)
    plt.ylim(0, 1.15)

    # panel B: fish rates
    ax2 = plt.subplot(132)
    format_panels(ax2)
    utop = np.array([u / ulo for u in u_fis])
    ubot = np.array([u / uhi for u in u_fis])
    uheight = utop-ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_fis, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('vertebrate', labelpad=2)
    plt.yticks(color='none', x=0.05)
    plt.ylim(0, 1.15)

    # panel C: CADD rates
    ax3 = plt.subplot(133)
    format_panels(ax3)
    utop = np.array([u / ulo for u in u_cad])
    ubot = np.array([u / uhi for u in u_cad])
    uheight = utop-ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_cad, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('CADD', labelpad=2)
    plt.yticks(color='none', x=0.05)
    plt.ylim(0, 1.15)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=6.5))
    # leg._legend_box.align = 'left'
    # plt.setp(leg.get_title(), fontsize='xx-small')
    # customize legend
    handles, labels = ax1.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc='lower center',
                     title='estimates from:', frameon=1, framealpha=1,
                     facecolor='white', prop=dict(size=8.5), ncol=2,
                     borderpad=0.2, labelspacing=0.2, handlelength=1.1,
                     handletextpad=0.4, bbox_to_anchor=(0.6, 0.195),
                     columnspacing=0.8)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=9),
    #                  ncol=2)
    leg._legend_box.align = 'left'
    leg.set_title('estimates from:', prop={'size':8.5})
    # plt.text(0.21, 0.92, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.62, 0.92, 'B', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/compare.urates.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


figure_4(s_ape, s_cad, s_fis, u_ape, u_cad, u_fis)
#%%

