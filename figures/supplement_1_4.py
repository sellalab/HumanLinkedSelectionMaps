__author__ = 'davidmurphy'


import os, re
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import simulations.cluster_sim_config as scfg
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import adjust_arrays, load_saved
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir
from figures.common_functions import format_panels


final_dir = root_dir + '/result/final_files'

#%%

# fdir = root_dir + '/result/final_files/ds_00_bound_bsx0.4_n2_sampled_mar2019/'
fldr1 = '/opttests/bs1_detrm/'
fldr2 ='/bound_bsx_n2_sampled_feb2019/'
#### fdir = final_dir + '/ds00rnd00bth0.65/'
fldr3 = '/opttests/bs1cs1/'

ddict_list = []
for fldr in [fldr1, fldr2, fldr3]:
    fdir = final_dir + fldr
    dr = defaultdict(list)
    for f in os.listdir(fdir):
        if f.endswith('.txt'):
            fname = fdir + f
            dist = re.search('rnd_\d', fname).group()
            dr[dist].append(ChromStruct(chrom='chr1', init=fname))
    ddict_list.append(dr)

#%% GET SIMULATION RUN PARAMS
def get_params(rlst):
    # get an indexed list of all likelihoods for run list
    ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
    # get indices of the top 3 best LH runs (after sorting on LH)
    top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

    pi0 = rlst[top3[0]].params[-1]
    tpi0 = rlst[top3[0]].stat.true_params[-1]

    # get the true params for the run and plot as the first variable
    tbprm = 10 ** rlst[0].stat.true_params[:6]
    tbpmf = tbprm * 1e8

    # plot results from the top3 indices of runs (averaged)
    bpmfs = []
    for i in top3:
        r = rlst[i]
        bpmfs.append(r.uvec[0])
    bpmf = np.average(np.array(bpmfs), axis=0) * 1e8

    clh, tclh = rlst[top3[0]].stat.best_lh, rlst[top3[0]].stat.true_clh

    if rlst[0].cnum:
        cpmfs = []
        for i in top3:
            r = rlst[i]
            cpmfs.append(r.avec[0])
        cpmf = np.average(np.array(cpmfs), axis=0)
        tcpmf = 10 ** rlst[0].stat.true_params[6:-1]

        return tbpmf, bpmf, tcpmf, cpmf, tpi0, pi0, tclh, clh

    else:

        return tbpmf, bpmf, tpi0, pi0, tclh, clh


#%% FOUR PANELS OF BS SIMULATION RESULTS
def four_panel_bs_simulation_results(keys, dr, rnd=False):
    """plot 4 sets of deterministic sim results"""
    # figure quadrants
    # quads = (221, 222, 223, 224)
    q1 = [(1,7), (9, 15), (17, 23), (25, 31)]
    q2 = [8, 16, 24, 32]
    lbls = ['a', 'b', 'c', 'd']
    x1, x2 = 0.0925, 0.55
    x_pos = [x1, x2, x1, x2]
    y1, y2 = 0.94, 0.49
    y_pos = [y1, y1, y2, y2]
    infcol = 'slategray'
    fig = plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(left=0.085, right=0.995, top=0.995, bottom=0.125,
                        hspace=0.075, wspace=0.15)
    # plot each quadrant of results
    for (i, k) in enumerate(keys):
        # modify panel format
        ax1 = plt.subplot(2, 16, q1[i])
        format_panels(ax1)

        # plot true and inferred params
        # tpmf, tpi0, tclh, mpmf, pi0, clh = get_params(dr[k])
        tpmf, mpmf, tpi0, pi0, tclh, clh = get_params(dr[k])

        # fix params to sum to 1
        tpmf /= sum(tpmf)
        mpmf /= sum(mpmf)
        xi = np.arange(6)
        if (i < 2) and (not rnd):
            mpmf = tpmf
        plt.bar(xi-0.2, tpmf, 0.4, color='k', label='true')
        plt.bar(xi+0.2, mpmf, 0.4, color=infcol, label='inferred')

        # conditional axis labels, legend
        ytck = np.arange(0, 0.46, 0.1)
        if i in [0, 2]:
            plt.yticks(ytck)
            # plt.ylabel(r'$\mu_{del}\ (\times 10^{-8})$')
        else:
            plt.yticks(ytck, ['']*len(ytck))
        plt.ylim(0, 0.47)
        if i > 1:
            x_exp = np.linspace(-4.5, -2.0, 6)
            plt.xticks(xi, [r'$10^{%.1f}$' % x for x in x_exp], fontsize=9)
            plt.xlabel('deleterious fitness effect')
        else:
            plt.xticks(xi, ['']*6)

        if i == 0:
            plt.legend(loc='upper center')

        # add letters to plots
        plt.text(x_pos[i], y_pos[i], lbls[i], transform=plt.gcf().transFigure,
                 fontweight='bold')

        # plot pi0
        ax2 = plt.subplot(2, 16, q2[i])
        format_panels(ax2)
        plt.bar(0, tpi0, 1,  color='k')
        plt.bar(1, pi0, 1, color=infcol)
        plt.xticks([])
        plt.yticks([])
        if i > 1:
            plt.xlabel(r'$\pi_0$')

        # # plot clh-true_clh
        # ax3 = plt.subplot(2, 16, q2[i]+1)
        # format_panels(ax3)
        # diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
        # col = 'green' if diff > 0 else 'red'
        # plt.bar(0, diff, 0.8, color=col)
        # plt.ylim(0.99 * diff, 1.01 * diff)
        # plt.xticks([])
    fig.add_subplot(111, frame_on=False)
    plt.yticks([])
    plt.xticks([])
    plt.ylabel(r'mutation rate ($u_d \times 10^{-8})$', labelpad=25)
    if rnd:
        f_save = final_dir + '/sfigs/fig_S5.optimization_random_n2.png'
    else:
        f_save = final_dir + '/sfigs/fig_S4.optimization_deterministic.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


try_keys = ['rnd_{}'.format(x) for x in [0, 3, 7, 1]]
four_panel_bs_simulation_results(try_keys, dr=ddict_list[1], rnd=True)
four_panel_bs_simulation_results(try_keys, dr=ddict_list[0])


#%% TWO PANELS OF BS+CS SIMULATION RESULTS
def two_panel_bscs_simulation_results(keys, ddict):
    """plot 4 sets of deterministic sim results"""
    # figure quadrants
    # quads = (221, 222, 223, 224)
    q1 = [(1, 7), (17, 23), (9, 15), (25, 31)]
    # q2 = [8, 16]
    q2 = [16, 32]
    lbls = ['a', 'b','(ii)', '(ii)']
    # pi0lbls =
    x1, x2 = 0.075, 0.545
    x_pos = [x1, x1, x2,  x2]
    y1, y2 = 0.95, 0.49
    y_pos = [y1, y2, y1, y2]
    infcol = 'slategray'
    fig = plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(left=0.07, right=0.995, top=0.995, bottom=0.115,
                        hspace=0.1, wspace=0.4)

    # plot each quadrant of results
    for (i, k) in enumerate(keys):
        # get true and inferred params for BS and CS
        tbpmf, bpmf, tcpmf, cpmf, tpi0, pi0, tclh, clh = get_params(ddict[k])

        # fix params to sum to 1
        tot_b = sum(tbpmf)
        tot_c = sum(tcpmf)
        tbpmf /= tot_b
        tcpmf /= (4.0*tot_c)
        bpmf /= tot_b
        cpmf /= (4.0*tot_c)

        # plot BS params
        ax1 = plt.subplot(2, 16, q1[i])
        format_panels(ax1)
        xi = np.arange(6)
        plt.bar(xi - 0.2, tbpmf, 0.4, color='k', label='true')
        plt.bar(xi + 0.2, bpmf, 0.4, color=infcol, label='inferred')

        # conditional axis labels, legend
        ytck = np.arange(0, 0.46, 0.1)
        plt.yticks(ytck, x=0.03)
        # plt.ylabel(r'$u_d\ (\times 10^{-8})$', labelpad=2)
        if i == 0:
            plt.legend(loc='upper right')
        # if i == 0:
        #     plt.yticks(ytck, x=0.03)
        #     plt.ylabel(r'$u_d\ (\times 10^{-8})$', labelpad=2)
        #     plt.legend(loc='upper right')
        # else:
        #     plt.yticks(ytck, [''] * len(ytck))
        plt.ylim(0, 0.47)
        if i == 1:
            x_exp = np.linspace(-4.5, -2.0, 6)
            plt.xticks(xi, [r'$10^{%.1f}$' % x for x in x_exp], fontsize=9, y=0.03)
            plt.xlabel('deleterious fitness effect', labelpad=2)
        else:
            plt.xticks(xi, ['' for _ in xi], fontsize=9, y=0.03)

        # add letters to plots
        plt.text(x_pos[i], y_pos[i], lbls[i], transform=plt.gcf().transFigure,
                 fontweight='bold')

        # plot CS params
        ax2 = plt.subplot(2, 16, q1[i+2])
        format_panels(ax2)
        xi = np.arange(6)
        plt.bar(xi - 0.2, tcpmf, 0.4, color='k', label='true')
        plt.bar(xi + 0.2, cpmf, 0.4, color=infcol, label='inferred')

        # conditional axis labels, legend
        ytck = np.arange(0, 0.15, 0.03)
        plt.yticks(ytck, x=0.03)
        # if i == 0:
        #     plt.yticks(ytck, x=0.03)
        #     plt.ylabel(r'$\alpha$', labelpad=2)
        # else:
        #     plt.yticks(ytck, [''] * len(ytck))
        plt.ylim(0, 0.14)
        if i == 1:
            x_exp = np.linspace(-4.5, -2.0, 6)
            plt.xticks(xi, [r'$10^{%.1f}$' % x for x in x_exp], fontsize=9, y=0.03)
            plt.xlabel('adaptive fitness effect', labelpad=2)
        else:
            plt.xticks(xi, ['' for _ in xi], fontsize=9, y=0.03)
        # # add letters to plots
        # plt.text(x_pos[i+2], y_pos[i+2], lbls[i+2],
        #          transform=plt.gcf().transFigure)

        # plot pi0
        # ax3 = plt.subplot(1, 16, q2[i])
        ax3 = plt.subplot(2, 16, q2[i])
        format_panels(ax3)
        plt.bar(0, tpi0, 1, color='k')
        plt.bar(1, pi0, 1, color=infcol)
        plt.xticks([])
        plt.yticks([])
        if i == 1:
            plt.xlabel(r'$\pi_0$')

    # udel y label
    fig.add_subplot(111, frame_on=False)
    plt.yticks([])
    plt.xticks([])
    plt.ylabel(r'mutation rate ($u_d \times 10^{-8})$', labelpad=18)
    # alpha y label
    fig.add_subplot(122, frame_on=False)
    plt.yticks([])
    plt.xticks([])
    plt.ylabel(r'proportion beneficial substitutions ($\alpha$)', labelpad=55)
    # if rnd:
    #     f_save = final_dir + '/sfigs/fig_S2.optimization_random_n2.png'
    # else:
    f_save = final_dir + '/sfigs/fig_S6.optimization_random_BSCS.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


try_keys = ['rnd_{}'.format(x) for x in [2, 1]]
two_panel_bscs_simulation_results(try_keys, ddict=ddict_list[2])


#%% CREATE RANDOM INITIAL CONDITIONS PLOT
def two_step_minimization_figure():
    """create figures of random initial conditions"""
    # PLOT 1: ALL 15 CONDITIONS, PICK 3 RANDOMLY AS "WINNERS" AND SAVE VALS
    tot = [0.05, 0.5, 5]
    picks = [4, 10, 14]
    ttls = [r'$u_d=$' + str(t) + r'$\times 10^{-8}$' for t in tot]
    ylbls = ['set 1: nearly neutral', 'set 2: moderate selection',
             'set 3: strong selection']
    xi = np.arange(6)
    winners = []
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.9, bottom=0.01,
                        wspace=0.2)
    for i in xrange(5):
        for j in xrange(3):
            q = 1 + i*3 + j
            ax = plt.subplot(5, 3, q)
            format_panels(ax)
            params = np.random.dirichlet([1]*6) * tot[j]
            if q in picks:
                col = 'fuchsia'
                winners.append(params)
            else:
                col = 'k'
            if i == 0:
                plt.title(ylbls[j], y=0.91)
            plt.bar(xi, params, color=col)
            plt.xticks([])
            plt.yticks([])
            if i == 2:
                plt.ylabel(ttls[j], labelpad=2)

    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_1.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

    # PLOT 2: 3 "WINNERS" FROM PLOT 1
    xi = np.arange(6)
    plt.figure(figsize=(2.16, 0.6 * 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.87, bottom=0.01,
                        wspace=0.2)

    for j in xrange(3):
        ax = plt.subplot(3, 1, j + 1)
        format_panels(ax)
        params = winners[j]
        if j == 0:
            plt.title('top 3 results', y=0.85)
        plt.bar(xi, params, color='fuchsia')
        plt.xticks([])
        plt.yticks([])

    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_2.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

    # PLOT 3: MEAN OF "WINNERS" FROM PLOT 2
    mparams = np.average(np.array(winners), axis=0)
    plt.figure(figsize=(2.16, 0.2 * 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.62, bottom=0.01,
                        wspace=0.2)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title('mean of top 3 results', y=0.84)
    plt.bar(xi, mparams, color='fuchsia')
    plt.xticks([])
    plt.yticks([])

    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_3.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


two_step_minimization_figure()


#%% CREATE RANDOM INITIAL CONDITIONS PLOT (updated for paper june 2021)
def two_step_minimization_figure_update2021():
    """create figures of random initial conditions"""
    # PLOT 1: ALL 15 CONDITIONS, PICK 3 RANDOMLY AS "WINNERS" AND SAVE VALS
    tot = [0.05, 0.5, 5]
    picks = [4, 10, 14]
    ttls = [r'$u_d=$' + str(t) + r'$\times 10^{-8}$' for t in tot]
    ylbls = ['set 1: nearly neutral', 'set 2: moderate selection',
             'set 3: strong selection']
    xi = np.arange(6)
    # winners = []
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.9, bottom=0.01,
                        wspace=0.2)

    # initial params (input 1)
    for i in xrange(5):
        for j in xrange(3):
            q = 1 + i*3 + j
            ax = plt.subplot(5, 3, q)
            format_panels(ax)
            params = np.random.dirichlet([1]*6) * tot[j]
            if q in picks:
                col = 'k'
                # winners.append(params)
            else:
                col = 'k'
            if i == 0:
                plt.title(ylbls[j], y=0.91)
            plt.bar(xi, params, color=col)
            plt.xticks([])
            plt.yticks([])
            if i == 2:
                plt.ylabel(ttls[j], labelpad=2)
    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_1.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

    # PLOT 2: 3 "WINNERS" FROM PLOT 1 (output 1)
    xi = np.arange(6)
    plt.figure(figsize=(2.16, 0.6 * 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.85, bottom=0.01,
                        wspace=0.2)

    input2 = []
    for j in xrange(3):
        ax = plt.subplot(3, 1, j + 1)
        format_panels(ax)
        # pick a NEW set for the output
        params = np.random.dirichlet([1]*6)
        input2.append(params)
        if j == 0:
            plt.title('top 3 results (output 1)', y=0.87)
        plt.bar(xi, params, color='fuchsia')
        plt.xticks([])
        plt.yticks([])
    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_2.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

    # PLOT 3: MEAN OF "WINNERS" FROM PLOT 2 (input 2)
    mparams = np.average(np.array(input2), axis=0)
    plt.figure(figsize=(2.16, 0.4 * 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.55, bottom=0.01,
                        wspace=0.2)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title('mean of top 3 results\n(input 2)', y=0.91)
    plt.bar(xi, mparams, color='k')
    plt.xticks([])
    plt.yticks([])
    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_3.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

    # PLOT 4: MEAN OF "WINNERS" FROM PLOT 2
    mparams = np.random.dirichlet([1]*6)
    plt.figure(figsize=(2.16, 0.4 * 2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.55, bottom=0.01,
                        wspace=0.2)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title('maximum likelihood\nparameters (output 2)', y=0.91)
    plt.bar(xi, mparams, color='fuchsia')
    plt.xticks([])
    plt.yticks([])

    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_4.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


two_step_minimization_figure_update2021()

#%% CREATE RANDOM INITIAL CONDITIONS PLOT
def step_two_initial():
    """create figures of random initial conditions"""
    xi = np.arange(6)
    plt.figure(figsize=(2.16, 0.6*2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.9, bottom=0.01,
                        wspace=0.2)
    all_params = []
    for j in xrange(3):
        ax = plt.subplot(3, 1, j+1)
        format_panels(ax)
        params = np.random.dirichlet([1]*6)
        all_params.append(params)
        if j == 0:
            plt.title('top 3 results', y=0.85)
        plt.bar(xi, params, color='fuchsia')
        plt.xticks([])
        plt.yticks([])

    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_2.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

    # plot mean params
    mparams = np.average(np.array(all_params), axis=0)
    plt.figure(figsize=(2.16, 0.2*2.16))
    plt.subplots_adjust(left=0.035, right=0.995, top=0.7, bottom=0.01,
                        wspace=0.2)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title('mean of top 3 results', y=0.84)
    plt.bar(xi, mparams, color='fuchsia')
    plt.xticks([])
    plt.yticks([])

    f_save = final_dir + '/sfigs/fig_S1_twostep_minimization_3.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

step_two_initial()
#%%
