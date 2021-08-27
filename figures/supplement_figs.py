__author__ = 'davidmurphy'

import os
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir
from figures.other_code.summary_slide import cst_from_fldr, rsq_from_fldr
from scipy.stats import linregress
from figures.common_functions import format_panels

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


#%% COMPARE SINGLE BS ANNOTATION RUNS


def single_bsanno(flist, llist, clist, sname, sfldr):
    # create lists for data to plot
    dfe_list = []
    pmf_list = []
    pi0_list = []
    clh_list = []
    rsq_list = []

    for fldr in flist:
        # get the rsq values
        r_log = root_dir + '/result/final_files/{}/rsq.log'.format(fldr)
        rsq = np.loadtxt(r_log)
        rsq_list.append(rsq)

        # initalize ChromStruct from final composite from folder
        cst = cst_from_fldr(fldr)

        # scale upmf by u0
        upmf = [u/u0 for u in cst.uvec]
        pmf_list.append(upmf)
        # get pi0 based on tau
        mpi0 = cst.params[-1] / cst.fixed.tau_init
        pi0_list.append(mpi0)
        # use representative DFE for all values
        dfe_list.append(cst.bdfe[0])
        # get best CLH
        clh_list.append(cst.stat.best_lh)

    # plot formatting parameters
    xi = np.arange(len(dfe_list[0]))
    xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]
    n = sum(len(pm) for pm in pmf_list)
    w = 0.8 / n
    s = -w * n / 2.0
    ncol = 2

    # create the plot
    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.06, wspace=1.75, right=1, bottom=0.15, top=0.92)
    axes_lbl_size = 9
    xtck_size = 7
    ytck_size = 8
    leg_size = 7.5
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize
    # plt.rc('figure', titlesize=10)  # fontsize of the figure title

    # DFE plot
    ax1 = plt.subplot(1, 8, (1, 3))
    format_panels(ax1)

    # plt.title('A (i)', loc='left', y=0.97)
    plt.title('(i)', loc='center', y=0.97)
    ymax = 0
    for (i, df) in enumerate(pmf_list):
        assert len(df) == 1
        lbl = llist[i]
        df = df[0]
        ymax = max(df.max(), ymax)
        plt.bar(xi + s, df, w, label=lbl, color=clist[i], align='edge')
        s += w

    plt.ylabel(r'mutation rate (units of $\mu_0$)', labelpad=3)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.04)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    plt.xlabel('deleterious fitness effect', labelpad=3)
    plt.legend(loc='upper left', ncol=ncol)

    # plot udel
    ax2 = plt.subplot(1,8,4)
    format_panels(ax2)

    plt.title('(ii)', loc='center', y=0.97)
    plt.ylabel(r'mutation rate (units of $\mu_0$)', labelpad=3)
    plt.yticks(x=0.04)
    s = -w / 2.0
    li = 0
    for (i, df) in enumerate(pmf_list):
        for d in df:
            udel = sum(d)
            plt.bar(0+s, udel, w, color=clist[li])
            li += 1
            s += w
    plt.xticks([])
    plt.yticks(x=0.2)

    # plot pi0
    ax3 = plt.subplot(1,8,5)
    format_panels(ax3)

    plt.title('(iii)', loc='center', y=0.97)
    plt.ylabel('average reduction\nin heterozygosity', labelpad=3)
    plt.yticks(x=0.2)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    # plot CLH
    ax4 = plt.subplot(1,8,6)
    format_panels(ax4)

    plt.title('(iv)', loc='center', y=0.97)
    plt.ylabel(r'$- \Delta CLLH\ (x10^5)$', labelpad=3)
    plt.yticks(x=0.2)
    s = -w * n / 2.0
    clh_list = [1e5 * (cl / min(clh_list) - 1) for cl in clh_list]
    for (i, cl) in enumerate(clh_list):
        print '{:<20} -CLLH {}'.format(flist[i], cl)
        plt.bar(0+s, cl, w, color=clist[i])
        s += w
    if max(clh_list) < 0.1:
        plt.ylim(0, 0.11)
    # plt.ylim(min(clh_list), max(clh_list))
    plt.xticks([])

    # R^2 PLOT
    ax5 = plt.subplot(1,8,(7,8))
    format_panels(ax5)

    # plt.title('B', loc='left', y=0.97)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=3)
    plt.yticks(x=0.06)
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    ws = ws[wmsk]
    for (i, rs) in enumerate(rsq_list):
        # save the 1Mb rsq vals
        rsq = rs[:,1][wmsk]
        plt.plot(ws, rsq, color=clist[i], marker='o',
                 linewidth=0, alpha=0.75, label=llist[i], ms=4)
    plt.xlabel('window size (bp)', labelpad=3)
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    # plt.legend(prop=dict(size=9), loc='upper left')

    # save the plot
    save_lbl = sname + '.parameters'
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)
    plt.text(0.01, 0.93, 'A', transform=plt.gcf().transFigure)
    plt.text(0.77, 0.93, 'B', transform=plt.gcf().transFigure)
    fsave = sdir + '/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=256)
    plt.close()


#%% COMPARE MULTI BS ANNOTATION RUNS


def multi_bsanno(flist, llist, subllist, clist, sname, sfldr):
    # create lists for data to plot
    dfe_list = []
    pmf_list = []
    pi0_list = []
    clh_list = []
    rsq_list = []

    for fldr in flist:
        # get the rsq values
        r_log = root_dir + '/result/final_files/{}/rsq.log'.format(fldr)
        rsq = np.loadtxt(r_log)
        rsq_list.append(rsq)

        # initalize ChromStruct from final composite from folder
        cst = cst_from_fldr(fldr)

        # scale upmf by u0
        upmf = [u/u0 for u in cst.uvec]
        pmf_list.append(upmf)
        # get pi0 based on tau
        mpi0 = cst.params[-1] / cst.fixed.tau_init
        pi0_list.append(mpi0)
        # use representative DFE for all values
        dfe_list.append(cst.bdfe[0])
        # get best CLH
        clh_list.append(cst.stat.best_lh)

    # plot formatting parameters
    xi = np.arange(len(dfe_list[0]))
    xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]
    n = sum(len(pm) for pm in pmf_list)
    w = 0.8 / n
    s = -w * n / 2.0
    ncol = 1
    hatches = ['.', '/', 'x']

    # create the plot
    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.06, wspace=1.75, right=1, bottom=0.15, top=0.92)
    axes_lbl_size = 9
    xtck_size = 7
    ytck_size = 8
    leg_size = 7.5
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize
    # plt.rc('figure', titlesize=10)  # fontsize of the figure title

    # DFE plot
    plt.subplot(1, 8, (1, 3))
    plt.title('(i)', loc='center', y=0.97)
    ymax = 0
    for (i, df) in enumerate(pmf_list):
        if len(df) == 1:
            lbl = llist[i]
            dd = df[0]
            ymax = max(dd.max(), ymax)
            plt.bar(xi + s, dd, w, label=lbl, color=clist[i], align='edge')
            s += w
        else:
            j = 0
            for dd in df:
                lbl = llist[i] + ': ' + subllist[j]
                ymax = max(dd.max(), ymax)
                if sum(dd) < 0.01:
                    continue
                plt.bar(xi + s, dd, w, label=lbl, color=clist[i],
                        hatch=hatches[j]*5, align='edge')
                s += w
                j += 1

    plt.ylabel(r'mutation rate (units of $\mu_0$)', labelpad=3)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.04)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    plt.xlabel('deleterious fitness effect', labelpad=3)
    plt.legend(loc='upper left', ncol=ncol)

    # plot udel
    plt.subplot(1,8,4)
    plt.title('(ii)', loc='center', y=0.97)
    plt.ylabel(r'mutation rate (units of $\mu_0$)', labelpad=3)
    plt.yticks(x=0.04)
    s = -w / 2.0
    for (i, df) in enumerate(pmf_list):
        if len(df) == 1:
            udel = sum(df[0])
            plt.bar(0+s, udel, w, color=clist[i])
            s += w
        else:
            j = 0
            for dd in df:
                udel = sum(dd)
                if udel < 0.01:
                    continue
                plt.bar(0 + s, udel, w, color=clist[i], hatch=hatches[j]*5)
                s += w
                j += 1

    plt.xticks([])
    plt.yticks(x=0.2)

    # plot pi0
    plt.subplot(1,8,5)
    plt.title('(iii)', loc='center', y=0.97)
    plt.ylabel('average reduction\nin heterozygosity', labelpad=3)
    plt.yticks(x=0.2)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    # plot CLH
    plt.subplot(1,8,6)
    plt.title('(iv)', loc='center', y=0.97)
    plt.ylabel(r'$- \Delta CLLH\ (x10^5)$', labelpad=3)
    plt.yticks(x=0.2)
    s = -w * n / 2.0
    clh_list = [1e5 * (cl / min(clh_list) - 1) for cl in clh_list]
    for (i, cl) in enumerate(clh_list):
        print '{:<20} -CLLH {}'.format(flist[i], cl)
        plt.bar(0+s, cl, w, color=clist[i])
        s += w
    if max(clh_list) < 0.1:
        plt.ylim(0, 0.11)
    # plt.ylim(min(clh_list), max(clh_list))
    plt.xticks([])

    # R^2 PLOT
    plt.subplot(1,8,(7,8))
    # plt.title('B', loc='left', y=0.97)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=3)
    plt.yticks(x=0.06)
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    ws = ws[wmsk]
    for (i, rs) in enumerate(rsq_list):
        # save the 1Mb rsq vals
        rsq = rs[:,1][wmsk]
        plt.plot(ws, rsq, color=clist[i], marker='o',
                 linewidth=0, alpha=0.75, label=llist[i], ms=4)
    plt.xlabel('window size (bp)', labelpad=3)
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    # plt.legend(prop=dict(size=9), loc='upper left')

    # save the plot
    save_lbl = sname + '.parameters'
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)
    plt.text(0.01, 0.93, 'A', transform=plt.gcf().transFigure)
    plt.text(0.77, 0.93, 'B', transform=plt.gcf().transFigure)
    fsave = sdir + '/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=256)
    plt.close()


#%% COMPARE BS TO CS RUN RESULTS


def bscs_joint(flist, llist, subllist, clist, sname, sfldr):
    # create lists for data to plot
    dfe_list = []
    pmf_list = []
    alf_list = []
    pi0_list = []
    clh_list = []
    rsq_list = []

    for fldr in flist:
        # get the rsq values
        r_log = root_dir + '/result/final_files/{}/rsq.log'.format(fldr)
        rsq = np.loadtxt(r_log)
        rsq_list.append(rsq)

        # initalize ChromStruct from final composite from folder
        cst = cst_from_fldr(fldr)
        cst.stat.calc_stats(cst)

        # scale upmf by u0
        upmf = [u/u0 for u in cst.uvec]
        if np.any(cst.uvec):
            pmf_list.append(upmf)
        # get alpha values
        if np.any(cst.avec):
            alf_list.append(cst.stat.avec)
        # get pi0 based on tau
        mpi0 = cst.params[-1] / cst.fixed.tau_init
        pi0_list.append(mpi0)
        # use representative DFE for all values
        if cst.bdfe:
            dfe_list.append(cst.bdfe[0])
        else:
            dfe_list.append(cst.cdfe[0])
        # get best CLH
        clh_list.append(cst.stat.best_lh)

    # plot formatting parameters
    xi = np.arange(len(dfe_list[0]))
    xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]
    n = sum(len(pm) for pm in pmf_list)
    w = 0.8 / n
    s = -w * n / 2.0
    ncol = 1
    hatches = ['.', '/', 'x']

    # create the plot
    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.06, wspace=1.75, right=1, bottom=0.18, top=0.92,
                        hspace=0.75)
    axes_lbl_size = 9
    xtck_size = 7
    ytck_size = 8
    leg_size = 7.5
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize
    # plt.rc('figure', titlesize=10)  # fontsize of the figure title

    # BS DFE plot
    plt.subplot(2, 8, (1, 3))
    plt.title('(i)', loc='center', y=0.94)
    ymax = 0
    for (i, df) in enumerate(pmf_list):
        if len(df) == 1:
            lbl = llist[i]
            dd = df[0]
            ymax = max(dd.max(), ymax)
            plt.bar(xi + s, dd, w, label=lbl, color=clist[i], align='edge')
            s += w
        else:
            j = 0
            for dd in df:
                lbl = llist[i]
                ymax = max(dd.max(), ymax)
                if sum(dd) < 0.01:
                    continue
                plt.bar(xi + s, dd, w, label=lbl, color=clist[i],
                        hatch=hatches[j]*5, align='edge')
                s += w
                j += 1

    plt.ylabel(r'$\mu_{del}$', labelpad=3)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.04)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    plt.xlabel('deleterious fitness effect', labelpad=3)
    plt.legend(loc='upper left', ncol=ncol)

    # CS DFE plot
    xi = np.arange(len(dfe_list[0]))
    n = sum(len(alf) for alf in alf_list)
    w = 0.8 / n
    s = -w * n / 2.0
    ncol = 2
    plt.subplot(2, 8, (9, 11))
    ymax = 0
    offset = len(pmf_list)
    j = 0
    for (i, df) in enumerate(alf_list):
        if len(df) == 1:
            lbl = llist[i+offset]
            dd = df[0]*100
            ymax = max(dd.max(), ymax)
            plt.bar(xi + s, dd, w, label=lbl, color=clist[i+offset],
                    align='edge')
            s += w
        else:
            if all(sum(dd) > 0.00001 for dd in df):
                use_hatch = True
            else:
                use_hatch = False
            for dd in df:
                dd *= 100
                lbl = llist[i+offset] + ': ' + sllist[j]
                ymax = max(dd.max(), ymax)
                if sum(dd) < 0.00001:
                    continue
                if use_hatch:
                    plt.bar(xi + s, dd, w, label=lbl, color=clist[i+offset],
                            align='edge', hatch=hatches[j-1]*5)
                    # j += 1

                else:
                    plt.bar(xi + s, dd, w, label=lbl, color=clist[i + offset],
                            align='edge')
                j += 1
                s += w

    plt.ylabel(r'$\alpha$ (%)', labelpad=3)
    plt.ylim(0, 1.5 * ymax)
    plt.yticks(x=0.04)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    plt.xlabel('adaptive fitness effect', labelpad=3)
    plt.legend(loc='upper left', ncol=ncol)

    # plot udel
    plt.subplot(2,8,4)
    plt.title('(ii)', loc='center', y=0.97)
    plt.ylabel(r'total $\mu_{del}$', labelpad=3)
    plt.yticks(x=0.04)
    s = -w / 2.0
    for (i, df) in enumerate(pmf_list):
        if len(df) == 1:
            udel = sum(df[0])
            plt.bar(0+s, udel, w, color=clist[i])
            s += w
        else:
            j = 0
            for dd in df:
                udel = sum(dd)
                if udel < 0.01:
                    continue
                plt.bar(0 + s, udel, w, color=clist[i], hatch=hatches[j]*5)
                s += w
                j += 1

    plt.xticks([])
    plt.yticks(x=0.2)

    # plot alpha
    plt.subplot(2,8,12)
    plt.ylabel(r'total $\alpha$ (%)', labelpad=3)
    plt.yticks(x=0.04)
    s = -w / 2.0
    j = 0
    for (i, df) in enumerate(alf_list):
        if len(df) == 1:
            alf = sum(df[0])*100
            plt.bar(0+s, alf, w, color=clist[i+offset])
            s += w
        else:

            if all(sum(dd) > 0.00001 for dd in df):
                use_hatch = True
            else:
                use_hatch = False
            for dd in df:
                alf = sum(dd)
                if alf < 0.000001:
                    continue
                if use_hatch:
                    plt.bar(0 + s, alf, w, color=clist[i + offset],
                            hatch=hatches[j]*5)
                    j += 1

                else:
                    plt.bar(0 + s, alf, w, color=clist[i+offset])
                s += w

    plt.xticks([])
    plt.yticks(x=0.2)

    # plot pi0
    plt.subplot(1,8,5)
    plt.title('(iii)', loc='center', y=0.97)
    plt.ylabel('average reduction\nin heterozygosity', labelpad=3)
    plt.yticks(x=0.2)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    # plot CLH
    plt.subplot(1,8,6)
    plt.title('(iv)', loc='center', y=0.97)
    plt.ylabel(r'$- \Delta CLLH\ (x10^5)$', labelpad=3)
    plt.yticks(x=0.2)
    s = -w * n / 2.0
    clh_list = [1e5 * (cl / min(clh_list) - 1) for cl in clh_list]
    for (i, cl) in enumerate(clh_list):
        print '{:<20} -CLLH {}'.format(flist[i], cl)
        plt.bar(0+s, cl, w, color=clist[i])
        s += w
    if max(clh_list) < 0.1:
        plt.ylim(0, 0.11)
    # plt.ylim(min(clh_list), max(clh_list))
    plt.xticks([])

    # R^2 PLOT
    plt.subplot(1,8,(7,8))
    # plt.title('B', loc='left', y=0.97)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=3)
    plt.yticks(x=0.06)
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    ws = ws[wmsk]
    for (i, rs) in enumerate(rsq_list):
        # save the 1Mb rsq vals
        rsq = rs[:,1][wmsk]
        plt.plot(ws, rsq, color=clist[i], marker='o',
                 linewidth=0, alpha=0.75, label=llist[i], ms=4)
    plt.xlabel('window size (bp)', labelpad=3)
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    # plt.legend(prop=dict(size=9), loc='upper left')

    # save the plot
    save_lbl = sname + '.parameters'
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)
    plt.text(0.01, 0.93, 'A', transform=plt.gcf().transFigure)
    plt.text(0.77, 0.93, 'B', transform=plt.gcf().transFigure)
    fsave = sdir + '/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=256)
    plt.close()


#%% COMPARE R^2 AT WINDOW SUBSET


def rsq_3_size(flist, llist, sname, sfldr, xlab, rotation=0, fletter='C'):
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)

    # create string format for rsq files
    fstr = root_dir + '/result/final_files/{}/rsq.log'

    # keep just 3 window sizes: 62.5K, 250K, 1M
    keep_wins = [62500, 2.5e5, 1e6]
    wlab = ['62.5Kb', '250Kb', '1Mb']
    rsq_list = []
    for fl in flist:
        # load rsq and window size
        win, rsq = np.loadtxt(fstr.format(fl)).T
        # keep rsq at the subset of windows needed
        idx = np.in1d(win, keep_wins)
        rsq_list.append(rsq[idx])

    # convert rsq values to array of columns for each size
    r = np.array(rsq_list)

    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.1, bottom=0.25, wspace=0.25, right=0.985,
                        top=0.9)
    axes_lbl_size = 9
    xtck_size = 8
    ytck_size = 8
    leg_size = 7.5
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize

    # plot subplot for each window size
    for pidx in xrange(1, 4):
        ax = plt.subplot(1, 3, pidx)
        format_panels(ax)
        # get column for current window size
        ri = r[:,pidx-1]
        print ri
        # rmin = int(1000*min(ri)) * 1e-3
        # rmax = int(1000*max(ri)+1) * 1e-3
        rmin = round(min(ri)-0.0005, 3)
        rmax = round(max(ri)+0.0005, 3)
        # rmin, rmax = round(min(ri), 3), round(max(ri), 3)
        r_interval = (rmax - rmin) / 5.0
        ytck = np.arange(rmin, rmax+r_interval, r_interval)
        # ytck = [round(rv, 3) for rv in ]

        xi = np.arange(len(ri))
        plt.title('window = ' + wlab[pidx-1], y=0.96)
        plt.plot(xi, ri, lw=0, marker='o', ms=6, color='k',
                 alpha=0.8)
        plt.xlabel(xlab)
        plt.xticks(xi, [l[:9] for l in llist], rotation=rotation, y=0.05,
                   ha='center')
        plt.yticks(ytck, x=0.05)
        plt.ylim(rmin - 0.5*r_interval, rmax + 0.5*r_interval)
        if pidx == 1:
            ylab = r'variance explained ($R^2$)'
            plt.ylabel(ylab, labelpad=3)

    # figure letter
    plt.text(0.01, 0.93, fletter, transform=plt.gcf().transFigure)

    f_save = sdir + '/{}.3range.rsq.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% COMRPARE MULTIPLE CHR1 PREDICTIONS


def combined_chr1(flist, llist, clist, sname, sfldr, fletter='C'):
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)

    obs_flag = False
    plt.figure(figsize=(6.5, 1.75))
    plt.subplots_adjust(left=0.07, bottom=0.2, right=1, top=0.9)
    axes_lbl_size = 9
    xtck_size = 8
    ytck_size = 8
    leg_size = 8
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize
    xi = None
    ax = plt.subplot(111)
    format_panels(ax)
    for (i, rdir) in enumerate(flist):
        fdir = root_dir + '/result/final_files/{}/'.format(rdir)
        chlist = [f for f in os.listdir(fdir) if ('chr1' in f) and
                 (f.endswith('.txt')) and 'TMRCA' not in f]
        f_name = fdir + chlist[0]
        prd, obs, num = np.loadtxt(f_name).T
        xi = np.arange(0, prd.size / 2.0, 0.5)
        pi_mean = np.nanmean(obs)
        prd /= pi_mean
        obs /= pi_mean
        obs[(obs < 0.25) | (obs > 1.75)] = np.nan

        if not obs_flag:
            plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
            obs_flag = True

        plt.plot(xi, prd, color=clist[i], lw=1.1, alpha=0.8, label=llist[i])

    plt.ylabel(r'heterozygosity ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(-0.4, 1.9)
    plt.ylim(0, 1.9)
    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    plt.xticks(xrange(25, 250, 50), y=0.05)
    plt.xlim(0, xi[-1])
    plt.legend(loc='lower center', ncol=len(llist)+1, frameon=1,
               framealpha=0.75, facecolor='white', handlelength=1.2)

    # figure letter
    plt.text(0.01, 0.93, fletter, transform=plt.gcf().transFigure)
    f_save = sdir + '/{}.chr1.combined.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% COLLATED PLOT AND MULTIPLE SORTED PREDICTIONS COMBINED


def collate_and_sort(flist, llist, clist, sname, sfldr, fletters=('D', 'E')):
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)

    obs_flag = False
    # plt.figure(figsize=(6.5, 2))
    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.08, bottom=0.2, right=0.99, top=0.98, wspace=0.25)
    axes_lbl_size = 9
    xtck_size = 8
    ytck_size = 8
    leg_size = 8
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize

    # COLLATED PLOT
    # plt.subplot(121)
    ax1 = plt.subplot(1, 3, (1,2))
    format_panels(ax1)
    for (i, rdir) in enumerate(flist):
        f_ape = root_dir + '/result/final_files/{}/nonsyn.collate.5.00e-03width.npy'.format(
            rdir)
        f_mcv = root_dir + '/result/collate/YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy'
        bins, div, pi, pr, cnts = np.load(f_ape).T
        _, _, _, mcv, mcv_cnts = np.load(f_mcv).T
        # bins, pi, div, pred, cnts = np.load(f_in).T
        obs = pi / cnts
        prd = pr / cnts
        mcv /= mcv_cnts
        # y_max = pi0
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        mcv /= y_max
        # plt.title(rdir)
        # pc = '{}%'.format(100-int(rdir.split('_')[1][4:]))
        if not obs_flag:
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.1)
            obs_flag = True
        # plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.1)
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.6, lw=0.8)

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'heterozygosity ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    plt.legend(prop={'size':7})
    # plt.ylim(0.83, 1.13)

    # PREDICTED VS. OBSERVED PLOT
    # plt.subplot(122)
    ax2 = plt.subplot(133)
    format_panels(ax2)
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    if 'McVicker' in llist:
        plt.text(0.5, 1.05, 'without linked selection', ha='center',
                 va='center')
    else:
        plt.text(0.75, 1.05, 'without linked selection', ha='center',
                 va='center')
    ymin = 1
    for (i, rdir) in enumerate(flist):
        # load data
        fdir = root_dir + '/result/final_files/{}/'.format(rdir)
        slist = [f for f in os.listdir(fdir) if ('predsort' in f) and
                 (f.endswith('.txt'))]
        if len(slist) == 1:
            f_name = fdir + slist[0]
        else:
            f_name = fdir + 'basic_sort_n100.txt'
        rst = cst_from_fldr(rdir)
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        f_in = f_name

        # normalize by pi0
        print f_in
        div, pi, pred = np.loadtxt(f_in).T
        pi /= pi0
        pred /= pi0
        ymin = min(ymin, pred.min())
        plt.plot(pred, pi, color='darkslategray', marker='o', ms=5,
                 markerfacecolor='None', markeredgecolor=clist[i],
                 markeredgewidth=0.9, lw=0, alpha=0.75)

    if 'McVicker' in llist:
        plt.ylim(0., 1.15)
        plt.xlim(0., 1.02)
        xtick = np.arange(0.1, 1.01, 0.1)
        ytick = np.arange(0.1, 1.11, 0.1)
    else:
        axmin = 0.5
        # plt.ylim(0.48, 1.15)
        # plt.xlim(0.48, 1.02)
        # xtick = np.arange(0.5, 1.01, 0.1)
        # ytick = np.arange(0.5, 1.11, 0.1)
        plt.ylim(axmin, 1.15)
        plt.xlim(axmin, 1.02)
        xtick = np.arange(axmin, 1.01, 0.1)
        ytick = np.arange(axmin, 1.11, 0.1)

    plt.xticks(xtick, y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick, x=0.02)
    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.legend(loc='lower right')

    plt.text(0.01, 0.93, fletters[0], transform=plt.gcf().transFigure)
    plt.text(0.68, 0.93, fletters[1], transform=plt.gcf().transFigure)

    f_save = sdir + '/{}.collated.sort.combined.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% R^2 MULTIPLE RUNS


def compare_rsq(flist, llist, clist, sname, sfldr, fletter='B'):
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(flist[0]))[:16,0])
    print w
    for fl in flist:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    plt.figure(figsize=(2.5,2.5))
    plt.subplots_adjust(right=1, top=1, left=0.25, bottom=0.17)
    for i in xrange(len(flist)):
        rs = rsq_list[i]
        lab = llist[i]
        col = clist[i]
        plt.plot(w, rs, label=lab, color=col, marker='o', lw=0, ms=6,
                 alpha=0.8)

    plt.xlabel('window size (log-scale)')
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4, 5, 6]])
    plt.ylim(0, 0.65)
    plt.ylabel('proportion of variance\n' + r'explained $(R^2)$')
    plt.legend(prop=dict(size=9), loc='upper left', handletextpad=0.1,
               borderaxespad=0.01)
    plt.text(0.01, 0.95, fletter, transform=plt.gcf().transFigure)

    f_save = final_dir + '/{}/{}.rsq.png'.format(sfldr, sname)
    plt.savefig(f_save, dpi=512)
    plt.close()

# Figure 2B
flist = ['ape_cons94_clean', 'ape_cons94_exonic']
clist = ['darkorange', 'fuchsia']
llist = ['all conserved', 'exonic conserved']
sname = 'fig_2B'
sfldr = 'mainfigs'
compare_rsq(flist, llist, clist, sname, sfldr, 'B')


#%% INDIVIDUAL COLLATED PLOT
def collated_plot(flist, llist, clist, sname, sfldr):
    obs_flag = False
    plt.figure(figsize=(4, 2.5))
    plt.subplots_adjust(left=0.12, bottom=0.15, right=1, top=1)

    # COLLATED PLOT
    for (i, fl) in enumerate(flist):
        f_ape = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(fl)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        obs = pi / cnts
        prd = pr / cnts
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        if not obs_flag:
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.1)
            obs_flag = True
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.1)

    plt.xlabel('distance to nearest amino acid substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'heterozygosity ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    plt.legend(prop=dict(size=9), loc='lower right', handletextpad=0.5,
               borderaxespad=0.05)
    f_save = final_dir + '/{}/{}.collated.png'.format(sfldr, sname)
    plt.savefig(f_save, dpi=512)
    plt.close()


#%% COLLATED AND R^2 COMPARISON TO MCVICKER

def collate_and_rsq_mcvicker():

    flist = ['cadd93', 'mcvicker']
    llist = ['our map', 'McVicker']
    clist = ['darkorange', 'mediumpurple']
    sname = 'fig_S20'
    sfldr = 'sfigs'
    # fletters = ('B', 'C')

    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)

    obs_flag = False
    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.06, bottom=0.2, right=1, top=0.98, wspace=0.3)
    axes_lbl_size = 9
    xtck_size = 8
    ytck_size = 8
    leg_size = 8
    fnt_size = 9

    plt.rc('font', size=fnt_size)  # controls default text sizes
    plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    plt.rc('legend', fontsize=leg_size)  # legend fontsize

    # COLLATED PLOT
    plt.subplot(1,3,(1,2))
    # plt.subplots_adjust(wspace=0.2, left=0.12)
    # plt.subplot(1, 3, (1,2))
    for (i, rdir) in enumerate(flist):
        f_ape = root_dir + '/result/final_files/{}/nonsyn.collate.5.00e-03width.npy'.format(
            rdir)
        # f_mcv = root_dir + '/result/collate/YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy'
        bins, div, pi, pr, cnts = np.load(f_ape).T
        # _, _, _, mcv, mcv_cnts = np.load(f_mcv).T
        # bins, pi, div, pred, cnts = np.load(f_in).T
        obs = pi / cnts
        prd = pr / cnts
        # mcv /= mcv_cnts
        # y_max = pi0
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        # mcv /= y_max
        # plt.title(rdir)
        # pc = '{}%'.format(100-int(rdir.split('_')[1][4:]))
        if not obs_flag:
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.1)
            obs_flag = True
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.1)

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'heterozygosity ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    # plt.ylim(0.83, 1.13)

    # R^2 PLOT
    dall = 'cadd93'
    f_mcv = root_dir + '/result/final_files/mcvicker/rsq.log'
    f_cad = root_dir + '/result/final_files/{}/rsq.log'.format(dall)
    w, mcv = np.loadtxt(f_mcv)[:16].T
    ape = np.loadtxt(f_cad)[:16, 1]

    w = np.log10(w)
    # plt.figure(figsize=(3.25, 2.5))
    plt.subplot(133)
    plt.plot(w, ape, label='our map', marker='o', lw=0,
             color='darkorange')
    plt.plot(w, mcv, label='McVicker', marker='o', lw=0, color='mediumpurple')
    # plt.plot(w, ex, label='our method (exon only)', marker='o', lw=0,
    #          color='fuchsia')

    plt.xlabel('window size (bp)')
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4, 5, 6]])
    plt.ylim(0, 0.65)
    plt.ylabel('variance explained ' + r'$(R^2)$')
    # plt.legend(prop=dict(size=9), loc='upper left')
    plt.text(0.01, 0.95, 'B', transform=plt.gcf().transFigure)
    plt.text(0.67, 0.95, 'C', transform=plt.gcf().transFigure)

    # f_save = root_dir + '/result/final_files/sfigs/fig_20.rsq.mcvicker.png'
    # plt.savefig(f_save, dpi=256)
    # plt.close()

    # plt.text(0.01, 0.93, fletters[0], transform=plt.gcf().transFigure)
    # plt.text(0.68, 0.93, fletters[1], transform=plt.gcf().transFigure)

    f_save = sdir + '/{}.collated.rsq.combined.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% PREDICTED VS. OBSERVED MCVICKER


def sorted_mcvicker_plot():
    """large formated observed vs. predicted plot comparing to mcvicker"""
    flist = ['cadd93', 'mcvicker']
    llist = ['our map', 'McVicker']
    clist = ['darkorange', 'mediumpurple']
    sname = 'fig_S20'
    sfldr = 'sfigs'
    sdir = root_dir + '/result/final_files/' + sfldr

    plt.figure(figsize=(3,3))
    plt.subplots_adjust(left=0.15, top=1, right=0.99)
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    if 'McVicker' in llist:
        plt.text(0.5, 1.05, 'without linked selection', ha='center',
                 va='center')
    else:
        plt.text(0.75, 1.05, 'without linked selection', ha='center',
                 va='center')
    ymin = 1
    for (i, rdir) in enumerate(flist):
        # load data
        fdir = root_dir + '/result/final_files/{}/'.format(rdir)
        slist = [f for f in os.listdir(fdir) if ('predsort' in f) and
                 (f.endswith('.txt'))]
        f_name = fdir + slist[0]
        rst = cst_from_fldr(rdir)
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        f_in = f_name

        # normalize by pi0
        div, pi, pred = np.loadtxt(f_in).T
        pi /= pi0
        pred /= pi0
        ymin = min(ymin, pred.min())
        plt.plot(pred, pi, color='darkslategray', marker='o', ms=5,
                 markerfacecolor='None', markeredgecolor=clist[i],
                 markeredgewidth=0.9, lw=0, alpha=0.75)

    if 'McVicker' in llist:
        plt.ylim(0., 1.15)
        plt.xlim(0., 1.02)
        xtick = np.arange(0.1, 1.01, 0.1)
        ytick = np.arange(0.1, 1.11, 0.1)
    else:
        plt.ylim(0.48, 1.15)
        plt.xlim(0.48, 1.02)
        xtick = np.arange(0.5, 1.01, 0.1)
        ytick = np.arange(0.5, 1.11, 0.1)

    plt.xticks(xtick, y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick, x=0.02)
    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.legend(loc='lower right')

    plt.text(0.01, 0.95, 'D', transform=plt.gcf().transFigure)
    # plt.text(0.68, 0.93, fletters[1], transform=plt.gcf().transFigure)

    f_save = sdir + '/{}.sort.vs.mcvicker.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


sorted_mcvicker_plot()
#%% R^2 MCVICKER


def rsq_mcvicker():
    dall= 'cadd93'
    # dex = 'ape_cons94_exonic'
    f_mcv = root_dir + '/result/final_files/mcvicker/rsq.log'
    f_cad = root_dir + '/result/final_files/{}/rsq.log'.format(dall)
    # f_exonly = root_dir + '/result/final_files/{}/rsq.log'.format(dex)
    w, mcv = np.loadtxt(f_mcv)[:16].T
    ape = np.loadtxt(f_cad)[:16, 1]
    # ex = np.loadtxt(f_exonly)[:16, 1]

    w = np.log10(w)
    plt.figure(figsize=(3.25,2.5))
    plt.subplots_adjust(top=0.98, right=0.99, left=0.2, bottom=0.18)
    plt.plot(w, ape, label='our method', marker='o', lw=0,
             color='darkorange')
    plt.plot(w, mcv, label='McVicker', marker='o', lw=0, color='mediumpurple')
    # plt.plot(w, ex, label='our method (exon only)', marker='o', lw=0,
    #          color='fuchsia')

    plt.xlabel('window size (bp)')
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4, 5, 6]])
    plt.ylim(0, 0.65)
    plt.ylabel('variance explained ' + r'$(R^2)$')
    plt.legend(prop=dict(size=9), loc='upper left')
    plt.text(0.01, 0.95, 'D', transform=plt.gcf().transFigure)

    f_save = root_dir + '/result/final_files/sfigs/fig_20.rsq.mcvicker.png'
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% R^2 COMPARED TO MEAN PI AND COEFFICIENT OF VARIATION


def rsq_vs_meanpi():
    """plot R^2 values at 1Mb vs. mean diversity per population"""
    f_save = root_dir + '/result/final_files/sfigs/rsq_vs_meanpi.png'
    pops = 'YRI CEU TSI CHB JPT MXL'.split()
    cols = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
     'rosybrown', 'darkturquoise']
    rsq, mpi = [], []
    plt.figure(figsize=(3, 3))
    plt.subplots_adjust(right=1, left=0.22, top=1, bottom=0.15)
    for (p, c) in zip(pops, cols):
        fldr = '{}_cadd93'.format(p)
        if p == 'YRI':
            fldr = 'cadd93'
        pi = cst_from_fldr(fldr).stat.meanpi * 1e3
        mpi.append(pi)
        w, r = rsq_from_fldr(fldr).T
        rmb = r[w==1e6][0]
        rsq.append(rmb)
        # plt.text(pi, rmb, p)
        plt.plot(pi, rmb, marker='o', ms=10, color=c, lw=0, alpha=0.75, label=p)
    xmin, xmax = min(mpi), max(mpi)
    xrng = xmax - xmin
    plt.xlim(xmin - 0.1 * xrng, xmax + 0.1 * xrng)
    ymin, ymax = min(rsq), max(rsq)
    yrng = ymax-ymin
    plt.ylim(ymin - 0.1 * yrng, ymax + 0.1 * yrng)
    plt.ylabel(r'$R^2$ (1Mb window)')
    plt.xlabel(r'$\bar{\pi}\ (\times 10^3)$')
    plt.legend()
    plt.savefig(f_save, dpi=512)
    plt.close()


def rsq_vs_meanpi_all():
    """plot R^2 values at 1Mb vs. mean diversity per population"""
    f_save = root_dir + '/result/final_files/sfigs/cadd93_rsq_vs_meanpi_all.png'
    pfile = root_dir + '/data/snps/grouped_pops.txt'
    pops, grps = [], []
    with open(pfile, 'r') as f:
        for line in f:
            pp, gp = line[:-1].split()
            pops.append(pp)
            grps.append(gp)

    cols = {'AFR':'darkorange', 'EUR':'steelblue', 'AMR':'fuchsia',
            'EAS':'purple', 'SAS':'darkturquoise'}
    rsq, mpi = [], []
    seen_grp = []
    plt.figure(figsize=(3, 3))
    plt.subplots_adjust(right=1, left=0.25, top=1, bottom=0.2)
    for (p, g) in zip(pops, grps):
        # fldr = '{}_ape_cons94_clean'.format(p)
        fldr = '{}_cadd93'.format(p)
        if p == 'YRI':
            fldr = 'cadd93'
        pi = cst_from_fldr(fldr).stat.meanpi * 1e3
        mpi.append(pi)
        w, r = rsq_from_fldr(fldr).T
        rmb = r[w==1e6][0]
        rsq.append(rmb)
        c = cols[g]
        if g in seen_grp:
            plt.plot(pi, rmb, marker='o', ms=5, color=c, lw=0, alpha=0.75)
        else:
            plt.plot(pi, rmb, marker='o', ms=5, color=c, lw=0, alpha=0.75,
                     label=g)
        seen_grp.append(g)
    xmin, xmax = min(mpi), max(mpi)
    xrng = xmax - xmin
    plt.xlim(xmin - 0.1 * xrng, xmax + 0.1 * xrng)
    # ymin, ymax = min(rsq), max(rsq)
    # yrng = ymax-ymin
    # plt.ylim(ymin - 0.1 * yrng, ymax + 0.1 * yrng)
    ytick = np.arange(0.425, 0.58, 0.025)
    plt.yticks(ytick)
    plt.ylim(0.405, 0.595)
    plt.ylabel(r'$R^2$ (1Mb window)')
    plt.xlabel(r'$\bar{\pi}\ (\times 10^3)$')
    plt.legend(loc='upper left', frameon=True, framealpha=0.5,
               facecolor='white')
    plt.savefig(f_save, dpi=512)
    plt.close()


def cv_vs_cd():
    """plot coefficient of variation of pi at 1Mb for each pop"""
    pfile = root_dir + '/data/snps/grouped_pops.txt'
    cvfile = root_dir + '/result/final_files/sfigs/cv_cadd93_all.txt'
    cols = {'AFR': 'darkorange', 'EUR': 'steelblue', 'AMR': 'fuchsia',
            'EAS': 'purple', 'SAS': 'darkturquoise'}

    # create list of pops ordered by group and pop:group dict
    pdict = {}
    pops, grps = [], []
    with open(pfile, 'r') as f:
        for line in f:
            pp, gp = line[:-1].split()
            pdict[pp] = gp
            pops.append(pp)
            grps.append(gp)

    # create dict of pop:CV
    cv_dict = {}
    with open(cvfile, 'r') as f:
        for line in f:
            pp, cv = line[:-1].split()
            cv_dict[pp] = float(cv)

    # plot results ordered by group
    plt.figure(figsize=(3, 3))
    plt.subplots_adjust(right=1, top=1, left=0.25, bottom=0.2)
    seen_groups = []
    cv_list = []
    rm_list = []
    for (i, pp) in enumerate(pops):
        gp = pdict[pp]
        col = cols[gp]
        cv = cv_dict[pp]
        if pp == 'YRI':
            # fldr = 'ape_cons94_clean'
            fldr = 'cadd93'

        else:
            # fldr = '{}_ape_cons94_clean'.format(pp)
            fldr = '{}_cadd93'.format(pp)

        w, r = rsq_from_fldr(fldr).T
        rmb_1 = r[w==1e6][0]
        # r_1.append(rmb_1)
        # print i, cv
        cv = 1.0 / (cv**2)
        cv_list.append(cv)
        rm_list.append(rmb_1)
        if gp in seen_groups:
            plt.plot(cv, rmb_1, color=col, marker='o', ms=5, lw=0, alpha=0.75)
        else:
            plt.plot(cv, rmb_1, color=col, label=gp, marker='o', ms=5, lw=0,
                     alpha=0.75)
        seen_groups.append(gp)

    lrg = linregress(cv_list, rm_list)

    # variance in B at 1Mb
    # bvar = 0.00793660243859
    bvar = 0.0329278992558
    # plt.plot([11, 19], [10*lrg.slope+lrg.intercept, 20*lrg.slope+lrg.intercept], ls='--', label=r'y=mx+b')
    plt.plot([11, 19], [10*bvar, 20*bvar], ls='--', label=r'y=(Var(B)/mean(B)^2)x')
    plt.ylabel(r'$R^2$ (1Mb window)')
    # ytick = np.arange(0.425, 0.58, 0.025)
    # plt.yticks(ytick)
    # plt.ylim(0.405, 0.595)
    # xtick = np.arange(0.23, 0.291, 0.01)
    # plt.xticks(xtick, ['', 0.24, '', 0.26, '', 0.28, ''])
    # plt.xlim(0.229, 0.291)
    # plt.xticks(range(len(pops)), pops, rotation=90)
    # plt.xlim(-0.5, 25.5)
    plt.xlabel(r'$1/(CoV)^2$ for $\bar{\pi}$ (1Mb window)')
    plt.legend(loc='upper left', frameon=True, framealpha=0.5,
               facecolor='white')

    f_save = final_dir + '/sfigs/cadd93_cv_vs_cd_1Mb_modified_bvar.png'
    plt.savefig(f_save, dpi=512)
    plt.close()

# rsq_vs_meanpi_all()
cv_vs_cd()

#%% R^2 COMPARING POP SPECIFIC VS. OUTSIDE B MAPS


def rsq_two_maps():
    """plot R^2 values at 1Mb using native vs. YRI map for each pop"""
    f_save = root_dir + '/result/final_files/sfigs/native_vs_YRI_map.png'
    pops = 'CEU TSI CHB JPT MXL'.split()
    cols = ['red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
    r_1, r_2 = [], []
    plt.figure(figsize=(4, 3))
    plt.subplots_adjust(right=1, top=1, left=0.17, bottom=0.17)
    for (p, c) in zip(pops, cols):
        fldr = '{}_cadd93'.format(p)
        w, r = rsq_from_fldr(fldr).T
        rmb_1 = r[w==1e6][0]
        r_1.append(rmb_1)
        f_r = root_dir + '/result/final_files/{}/rsq.YRI.map.log'.format(fldr)
        w, r = np.loadtxt(f_r).T
        rmb_2 = r[w == 1e6][0]
        r_2.append(rmb_2)
        plt.plot(rmb_1, rmb_2, marker='o', ms=10, color=c, lw=0, alpha=0.75,
                 label=p)

    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    xmin, xmax = min(r_1), max(r_1)
    xrng = xmax - xmin
    # plt.xlim(xmin - 0.1 * xrng, xmax + 0.1 * xrng)
    plt.xlim(0.42, 0.4625)
    # plt.xticks([0.425, 0.43, 0.435, 0.44, 0.445, 0.45, 0.455, 0.46])
    ymin, ymax = min(r_2), max(r_2)
    yrng = ymax-ymin
    # plt.ylim(ymin - 0.1 * yrng, ymax + 0.1 * yrng)
    plt.ylim(0.42, 0.4625)
    # plt.yticks([0.425, 0.43, 0.435, 0.44])
    plt.ylabel(r'$R^2$ YRI B-map (1Mb window)')
    plt.xlabel(r'$R^2$ population B-map (1Mb window)')
    plt.legend()
    plt.savefig(f_save, dpi=512)
    plt.close()


def rsq_two_maps_all(refmap='YRI'):
    """plot R^2 values at 1Mb using native vs. ref map map for each pop"""
    save_str = '/result/final_files/sfigs/native_vs_{}_map_all.png'
    f_save = root_dir + save_str.format(refmap)
    pfile = root_dir + '/data/snps/grouped_pops.txt'
    pops, grps = [], []
    with open(pfile, 'r') as f:
        for line in f:
            pp, gp = line[:-1].split()
            pops.append(pp)
            grps.append(gp)

    cols = {'AFR':'darkorange', 'EUR':'steelblue', 'AMR':'fuchsia',
            'EAS':'purple', 'SAS':'darkturquoise'}
    r_1, r_2 = [], []
    seen_grp = []
    plt.figure(figsize=(4, 3))
    plt.subplots_adjust(right=1, top=1, left=0.18, bottom=0.17)
    for (p, g) in zip(pops, grps):
        if p == refmap:
            continue
        if p == 'YRI':
            # fldr = 'ape_cons94_clean'
            fldr = 'cadd93'
        else:
            # fldr = '{}_ape_cons94_clean'.format(p)
            fldr = '{}_cadd93'.format(p)

        w, r = rsq_from_fldr(fldr).T
        rmb_1 = r[w==1e6][0]
        r_1.append(rmb_1)
        f_r = root_dir + '/result/final_files/{}/rsq.{}.map.log'.format(fldr,
                                                                        refmap)
        w, r = np.loadtxt(f_r).T
        rmb_2 = r[w == 1e6][0]
        r_2.append(rmb_2)
        c = cols[g]
        if g in seen_grp:
            plt.plot(rmb_1, rmb_2, marker='o', ms=5, color=c, lw=0, alpha=0.75)
        else:
            plt.plot(rmb_1, rmb_2, marker='o', ms=5, color=c, lw=0, alpha=0.75,
                     label=g)
        seen_grp.append(g)

    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    xy_lim = (0.415, 0.585)
    xy_ticks = np.arange(0.425, 0.58, 0.025)
    plt.xticks(xy_ticks)
    plt.xlim(xy_lim)
    plt.yticks(xy_ticks)
    plt.ylim(xy_lim)
    plt.ylabel(r'$R^2$ {} B-map (1Mb window)'.format(refmap))
    plt.xlabel(r'$R^2$ population B-map (1Mb window)')
    plt.legend()
    plt.savefig(f_save, dpi=512)
    plt.close()


def rsq_two_maps_all_2(refmap='YRI'):
    """plot R^2 values at 1Mb using native vs. ref map map for each pop"""
    save_str = '/result/final_files/sfigs/native_vs_{}_map_all_2.png'
    f_save = root_dir + save_str.format(refmap)
    pfile = root_dir + '/data/snps/grouped_pops.txt'
    pops, grps = [], []
    with open(pfile, 'r') as f:
        for line in f:
            if refmap in line:
                continue
            pp, gp = line[:-1].split()
            pops.append(pp)
            grps.append(gp)

    cols = {'AFR':'rosybrown', 'EUR':'steelblue', 'AMR':'fuchsia',
            'EAS':'purple', 'SAS':'darkturquoise'}
    r_1, r_2 = [], []
    seen_grp = []
    plt.figure(figsize=(5, 3.33))
    plt.subplots_adjust(right=1, top=0.99, left=0.18, bottom=0.2)
    i = 0
    for (p, g) in zip(pops, grps):
        if p == refmap:
            continue
        if p == 'YRI':
            fldr = 'ape_cons94_clean'.format(p)
        else:
            fldr = '{}_ape_cons94_clean'.format(p)
        w, r = rsq_from_fldr(fldr).T
        rmb_1 = r[w==1e6][0]
        r_1.append(rmb_1)
        f_r = root_dir + '/result/final_files/{}/rsq.{}.map.log'.format(fldr,
                                                                        refmap)
        w, r = np.loadtxt(f_r).T
        rmb_2 = r[w == 1e6][0]
        r_2.append(rmb_2)
        c = cols[g]
        if g in seen_grp:
            plt.plot(i, rmb_1-rmb_2, marker='o', ms=5, color=c, lw=0, alpha=0.75)
        else:
            plt.plot(i, rmb_1-rmb_2, marker='o', ms=5, color=c, lw=0, alpha=0.75,
                     label=g)
        seen_grp.append(g)
        i += 1

    plt.plot([-1, 25], [0, 0], color='darkslategray', ls='--', alpha=0.65, lw=0.75)
    # xmin, xmax = min(r_1), max(r_1)
    # xrng = xmax - xmin
    # plt.xlim(xmin - 0.1 * xrng, xmax + 0.1 * xrng)
    # # plt.xlim(0.42, 0.4625)
    # plt.xticks([0.425, 0.43, 0.435, 0.44, 0.445, 0.45, 0.455, 0.46])
    # diff = [r1-r2 for r1, r2 in zip(r_1, r_2)]
    # ymin, ymax = min(diff), max(diff)
    # yrng = ymax-ymin
    # plt.ylim(ymin - 0.1 * yrng, ymax + 0.1 * yrng)
    plt.ylim(-0.00125, 0.00475)
    # plt.yticks([0.425, 0.43, 0.435, 0.44])
    plt.ylabel(r'population - {} $R^2$ (1Mb window)'.format(refmap))
    plt.xticks(range(len(pops)), pops, rotation=90)
    plt.xlim(-0.5, 24.5)
    plt.xlabel(r'population')
    plt.legend(loc='upper right')
    plt.savefig(f_save, dpi=512)
    plt.close()


rsq_two_maps_all('YRI')

#%% COEFFICIENT OF VARIATION IN DIVERSITY AT 1MB WINDOWS


def coefficient_of_variation_mb():
    """plot coefficient of variation of pi at 1Mb for each pop"""
    f_save = root_dir + '/result/final_files/sfigs/coefficient_of_var_1Mb.png'
    pfile = root_dir + '/data/snps/grouped_pops.txt'
    cvfile = root_dir + '/result/final_files/sfigs/cv_all.txt'
    cols = {'AFR': 'rosybrown', 'EUR': 'steelblue', 'AMR': 'fuchsia',
            'EAS': 'purple', 'SAS': 'darkturquoise'}

    # create list of pops ordered by group and pop:group dict
    pdict = {}
    pops, grps = [], []
    with open(pfile, 'r') as f:
        for line in f:
            pp, gp = line[:-1].split()
            pdict[pp] = gp
            pops.append(pp)
            grps.append(gp)

    # create dict of pop:CV
    cv_dict = {}
    with open(cvfile, 'r') as f:
        for line in f:
            pp, cv = line[:-1].split()
            cv_dict[pp] = float(cv)

    # plot results ordered by group
    plt.figure(figsize=(5, 3.33))
    plt.subplots_adjust(right=0.98, top=0.99, left=0.18, bottom=0.2)
    seen_groups = []
    for (i, pp) in enumerate(pops):
        gp = pdict[pp]
        col = cols[gp]
        cv = cv_dict[pp]
        # print i, cv
        if gp in seen_groups:
            plt.plot(i, cv, color=col, marker='o', ms=5, lw=0)
        else:
            plt.plot(i, cv, color=col, label=gp, marker='o', ms=5, lw=0)
        seen_groups.append(gp)

    plt.ylabel(r'coefficient of variation for $\pi$ (1Mb window)')
    plt.xticks(range(len(pops)), pops, rotation=90)
    plt.xlim(-0.5, 25.5)
    plt.xlabel(r'population')
    plt.legend(loc='upper left', frameon=True, framealpha=0.5,
               facecolor='white')
    plt.savefig(f_save, dpi=512)
    plt.close()


#%% FIGURE S8 A-C
spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()
flist = ['{}_cons94_gmask_fixed'.format(sp) for sp in spec]
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
         'rosybrown', 'darkturquoise']
llist = ['ape (4)', 'prim (8)', 'pros (12)', 'supr (25)', 'laur (50)',
         'mamm (61)', 'vert (98)']
sname = 'fig_S8'
sfldr = 'sfigs'
xlab = 'phylogenetic depth'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, llist, sname, sfldr, xlab, rotation=45)
#%% FIGURE S8 D-F
spec = 'ape euarchontoglires fish'.split()
flist = ['{}_cons94_gmask_fixed'.format(sp) for sp in spec]
clist = ['darkorange', 'fuchsia', 'darkturquoise']
llist = ['ape (4)', 'supr (25)', 'vert (98)']
sname = 'fig_S8_v2'
sfldr = 'sfigs'
# chr1_diversity(flist, llist, clist, sname, sfldr, fletter='D')
collate_and_sort(flist, llist, clist, sname, sfldr, fletters=('E', 'F'))
#%% FIGURE S10 A-C
pcons = range(98, 89, -1)
flist = ['ape_cons{}_clean'.format(pc) for pc in pcons]
clist = 'darkorange crimson lightcoral indianred fuchsia firebrick darkred ' \
        'red darkturquoise'.split()
llist = ['{}%'.format(100-pc) for pc in pcons]
sname = 'fig_S10'
sfldr = 'sfigs'
xlab = 'fraction of conserved sites'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, llist, sname, sfldr, xlab, 45)
#%% FIGURE S10 D-F
pcons = [98, 94, 90]
flist = ['ape_cons{}_clean'.format(pc) for pc in pcons]
clist = 'darkorange fuchsia darkturquoise'.split()
llist = ['{}%'.format(100-pc) for pc in pcons]
sname = 'fig_S10'
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S11 A-E
flist = ['ape_cons94_clean', 'genic']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'genic']
sllist = ['cds', 'splicing']
sname = 'fig_S11'
sfldr = 'sfigs'
multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S12 A-E: CONSERVED STRATIFIED EX/NEX
flist = ['ape_cons94_clean', 'ape_exnex_94']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'conserved (ex/nex)']
# sllist = ['exonic', 'nonexonic']
sname = 'fig_S12'
sfldr = 'sfigs'
# multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S13 A-E: CONSERVED EX/NEX INDIVIDUALLY
flist = ['ape_cons94_clean', 'ape_cons94_exonic', 'ape_cons94_nonexonic']
clist = ['darkorange', 'dodgerblue', 'firebrick']
llist = ['conserved', 'conserved exons', 'conserved other']
sname = 'fig_S13'
sfldr = 'sfigs'
# multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S15 A-E: COMPARISONS TO CHROMHMM
flist = ['ape_cons94_clean', 'huvec_chromhmm']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'Huvec']
sllist = ['transcribed', 'promoter', 'enhancer']
sname = 'fig_S15'
sfldr = 'sfigs'
# multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S16 A-E: COMPARISONS TO CADD
flist = ['ape_cons94_clean', 'cadd94']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'CADD']
sllist = ['transcribed', 'promoter', 'enhancer']
sname = 'fig_S16'
sfldr = 'sfigs'
# multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S18 A-E: CS MODELS
flist = ['cadd93', 'nonsyn', 'cadd90_CS_only', 'cadd98_CS_only']
clist = ['darkorange', 'dodgerblue', 'firebrick', 'mediumpurple']
llist = ['BS', 'NS', '10%', '2%']
sllist = ['other', 'NS', 'other']
sname = 'fig_S18'
sfldr = 'sfigs'
bscs_joint(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE S20 A-E: COMPARISONS TO MCVICKER
flist = ['cadd93', 'mcvicker']
clist = ['darkorange', 'mediumpurple']
llist = ['our map', 'McVicker']
sname = 'fig_S20'
sfldr = 'sfigs'
# bscs_joint(flist, llist, sllist, clist, sname, sfldr)
combined_chr1(flist, llist, clist, sname, sfldr, fletter='A')
#%% FIGURE S23 A-C: COMPARING POPULATIONS
pops = 'CEU GIH JPT MXL'.split()
flist = ['cadd93'] + ['{}_ape_cons94_clean'.format(p) for p in pops]
clist = ['darkorange', 'steelblue', 'darkturquoise', 'purple', 'fuchsia']
llist = ['YRI'] + pops
sname = 'fig_S23'
sfldr = 'sfigs'
xlab = 'population'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, llist, sname, sfldr, xlab, rotation=45)
#%% FIGURE SX 1-5: COMPARING POPULATIONS
pops = 'CEU TSI CHB JPT MXL'.split()
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
         'rosybrown', 'darkturquoise']
sfldr = 'sfigs'
for (i, p) in enumerate(pops):
    f_name = '{}_cadd93'.format(p)
    s_name = 'fig_SX_{}'.format(p)
    combined_chr1([f_name], [p], [clist[i]], s_name, sfldr)
    collate_and_sort([f_name], [p], [clist[i]], s_name, sfldr)
#%% FIGURE SY A-C: COMPARING CADD %
pct = range(98, 89, -1)
flist = ['cadd{}'.format(p) for p in pct]
clist = 'crimson steelblue mediumvioletred fuchsia purple darkorange darkred ' \
        'gold darkturquoise'.split()
llist = ['{}%'.format(100-p) for p in pct]
sname = 'fig_SY'
sfldr = 'sfigs'
xlab = 'CADD threshold'
# rsq_3_size(flist, llist, sname, sfldr, xlab, rotation=45)
single_bsanno(flist, llist, clist, sname, sfldr)
#%% FIGURE SX D-F: COMPARING CADD % (SUBSET)
pct = [98, 93, 90]
flist = ['cadd{}'.format(p) for p in pct]
clist = ['crimson', 'darkorange', 'darkturquoise']
llist = ['{}%'.format(100-p) for p in pct]
sname = 'fig_SY'
sfldr = 'sfigs'
combined_chr1(flist, llist, clist, sname, sfldr, fletter='D')
collate_and_sort(flist, llist, clist, sname, sfldr, fletters=('E', 'F'))
#%% FIGURE SZ A-E
flist = ['cadd93', 'cadd93_CpG_filter']
clist = ['darkorange', 'dodgerblue']
llist = ['CADD 7%', 'CADD 7% non-CpG']
sllist = ['nan', 'nan']
sname = 'fig_SZ'
sfldr = 'sfigs'
multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
combined_chr1(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE ST C
flist = ['cadd93', 'cadd93_new_nu']
clist = ['darkorange', 'dodgerblue']
llist = ['CADD 7%', 'CADD 7% new nu']
sllist = ['nan', 'nan']
sname = 'fig_ST'
sfldr = 'sfigs'
multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
combined_chr1(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE SV C
flist = ['cadd93', 'cadd93_new_nu', 'cadd93_CpG', 'cadd93_BGC']
clist = ['darkorange', 'dodgerblue', 'purple', 'rosybrown']
llist = ['CADD 7%', 'CADD 7% new nu', 'CADD 7% non-CpG', 'CADD 7% non-BGC']
sllist = ['nan', 'nan']
sname = 'fig_SV'
sfldr = 'sfigs'
multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
combined_chr1(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
# %% FIGURE COMPARING FIXES TO COMPRESS WITHIN CADD
flist = ['cadd93', 'cadd93_updated_compress']
clist = ['darkorange', 'dodgerblue']
llist = ['old', 'new']
sllist = []
sname = 'update_compress'
sfldr = 'sfigs'
multi_bsanno(flist, llist, sllist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr)
collate_and_sort(flist, llist, clist, sname, sfldr)
#%% FIGURE SGG A-E: different mutation type filters
flist = ['cadd93', 'cadd93_align8', 'cadd93_cg_filter']
clist = ['darkorange', 'dodgerblue', 'mediumvioletred' ]
llist = ['all sites', 'aligned', 'filter C>G']
sname = 'fig_SGG'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, ['all', 'aligned', 'C>G'], sname, sfldr, 'filter')
combined_chr1(flist, llist, clist, sname, sfldr, 'D')
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
# %% FIGURE SGG A-E: different mutation type filters
flist = ['cadd93', 'cadd93_gmask_2']
clist = ['darkorange', 'dodgerblue']
llist = ['all sites', 'gmap mask']
sname = 'fig_SX_gmask'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, ['all', 'gmask'], sname, sfldr, 'filter')
combined_chr1(flist, llist, clist, sname, sfldr, 'D')
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
# %% Figure 2A
flist = ['ape_cons94_clean']
clist = ['darkorange']
llist = ['predicted']
sname = 'fig_2A'
sfldr = 'mainfigs'
combined_chr1(flist, llist, clist, sname, sfldr, 'A')
# %% Figure 2B
flist = ['ape_cons94_clean', 'ape_cons94_exonic']
clist = ['darkorange', 'fuchsia']
llist = ['all conserved', 'exonic conserved']
sname = 'fig_2B'
sfldr = 'mainfigs'
compare_rsq(flist, llist, clist, sname, sfldr, 'B')
#%%
flist = ['cadd93_gmask_000', 'cadd93_gmask_001', 'cadd93_gmask_005',
         'cadd93_gmask_01', 'cadd93_gmask_05', 'cadd93_gmask_10',]
clist = ['darkorange', 'dodgerblue', 'fuchsia', 'mediumvioletred', 'firebrick',
         'darkturquoise']
llist = '0.0 0.01 0.05 0.1 0.5 1.0'.split()
sname = 'fig_X_gmask_thresh'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, llist, sname, sfldr, 'cM masked', fletter='B', rotation=45)
# chr1_diversity(flist, llist, clist, sname, sfldr, 'A')
# collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))

#%%
flist = ['cadd93_bt25', 'cadd93_bt45', 'cadd93_bt50',
         'cadd93_bt55']
clist = ['darkorange', 'dodgerblue', 'fuchsia', 'mediumvioletred']
llist = ['B>0.25', 'B>0.45', 'B>0.50', 'B>0.55']
sname = 'fig_X_B_thresh'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, llist, sname, sfldr, 'B threshold')
combined_chr1(flist, llist, clist, sname, sfldr, 'D')
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
#%%
flist = ['cadd93_bt25', 'cadd93_bt45', 'cadd93_bt50',
         'cadd93_bt55']
clist = ['darkorange', 'dodgerblue', 'fuchsia', 'mediumvioletred']
llist = ['B>0.25', 'B>0.45', 'B>0.50', 'B>0.55']
sname = 'fig_X_B_thresh'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
rsq_3_size(flist, llist, sname, sfldr, 'B threshold')
combined_chr1(flist, llist, clist, sname, sfldr, 'D')
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
#%%
mnb_tkn = '603 520 444 375'.split()
bth_tkn = '000 200 500 600'.split()
clist = ['darkorange', 'dodgerblue', 'mediumvioletred', 'darkturquoise']
sfldr = 'sfigs'

# llist = [r'LLH$\geq$0', 'LLH>0.2', 'LLH>0.5', 'LLH>0.6']
llist = [r'B$\geq$0', r'B$\geq$0.2', r'B$\geq$0.5', r'B$\geq$0.6']
flist = ['cadd93_bth{}'.format(b) for b in bth_tkn]
sname = 'cadd_llh_bth'
single_bsanno(flist, llist, clist, sname, sfldr)
# chr1_diversity(flist, llist, clist, sname, sfldr, '')

#%%
mnb_tkn = '140 603 444'.split()
llist = [r'no thresh.', r'prec$\geq$0.2', r'prec$\geq$0.5',
         r'prec$\geq$0.6']
flist = ['cadd93_bth000'] + ['cadd93_mnb{}'.format(b) for b in mnb_tkn]
sname = 'cadd_mnb_bth'
single_bsanno(flist, llist, clist, sname, sfldr)
#%%
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
combined_chr1(flist, llist, clist, sname, sfldr, 'D')

flist = ['cadd93_mnb{}'.format(b) for b in mnb_tkn]
llist = ['prec B>0.50', 'prec B>0.55', 'prec B>0.60', 'prec B>0.65']
sname = 'cadd_prec_bth'
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
combined_chr1(flist, llist, clist, sname, sfldr, 'D')

flist = ['ape_cons94_bth{}'.format(b) for b in bth_tkn]
llist = ['LLH B>0.50', 'LLH B>0.55', 'LLH B>0.60', 'LLH B>0.65']
sname = 'ape_llh_bth'
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
combined_chr1(flist, llist, clist, sname, sfldr, 'D')

flist = ['ape_cons94_mnb{}'.format(b) for b in mnb_tkn]
llist = ['prec B>0.50', 'prec B>0.55', 'prec B>0.60', 'prec B>0.65']
sname = 'ape_prec_bth'
collate_and_sort(flist, llist, clist, sname, sfldr, ('E', 'F'))
combined_chr1(flist, llist, clist, sname, sfldr, 'D')

#%%
mnb_all = '119 890 678 590 512 442 378 318 263 165'.split()
bth_all = '200 300 400 450 500 550 600 650 700 800'.split()
mnb_tkn_all = [r'B$\geq${:.2f}'.format(float(b)*0.001) for b in bth_all]
llh_tkn_all = [r'LLH$\geq${:.2f}'.format(float(b)*0.001) for b in bth_all]

# bth_tkn = '000 200 500 600'.split()
# mnb_tkn = '119 512 378'.split()
clist = ['darkorange', 'dodgerblue', 'mediumvioletred', 'darkturquoise',
         'fuchsia', 'orangered', 'cyan', 'purple', 'rosybrown']

#
idx = [0, 4, 6]
# anno = 'ape_cons94'
anno = 'cadd94'
# anno = 'fish_cons94'

bthlist = ['{}_gmask_bth_000'.format(anno)]
bthlist += ['{}_gmask_bth_{}'.format(anno, bth_all[i]) for i in idx]
mnblist = ['{}_gmask_bth_000'.format(anno)]
mnblist += ['{}_gmask_mnb_{}'.format(anno, mnb_all[i]) for i in idx]
llist1 = ['none'] + [llh_tkn_all[i] for i in idx]
llist2 = ['none'] + [mnb_tkn_all[i] for i in idx]
sname1 = anno + '_bth'
sname2 = anno + '_mnb'
single_bsanno(mnblist, llist2, clist, sname2, 'sfigs')
#%%