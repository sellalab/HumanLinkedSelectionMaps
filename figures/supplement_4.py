__author__ = 'davidmurphy'

import os
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, cst_from_fldr
from figures.common_functions import format_panels, predict_loess
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


def get_loess_line(fldr, span, return_points=False):
    # load results in 2000 bins for LOESS plots
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 2000)
    div, pi, pred = np.loadtxt(sort_file).T
    # f_sort = fdir + 'sort_gc_cm_cn_il_n2000.txt'
    # div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    # # load sorted conserved separately (done using 0.05cM radius of sites)
    # if load_con:
    #     fcon = final_dir + '/{}/sorted.cons.n2000.txt'.format(fldr)
    #     cn = np.loadtxt(fcon)[:2000,2]

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # points = [pi, gc, cn, cm, il]
    points = [pi]
    wts = np.ones(shape=len(pred))
    loess_lines = []
    for a in points:
        # use this try/except sequence  for rare LOESS bugs where slightly
        # increasing the span parameter fixes the problem
        try:
            lo_line = predict_loess(pred, a, wts, span, pred)
            loess_lines.append(lo_line)
        except ValueError:
            try:
                lo_line = predict_loess(pred, a, wts, span*1.1, pred)
                loess_lines.append(lo_line)
            except ValueError:
                lo_line = predict_loess(pred, a, wts, span * 1.5, pred)
                loess_lines.append(lo_line)


    if return_points:
        return pred, loess_lines, points
    else:
        return pred, loess_lines


#%%
def print_stats(flist):
    for fldr in flist:
        print '\n{}\n---------'.format(fldr)
        # get 100 points data
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 100)
        div, pi, pred = np.loadtxt(sort_file).T

        rst = cst_from_fldr(fldr)
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        # mean reduction
        print 'div reduction: {:.2f}'.format((1-(np.mean(pi)/pi0))*100)
        # reduction in bottom 10%
        print 'bottom 10% reduction: {:.2f}'.format((1-np.mean(pi[0:10]/pi0))*100)
        # reduction in top 10%
        print 'top 10% reduction: {:.2f}'.format((1-np.mean(pi[90:]/pi0))*100)
        print 'top 10% reduction (trimmed): {:.2f}'.format((1-np.mean(pi[90:98]/pi0))*100)
        # print pi[90:]/pi0


flist = ['fish_cons94_new', 'cadd94_gmask_v1.6_without_bstat']
print_stats(flist)
#%% COMPARE SINGLE BS ANNO RESULTS
def single_bsanno(flist, llist, clist, sname, sfldr, letters=('a', 'b')):
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
    ncol = 1 + len(flist) / 5

    # create the plot
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.066, wspace=1.75, right=0.985, bottom=0.15,
                        top=0.91)

    ### 1. DFE plot
    ax1 = plt.subplot(1, 9, (1,4))  #1, 9, (1,4)
    format_panels(ax1)
    #
    # plt.title('{} (i)'.format(letters[0]), loc='center', y=0.985,
    #           fontweight='bold')
    plt.title('{}'.format(letters[0]), loc='left', y=0.985,
              fontweight='bold')
    plt.title('(i)', loc='center', y=0.985, fontweight='bold')
    ymax = 0
    for (i, df) in enumerate(pmf_list):
        assert len(df) == 1
        lbl = llist[i]
        df = df[0]
        ymax = max(df.max(), ymax)
        plt.bar(xi + s, df, w, label=lbl, color=clist[i], align='edge')
        s += w

    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=0)
    plt.ylim(0, 1.8*ymax)
    plt.yticks(x=0.03)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.035, fontsize=9)
    plt.xlabel('deleterious fitness effect', labelpad=1)
    plt.legend(loc='upper left', ncol=ncol)
    # make small ticks between selection coefs
    for xx in xi[:-1]:
        plt.axvline(xx+0.5, 0, 0.2, lw=1, ls='--', color='k')

    ### 2. plot udel
    ax2 = plt.subplot(1,9,5)
    format_panels(ax2)

    plt.title('(ii)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=0)
    # plt.ylabel('mutation rate per site ' +r'($\mathrm{\times u_0^{-1}}$)', labelpad=0)
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
    plt.yticks(x=0.3)

    ### 3. plot pi0
    ax3 = plt.subplot(1,9,6)
    format_panels(ax3)

    plt.title('(iii)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=1)
    plt.yticks(x=0.3)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    ### 4. plot CLL
    ax4 = plt.subplot(1,9,7)
    format_panels(ax4)

    plt.title('(iv)', loc='center', y=0.985, fontweight='bold')
    # plt.ylabel(r'$\mathrm{-\Delta CLL\ per\ site\ \times10^5}$', labelpad=3)
    plt.ylabel(r'$\mathrm{-\Delta}$' + 'CLL per site ' + r'$\times10^5$',
               labelpad=3)

    plt.yticks(x=0.3)
    s = -w * n / 2.0
    clh_list = [1e5 * (cl / min(clh_list) - 1) for cl in clh_list]
    for (i, cl) in enumerate(clh_list):
        print '{:<20} -CLL {}'.format(flist[i], cl)
        plt.bar(0+s, cl, w, color=clist[i])
        s += w
    if max(clh_list) < 0.1:
        plt.ylim(0, 0.11)
    # plt.ylim(min(clh_list), max(clh_list))
    plt.xticks([])

    ### 5. R^2 PLOT
    ax5 = plt.subplot(1,9, (8,9))
    format_panels(ax5)

    # plt.title('B', loc='left', y=0.97)
    plt.title(letters[1], loc='left', y=0.985, fontweight='bold')
    plt.ylabel(r'variance explained ($R^2$)', labelpad=1)
    plt.yticks(x=0.08)
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    ws = ws[wmsk]
    for (i, rs) in enumerate(rsq_list):
        # save the 1Mb rsq vals
        rsq = rs[:,1][wmsk]
        plt.plot(ws, rsq, color=clist[i], marker='o',
                 linewidth=0, alpha=0.75, label=llist[i], ms=4)
    plt.xlabel('window size (bp)', labelpad=1)
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in xtck], y=0.035)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    plt.ylim(0.05, 0.65)

    # plt.legend(prop=dict(size=9), loc='upper left')

    # save the plot
    save_lbl = sname + '.parameters'
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)
    fsave = sdir + '/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=512)
    plt.close()


#%% R^2 ACROSS 3 SCALES
def rsq_3_size(flist, llist, sname, xlab, rotation=0, fletter='c'):
    # create string format for rsq files
    fstr = root_dir + '/result/final_files/{}/rsq.log'

    # keep just 3 window sizes: 62.5K, 250K, 1M
    keep_wins = [62500, 2.5e5, 1e6]
    wlab = ['62.5 kb', '250 kb', '1 Mb']
    rsq_list = []
    for fl in flist:
        # load rsq and window size
        win, rsq = np.loadtxt(fstr.format(fl)).T
        # keep rsq at the subset of windows needed
        idx = np.in1d(win, keep_wins)
        rsq_list.append(rsq[idx])

    # convert rsq values to array of columns for each size
    r = np.array(rsq_list)

    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.105, bottom=0.21, wspace=0.3, right=0.99,
                        top=0.88)

    # plot subplot for each window size
    for pidx in range(1, 4):
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
        plt.title(wlab[pidx-1]+ ' windows', y=0.96)
        plt.plot(xi, ri, lw=0, marker='o', ms=6, color='k',
                 alpha=0.8)
        plt.xlabel(xlab, labelpad=1)
        plt.xticks(xi, [l[:9] for l in llist], rotation=rotation, y=0.06,
                   ha='center')
        plt.yticks(ytck, x=0.05)
        plt.ylim(rmin - 0.5*r_interval, rmax + 0.5*r_interval)
        if pidx == 1:
            ylab = r'variance explained ($R^2$)'
            plt.ylabel(ylab, labelpad=1)

    # figure letter
    plt.text(0.01, 0.93, fletter, transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/{}.3range.rsq.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% COMRPARE MULTIPLE CHR1 PREDICTIONS
def combined_chr1(flist, llist, clist, sname, fletter='d'):
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/sfigs'
    obs_flag = False
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.07, bottom=0.155, right=0.995, top=0.995)
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
            plt.plot(xi, obs, label='observed', color='darkslategray', lw=2)
            obs_flag = True

        plt.plot(xi, prd, color=clist[i], lw=1.5, alpha=0.8, label=llist[i])

    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(-0.4, 1.9)
    plt.ylim(0, 1.9)
    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x / 25) % 2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])
    loc = 'lower center'
    ncol = len(llist)+1
    if 'McVicker' in llist:
        loc = 'upper center'
        plt.ylim(0, 2)
        # ncol = 2
    plt.legend(loc=loc, ncol=ncol, frameon=1,
               framealpha=0.75, facecolor='white', handlelength=0.8,
               borderaxespad=0.3, columnspacing=0.8)

    # figure letter
    plt.text(0.01, 0.93, fletter, transform=plt.gcf().transFigure,
             fontweight='bold')
    f_save = sdir + '/{}.chr1.combined.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% COLLATED PLOT AND MULTIPLE SORTED PREDICTIONS COMBINED
def collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'),
                     legend_on=False):

    obs_flag = False
    # plt.figure(figsize=(6.5, 2))
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.062, bottom=0.152, right=0.99,
                        top=0.98, wspace=0.25)

    # 1. COLLATED PLOT
    # plt.subplot(121)
    ax1 = plt.subplot(1, 3, (1,2))
    format_panels(ax1)
    for (i, rdir) in enumerate(flist):
        f_ape = root_dir + '/result/final_files/{}/nonsyn.collate.5.00e-03width.npy'.format(
            rdir)
        f_mcv = root_dir + '/result/collate/YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy'
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
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=2)
            obs_flag = True
        # plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.1)
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.5)

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=1)
    plt.xticks(y=0.04)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=1)
    plt.yticks(x=0.02)
    if legend_on:
        plt.legend(loc='lower right', ncol=1, frameon=1,
                   framealpha=0.75, facecolor='white', handlelength=0.8,
                   borderaxespad=0.3, columnspacing=0.8, prop=dict(size=9))
    # plot inset
    obs_flag = False
    axins = inset_axes(ax1, width="100%", height="100%", loc=6,
                       bbox_to_anchor=(0.021, 0.1, 0.35, 0.35),
                       bbox_transform=ax1.transAxes)
    format_panels(axins)
    for (i, fl) in enumerate(flist):
        f_ape = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(fl)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        msk = (bins >= -0.05) & (bins <= 0.05)
        obs = pi / cnts
        prd = pr / cnts
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        if not obs_flag:
            plt.plot(bins[msk], obs[msk], label='observed',
                     color='darkslategray', lw=1.5)
            obs_flag = True
        plt.plot(bins[msk], prd[msk], label=llist[i], color=clist[i], alpha=0.9,
                 lw=1.1, ls='-')
    plt.xticks([-0.05, 0.05], y=0.07, fontsize=8)
    plt.xlim(-0.05, 0.05)
    axins.set_yticks([])
    plt.ylim(0.81, 1.05)

    # 2. PREDICTED VS. OBSERVED PLOT
    # plt.subplot(122)
    ymin, xmin = 1, 1
    xmax, ymax = 1, 1
    ax2 = plt.subplot(133)
    format_panels(ax2)
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='k',
             ls='--', alpha=1)
    for (i, rdir) in enumerate(flist):
        span = 0.1
        fldr = flist[i]
        color = clist[i]
        # get 100 points data
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 100)
        div, pi, pred = np.loadtxt(sort_file).T

        rst = cst_from_fldr(fldr)
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        pi /= pi0
        pred /= pi0

        # get loess line and original points
        prlo, lolist = get_loess_line(fldr, span)
        # pilo, gclo, cnlo, cmlo, illo = lolist
        pilo = lolist[0]

        # plot horizontal line at y=1
        plt.axhline(y=1, color='k', alpha=0.8, ls='-')

        # plot predicted vs. observed
        plt.plot(pred, pi, marker='o', ms=2, color='white', lw=0)
        plt.plot(pred, pi, marker='o', ms=2, color=color, lw=0, alpha=0.5)

        # plot LOESS line
        # plt.plot(prlo, pilo, lw=2, color='white')
        plt.plot(prlo, pilo, lw=2, color=color, alpha=0.8)
        xmin = min(xmin, prlo.min(), pred.min())
        xmax = max(xmax, prlo.max(), pred.max())
        ymin = min(ymin, pilo.min(), pi.min())
        ymax = max(ymax, pilo.max(), pi.max())

        # if color == 'darkorange':
        #     plt.plot(prlo, pilo, lw=1.5, color=color)
    print('ylims {} {}'.format(ymin, ymax))
    print('xlims {} {}'.format(xmin, xmax))
    if 'McVicker' in llist:
        plt.ylim(0., 1.15)
        plt.xlim(0., 1.02)
        xtick = np.arange(0.1, 1.01, 0.1)
        ytick = np.arange(0.1, 1.11, 0.1)
    else:
        axmin = 0.48
        tckmin = 0.5
        # plt.ylim(0.48, 1.15)
        # plt.xlim(0.48, 1.02)
        # xtick = np.arange(0.5, 1.01, 0.1)
        # ytick = np.arange(0.5, 1.11, 0.1)
        plt.ylim(ymin-0.02, ymax+0.02)
        plt.xlim(xmin-0.02, xmax+0.02)
        xtick = np.arange(0.1, 1.01, 0.1)
        ytick = np.arange(0.1, 1.21, 0.1)

    plt.xticks(xtick, y=0.04)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=1)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=1)
    plt.yticks(ytick, x=0.04)
    plt.ylim(ymin - 0.02, ymax + 0.02)
    plt.xlim(xmin - 0.02, xmax + 0.02)
    plt.text(xmin, 1.03, 'without linked selection', ha='left',
             va='center', fontsize=10)
    # solve y=x rotation
    adlen = xmax+0.02-(ymin-0.02)
    oplen = adlen * (xmax-xmin+0.04) / (ymax-ymin+0.04)
    rot = np.arctan((oplen/adlen)) * (180.0 / np.pi)
    # rot = 45 * oplen/adlen
    print('rotation = {}'.format(rot))
    plt.text(0.76, 0.64, r'$y=x$', rotation=rot, ha='left', va='bottom',
             color='k', alpha=1)
    plt.text(0.005, 0.93, fletters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.665, 0.93, fletters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/{}.collated.sort.combined.png'.format(sname)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% MULTI BS ANNO
def multi_bsanno(flist, llist, subllist, clist, sname, sfldr, ytck):
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
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.065, wspace=1.75, right=0.985, bottom=0.15,
                        top=0.91)
    ### 1. DFE plot
    ax1 = plt.subplot(1, 9, (1, 4))
    format_panels(ax1)
    plt.title('a', loc='left', y=0.985,
              fontweight='bold')
    plt.title('(i)', loc='center', y=0.985, fontweight='bold')
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
                # zpad the DFE
                if len(dd) < 6:
                    num_z = 6-len(dd)
                    zpad = [0] * num_z
                    dd = np.concatenate((zpad, dd))
                # lbl = llist[i] + ' ' + subllist[j]
                lbl = subllist[j]
                ymax = max(dd.max(), ymax)
                if sum(dd) < 0.01:
                    continue
                plt.bar(xi + s, dd, w, label=lbl, color=clist[i],
                        hatch=hatches[j]*5, align='edge')
                s += w
                j += 1

    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=0)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(ytck, x=0.03)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.035, fontsize=9)
    plt.xlabel('deleterious fitness effect', labelpad=1)
    plt.legend(loc='upper left', ncol=ncol)
    # make small ticks between selection coefs
    for xx in xi[:-1]:
        plt.axvline(xx + 0.5, 0, 0.2, lw=1, ls='--', color='k')

    ### 2. plot udel
    ax2 = plt.subplot(1, 9, 5)
    format_panels(ax2)
    plt.title('(ii)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=0)
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
    plt.yticks(ytck, x=0.3)

    ### 3. plot pi0
    ax3 = plt.subplot(1,9,6)
    format_panels(ax3)
    plt.title('(iii)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=1)
    plt.yticks(x=0.3)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    ### 4. plot CLH
    ax4 = plt.subplot(1, 9, 7)
    format_panels(ax4)
    plt.title('(iv)', loc='center', y=0.985, fontweight='bold')
    # plt.ylabel(r'$\mathrm{-\Delta CLL\ per\ site\ \times10^5}$', labelpad=3)
    plt.ylabel(r'$\mathrm{-\Delta}$' + 'CLL per site ' + r'$\times10^5$',
               labelpad=3)
    plt.yticks(x=0.3)
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

    ### 5. R^2 PLOT
    ax5 = plt.subplot(1,9, (8,9))
    format_panels(ax5)
    plt.title('b', loc='left', y=0.985, fontweight='bold')
    plt.ylabel(r'variance explained $(R^2)$', labelpad=3)
    plt.yticks(x=0.08)
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    ws = ws[wmsk]
    for (i, rs) in enumerate(rsq_list):
        # save the 1Mb rsq vals
        rsq = rs[:,1][wmsk]
        plt.plot(ws, rsq, color=clist[i], marker='o',
                 linewidth=0, alpha=0.75, label=llist[i], ms=4)
    plt.xlabel('window size (bp)', labelpad=1)
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in xtck], y=0.035)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    plt.ylim(0.05, 0.65)
    # save the plot
    save_lbl = sname + '.parameters'
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)
    # plt.text(0.01, 0.93, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.77, 0.93, 'B', transform=plt.gcf().transFigure)
    fsave = sdir + '/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=512)
    plt.close()


#%% BS+CS PARAMETERS PLOT
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
    hatches = ['.', '/']*2

    # create the plot
    plt.figure(figsize=(6.5, 2.5))
    # plt.subplots_adjust(left=0.065, wspace=1.75, right=0.98, bottom=0.15,
    #                     top=0.91, hspace=0.5)
    plt.subplots_adjust(left=0.066, wspace=1.75, right=0.985, bottom=0.15,
                        top=0.91, hspace=0.5)

    ### 1a. BS DFE plot
    ax1 = plt.subplot(2, 9, (1, 4))
    format_panels(ax1)
    plt.title('a', loc='left', y=0.985, fontweight='bold')
    plt.title('(i)', loc='center', y=0.985, fontweight='bold')
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

    plt.ylabel(r'$u_d/u_0$', labelpad=0)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.03)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.05, fontsize=9)
    plt.xlabel('deleterious fitness effect', labelpad=1)
    plt.legend(loc='upper left', ncol=ncol)
    # make small ticks between selection coefs
    for xx in xi[:-1]:
        plt.axvline(xx + 0.5, 0, 0.2, lw=1, ls='--', color='k')

    ### 1b. CS DFE plot
    xi = np.arange(len(dfe_list[0]))
    n = sum(len(alf) for alf in alf_list)
    w = 0.8 / n
    s = -w * n / 2.0
    ncol = 2
    # plt.subplot(2, 8, (9, 11))
    ax1 = plt.subplot(2, 9, (10, 13))
    format_panels(ax1)
    ymax = 0
    offset = len(pmf_list)
    j = 0
    for (i, df) in enumerate(alf_list):
        if len(df) == 1:
            lbl = llist[i+offset]
            dd = np.concatenate(([0, 0], df[0]*100))
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
                dd = np.concatenate(([0, 0], dd * 100))
                # dd *= 100
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

    plt.ylabel(r'$\alpha$ (%)', labelpad=1.5)
    plt.ylim(0, 1.5 * ymax)
    plt.yticks([0, 1, 2], x=0.03)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.05, fontsize=9)
    plt.xlabel('adaptive fitness effect', labelpad=1)
    plt.legend(loc='upper left', ncol=ncol, columnspacing=0.1,
               handlelength=0.8, handletextpad=0.4)
    # make small ticks between selection coefs
    for xx in xi[:-1]:
        plt.axvline(xx + 0.5, 0, 0.2, lw=1, ls='--', color='k')

    ### 2. plot udel
    ax2 = plt.subplot(2, 9, 5)
    format_panels(ax2)
    plt.title('(ii)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'total $u_d/u_0$', labelpad=1)
    plt.yticks([0, 0.2, 0.4, 0.6], x=0.3)
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
    # plt.yticks(x=0.3)

    # plot alpha
    # plt.subplot(2,8,12)
    ax2 = plt.subplot(2, 9, 14)
    format_panels(ax2)
    plt.ylabel(r'total $\alpha$ (%)', labelpad=1)
    plt.yticks([0, 1, 2], x=0.3)
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
                alf = sum(dd)*100
                if alf < 0.000001:
                    continue
                if use_hatch:
                    plt.bar(0 + s, alf, w, color=clist[i + offset],
                            hatch=hatches[j-1]*5)
                    j += 1

                else:
                    plt.bar(0 + s, alf, w, color=clist[i+offset])
                s += w

    plt.xticks([])
    # plt.yticks(x=0.3)

    ### 3. plot pi0
    ax3 = plt.subplot(1, 9, 6)
    format_panels(ax3)
    plt.title('(iii)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=1)
    plt.yticks(x=0.3)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    ### 4. plot CLH
    ax4 = plt.subplot(1, 9, 7)
    format_panels(ax4)
    plt.title('(iv)', loc='center', y=0.985, fontweight='bold')
    plt.ylabel(r'$\mathrm{-\Delta}$' + 'CLL per site ' + r'$\times10^5$',
               labelpad=0)
    plt.yticks(x=0.3)
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

    ### 5. R^2 PLOT
    ax5 = plt.subplot(1,9, (8,9))
    format_panels(ax5)
    plt.title('b', loc='left', y=0.985, fontweight='bold')
    plt.ylabel(r'variance explained $(R^2)$', labelpad=1)
    plt.yticks(np.arange(0, 0.61, 0.1), x=0.08)
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    ws = ws[wmsk]
    for (i, rs) in enumerate(rsq_list):
        # save the 1Mb rsq vals
        rsq = rs[:,1][wmsk]
        plt.plot(ws, rsq, color=clist[i], marker='o',
                 linewidth=0, alpha=0.75, label=llist[i], ms=4)
    plt.xlabel('window size (bp)', labelpad=1)
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in xtck], y=0.035)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    plt.ylim(0.05, 0.65)

    # save the plot
    save_lbl = sname + '.parameters'
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)
    fsave = sdir + '/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=256)
    plt.close()


#%% MCVICKER COMPARISON PLOTS
def collate_and_rsq_mcvicker(fldr):
    flist = [fldr, 'mcvicker']
    llist = ['our map', 'McVicker']
    clist = ['darkorange', 'mediumpurple']
    sfldr = 'sfigs'

    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/' + sfldr
    if not os.path.isdir(sdir):
        os.mkdir(sdir)

    obs_flag = False
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.07, bottom=0.18, right=0.995, top=0.995,
                        wspace=0.3)

    # COLLATED PLOT
    ax1 = plt.subplot(1,3,(1,2))
    format_panels(ax1)
    for (i, rdir) in enumerate(flist):
        f_ape = root_dir + '/result/final_files/{}/nonsyn.collate.5.00e-03width.npy'.format(
            rdir)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        obs = pi / cnts
        prd = pr / cnts
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        if not obs_flag:
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=2)
            obs_flag = True
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.5)

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    plt.ylim(0.68, 1.15)

    # R^2 PLOT
    f_mcv = root_dir + '/result/final_files/mcvicker/rsq.log'
    f_cad = root_dir + '/result/final_files/{}/rsq.log'.format(fldr)
    w, mcv = np.loadtxt(f_mcv)[:16].T
    ape = np.loadtxt(f_cad)[:16, 1]

    w = np.log10(w)
    # plt.figure(figsize=(3.25, 2.5))
    ax2 = plt.subplot(133)
    format_panels(ax2)
    plt.plot(w, ape, label='our map', marker='o', lw=0,
             color='darkorange', ms=5)
    plt.plot(w, mcv, label='McVicker', marker='o', lw=0, color='mediumpurple',
             ms=5)
    # plt.plot(w, ex, label='our method (exon only)', marker='o', lw=0,
    #          color='fuchsia')

    plt.xlabel('window size (bp)', labelpad=1)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=1)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.04)
    plt.ylim(0.01, 0.65)
    # plt.legend(prop=dict(size=9), loc='upper left')
    plt.text(0.073, 0.92, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.742, 0.92, 'c', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)
    # f_save = root_dir + '/result/final_files/sfigs/fig_20.rsq.mcvicker.png'
    # plt.savefig(f_save, dpi=256)
    # plt.close()

    # plt.text(0.01, 0.93, fletters[0], transform=plt.gcf().transFigure)
    # plt.text(0.68, 0.93, fletters[1], transform=plt.gcf().transFigure)

    f_save = sdir + '/{}.collated.rsq.mcvicker.png'.format(fldr)
    plt.savefig(f_save, dpi=256)
    plt.close()


# PREDICTED VS. OBSERVED MCVICKER
def sorted_mcvicker_plot(fldr):
    """large formated observed vs. predicted plot comparing to mcvicker"""
    flist = [fldr, 'mcvicker']
    llist = ['our map', 'McVicker']
    clist = ['darkorange', 'mediumpurple']
    sfldr = 'sfigs'
    sdir = root_dir + '/result/final_files/' + sfldr

    plt.figure(figsize=(3.25,3.25))
    ax = plt.subplot(111)
    format_panels(ax)
    plt.subplots_adjust(left=0.15, top=0.995, right=0.985)
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(0.5, 1.05, 'without linked selection', ha='center',
             va='center')
    plt.plot([0, 1], [0, 1], color='k', ls='--')
    for (i, rdir) in enumerate(flist):
        span = 0.1
        fld = flist[i]
        color = clist[i]
        lbl = llist[i]
        # get 100 points data
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fld, 100)
        div, pi, pred = np.loadtxt(sort_file).T

        rst = cst_from_fldr(fld)
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        pi /= pi0
        pred /= pi0

        # get loess line and original points
        prlo, lolist = get_loess_line(fld, span)
        pilo = lolist[0]

        # plot horizontal line at y=1
        plt.axhline(y=1, color='k', alpha=0.8, ls='-')

        # plot predicted vs. observed
        plt.plot(pred, pi, marker='o', ms=3, color='white', lw=0)
        plt.plot(pred, pi, marker='o', ms=3, color=color, lw=0, alpha=0.5)

        # plot LOESS line
        plt.plot(prlo, pilo, lw=2, color=color, alpha=0.8, label=lbl)

    plt.ylim(-0.01, 1.25)
    plt.xlim(-0.01, 1.02)
    xtick = np.arange(0., 1.01, 0.1)
    ytick = np.arange(0., 1.11, 0.1)

    plt.xticks(xtick, y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick, x=0.02)
    # solve y=x rotation
    adlen = 1.03
    oplen = adlen**2 / 1.26
    rot = np.arctan((oplen/adlen)) * (180.0 / np.pi)
    # rot = 45 * oplen/adlen
    print 'rotation = {}'.format(rot)
    plt.text(0.5, 0.38, r'$y=x$', rotation=rot, ha='left', va='bottom',
             color='k', alpha=1)
    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)

    plt.text(0.01, 0.95, 'd', transform=plt.gcf().transFigure,
             fontweight='bold')
    # plt.text(0.68, 0.93, fletters[1], transform=plt.gcf().transFigure)

    f_save = sdir + '/{}.sort.vs.mcvicker.png'.format(fldr)
    plt.savefig(f_save, dpi=256)
    plt.close()


#%% CADD U GAMETIC
def cadd_u_gamete(pcons, ugam):
    plt.figure(figsize=(2.16, 2.16))
    plt.subplots_adjust(left=0.17, right=0.995, top=0.995, bottom=0.145)
    ax3 = plt.subplot(111)
    format_panels(ax3)

    pc_cons = [100-p for p in pcons]
    plt.bar(pc_cons, ugam, color='firebrick')
    plt.ylim(0, 2.3)
    plt.yticks(x=0.04)
    plt.xticks(pc_cons, y=0.04)
    plt.ylabel(r'$u_d$ per gamete', labelpad=1)
    plt.xlabel('% CADD threshold', labelpad=1)
    plt.text(0.175, 0.92, 'G', transform=plt.gcf().transFigure)

    f_save = root_dir + '/result/final_files/sfigs/cadd_u_gamete.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


#%% GAMETIC U FOR SELECT ANNOTATIONS
adict = {'fishex': 27085413,
         'fishnex': 145872523,
         'fishall': 172957936,
         'caddnex': 142838349,
         'caddex': 30054071,
         'caddall': 172892420}
ufish = cst_from_fldr('fish_cons94_new').stat.utot[0]
ufishex, ufishnex = cst_from_fldr('fish_cons94_new_exnex').stat.utot
ucadd = cst_from_fldr('cadd94_gmask_v1.6_without_bstat').stat.utot[0]
ucaddex, ucaddnex = cst_from_fldr('cadd94_gmask_v1.6_without_bstat_exnex').stat.utot
print 'fish all: {}'.format(adict['fishall'] * ufish)
print 'cadd all: {}'.format(adict['caddall'] * ucadd)
print 'fish ex: {}'.format(adict['fishex'] * ufishex)
print 'cadd ex: {}'.format(adict['caddex'] * ucaddex)
print 'fish nex: {}'.format(adict['fishnex'] * ufishnex)
print 'cadd nex: {}'.format(adict['caddnex'] * ucaddnex)


#%% FIGURE S7 COMPARE THRESHOLDS
mnblist = ['cadd94_gmask_v1.6_without_bstat_bth_000',
           'cadd94_gmask_v1.6_without_bstat_minb_119',
           'cadd94_gmask_v1.6_without_bstat_minb_513',
           'cadd94_gmask_v1.6_without_bstat']
llist = ['none'] + [r'B = {:.2f}'.format(b) for b in 0.2, 0.5, 0.6]
clist = 'Maroon dodgerblue deeppink darkorange'.split()
sname = 'fig_S7c'
sfldr = 'sfigs'
single_bsanno(mnblist, llist, clist, sname, sfldr, letters=('c', 'd'))


#%% FIGURE S14 COMPARE PHYLO DEPTHS
spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()
flist = []
for sp in spec:
    flist.append('{}_cons94_new'.format(sp))
clist = 'purple red lightcoral darkturquoise' \
        ' firebrick fuchsia darkorange'.split()
llist = ['4-ape', '8-prim', '12-pros', '25-supr', '50-laur',
         '61-mamm', '99-vert']
sname = 'fig_S14'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)

xlab = 'phylogenetic depth'
llist = '4 8 12 25 50 61 99'.split()
# rsq_3_size(flist, llist, sname, xlab, rotation=45)
flist = ['ape_cons94_new', 'euarchontoglires_cons94_new',
         'fish_cons94_new']
clist = 'purple darkturquoise darkorange'.split()
llist = ['ape (4)', 'supr (25)', 'vert (99)']
# combined_chr1(flist, llist, clist, sname, 'd')
# collate_and_sort(flist, llist, clist, sname, fletters=('e', 'f'))


#%% FIGURE S16 COMPARE PHASTCONS CUTOFFS IN FISH
pcons = range(98, 90, -1)
flist = ['fish_cons{}_new'.format(p) for p in pcons]
clist = 'purple crimson lightcoral indianred darkorange' \
        ' firebrick darkred darkturquoise'.split()
llist = ['{}%'.format(100-pc) for pc in pcons]
sname = 'fig_S16'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
xlab = 'fraction conserved'
# rsq_3_size(flist, llist, sname, xlab, rotation=45)
pcons = [98, 94, 91]
flist = ['fish_cons{}_new'.format(p) for p in pcons]
clist = 'purple  darkorange darkturquoise'.split()
llist = ['{}%'.format(100-pc) for pc in pcons]
# combined_chr1(flist, llist, clist, sname, 'd')
# collate_and_sort(flist, llist, clist, sname, fletters=('e', 'f'))


#%% FIGURE S17 COMPARE CONS TO GENIC
flist = ['fish_cons94_new', 'genic_gmask']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'genic']
sllist = ['genic: splice', 'genic: CDS']
sname = 'fig_S17'
sfldr = 'sfigs'
ytck = [0, 0.5, 1, 1.5, 2]
multi_bsanno(flist, llist, sllist, clist, sname, sfldr, ytck)
# combined_chr1(flist, llist, clist, sname, 'c')
# collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))


#%% FIGURE S19 CONSERVED STRATIFIED EX/NEX
flist = ['fish_cons94_new', 'fish_cons94_new_exnex']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'conserved: exonic']
sllist = ['conserved: exonic', 'conserved: non-exonic']
sname = 'fig_S19'
sfldr = 'sfigs'
ytck = [0, 0.2, 0.4, 0.6, 0.8, 1]
multi_bsanno(flist, llist, sllist, clist, sname, sfldr, ytck)
# llist = ['conserved', 'conserved (exonic/non-exonic)']
# combined_chr1(flist, llist, clist, sname, 'c')
# collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))


#%% FIGURE S20 COMPARE CONS EX-CONS/NEX-CONS SEPARATED
flist = ['fish_cons94_new', 'fish_cons94_new_exonic',
         'fish_cons94_new_nonexonic']
clist = ['darkorange', 'dodgerblue', 'fuchsia']
llist = ['conserved', 'conserved: exonic', 'conserved: non-exonic']
sname = 'fig_S20'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
# combined_chr1(flist, llist, clist, sname, 'c')
# collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))


#%% FIGURE S22 COMPARE CONS TO ENCODE cCREs
flist = ['fish_cons94_new', 'ccre_split']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'coding+cCRE']
sllist = ['CDS', r'$\mathrm{ELS_1}$', r'$\mathrm{CTCF_1}$']
sname = 'fig_S22'
sfldr = 'sfigs'
ytck = [0, 0.5, 1, 1.5]

multi_bsanno(flist, llist, sllist, clist, sname, sfldr, ytck)
# combined_chr1(flist, llist, clist, sname, 'c')
# collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))


#%% FIGURE S24 COMPARE CONS AND CADD
flist = ['fish_cons94_new', 'cadd94_gmask_v1.6_without_bstat']
clist = ['darkorange', 'dodgerblue']
llist = ['conserved', 'CADD']
sname = 'fig_S24'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
# combined_chr1(flist, llist, clist, sname, 'c')
# collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))
# print_stats(flist)

#%% FIGURE S26 CADD COMPARE % SCORE THRESHOLDS
pcons = range(98, 90, -1)
flist = ['cadd{}_gmask_v1.6_without_bstat'.format(p) for p in pcons]
# for pc in pcons:
#     if pc == 94:
#         flist.append('cadd94_gmask_mnb_378')
#     else:
#         flist.append('cadd{}_gmask'.format(pc))
clist = 'purple crimson lightcoral indianred darkorange' \
        ' firebrick darkred darkturquoise'.split()
llist = ['{}%'.format(100-pc) for pc in pcons]
sname = 'fig_S26'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
# utots = [cst_from_fldr(fl).stat.utot[0] for fl in flist]
# ugams = [0.01*(100-p)*2.881e9*u for (u, p) in zip(utots, pcons)]
# cadd_u_gamete(pcons, ugams)
xlab = 'CADD threshold'
# rsq_3_size(flist, llist, sname, xlab, 35)
pcons = [98, 94, 91]
flist = ['cadd{}_gmask_v1.6_without_bstat'.format(p) for p in pcons]
clist = 'purple  darkorange darkturquoise'.split()
llist = ['{}%'.format(100-pc) for pc in pcons]
# combined_chr1(flist, llist, clist, sname, 'd')
# collate_and_sort(flist, llist, clist, sname, fletters=('e', 'f'))


#%% FIGURE S27 BS+CS PLOT
flist = ['fish_cons94_new', 'YRI_nonsyn_s1',
         'fish_cons91_gmask_hc_subs', 'fish_cons98_gmask_hc_subs']
clist = ['darkorange', 'royalblue', 'mediumpurple', 'firebrick']
llist = ['conserved', 'NS', '9%', '2%']
sllist = ['NS', 'other']*2
sname = 'fig_S27'
sfldr = 'sfigs'
bscs_joint(flist, llist, sllist, clist, sname, sfldr)
# combined_chr1(flist, llist, clist, sname, 'c')
# collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))


#%% FIGURE S29 COMPARISON TO MCVICKER
flist = ['cadd94_gmask_v1.6_without_bstat', 'mcvicker']
clist = ['darkorange', 'mediumpurple']
llist = ['our map', 'McVicker']
sname = 'fig_S29'
sfldr = 'sfigs'
# combined_chr1(flist, llist, clist, sname, fletter='a')
fldr = 'cadd94_gmask_v1.6_without_bstat'
collate_and_rsq_mcvicker(fldr)
# sorted_mcvicker_plot(fldr)


#%% FOR SECTION 8: COMPARE STANDARD CADD94 W/ CADD94 FILTER CG
flist = ['cadd94_gmask_mnb_378', 'cadd94_gmask_filter_CG']
clist = ['darkorange', 'dodgerblue']
llist = ['CADD', 'CADD (filter CG)']
sname = 'cadd_vs_cadd_filter_cg'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
combined_chr1(flist, llist, clist, sname, 'c')
collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))
#%% COLLATED CADD VS. BEST CS MODEL
flist = ['cadd94_gmask_mnb_378', 'YRI_nonsyn_s1']
llist = ['background selection', 'selective sweeps']
clist = ['darkorange', 'darkturquoise']
sname = 'cadd_vs_sweeps'
sfldr = 'sfigs'
# single_bsanno(flist, llist, clist, sname, sfldr)
# combined_chr1(flist, llist, clist, sname, 'c')
collate_and_sort(flist, llist, clist, sname, fletters=('', ''), legend_on=1)

#%% COMPARE NEW GENETIC MAPS
# flist = ['cadd94_gmask_mnb_378', 'cadd94_gmask_YRI_LD',
#          'cadd94_gmask_deCODE_2019']
# clist = ['darkorange', 'dodgerblue', 'darkturquoise']
# llist = ['CADD AA', 'CADD YRI LD', 'CADD deCODE']
# sname = 'CADD_compare_gmaps'
flist = ['fish_cons94_gmask_mnb_378', 'fish_cons94_gmask_YRI_LD',
         'fish_cons94_gmask_deCODE_2019']
clist = ['darkorange', 'dodgerblue', 'darkturquoise']
llist = ['cons AA', 'cons YRI LD', 'cons deCODE']
sname = 'fish_compare_gmaps'
sfldr = 'sfigs'
single_bsanno(flist, llist, clist, sname, sfldr)
combined_chr1(flist, llist, clist, sname, 'c')
collate_and_sort(flist, llist, clist, sname, fletters=('d', 'e'))
#%% COMPARE DIFFERENT CADD SCORES
spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()
flist = ['cadd94_gmask_mnb_378', 'cadd94_gmask_v1.6',
         'cadd94_gmask_v1.6_without_bstat']
clist = ['darkorange', 'dodgerblue', 'purple']
llist = ['CADD v1.4 (used in main)', 'CADD v1.6', 'CADD v1.6 without B']
sname = 'compare_CADD_scores'
sfldr = 'sfigs'
# single_bsanno(flist, llist, clist, sname, sfldr)
xlab = 'CADD score version'
llist2 = ['1.4', '1.6', '1.6 noB']
rsq_3_size(flist, llist2, sname, xlab, rotation=0)
# combined_chr1(flist, llist, clist, sname)
# collate_and_sort(flist, llist, clist, sname, fletters=('E', 'F'))
#%% FIGURE A48: COMPARE DIFFERENT LEAVE ONE OUT TO ALL DATA
flist = ['cadd94_gmask_v1.6_without_bstat', 'cadd94_gmask_v1.6_without_bstat_jackknife_results',
         'dropped_chrom_results']
llist = ['all data', 'leave 2 Mb out', 'leave autosome out']
clist = ['darkorange', 'deepskyblue', 'fuchsia']
sname = 'compare_LOO'

# rsq_3_size(flist, llist2, sname, xlab, rotation=0)
# combined_chr1(flist, llist, clist, sname)
collate_and_sort(flist, llist, clist, sname, fletters=('C', 'D'))

#%%