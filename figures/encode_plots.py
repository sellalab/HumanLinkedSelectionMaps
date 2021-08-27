__author__ = 'davidmurphy'


import os
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir
from figures.other_code.summary_slide import cst_from_fldr
from figures.common_functions import format_panels
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


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
    plt.subplots_adjust(left=0.065, wspace=1.75, right=0.98, bottom=0.15,
                        top=0.91)
    # DFE plot
    ax1 = plt.subplot(1, 9, (1, 4))
    format_panels(ax1)
    plt.title('A (i)', loc='center', y=0.985)
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
                # add a zero pad:
                if len(dd) < 6:
                    zpad = [0] * (6-len(dd))
                    dd = np.concatenate((zpad, dd))
                lbl = llist[i] + ' ' + subllist[j]
                ymax = max(dd.max(), ymax)
                if sum(dd) < 0.01:
                    continue
                plt.bar(xi + s, dd, w, label=lbl, color=clist[i],
                        hatch=hatches[j]*5, align='edge')
                s += w
                j += 1

    plt.ylabel(r'mutation rate (units of $\mu_0$)', labelpad=1.5)
    plt.ylim(0, 7)
    # plt.ylim(0, 1.5*ymax)
    plt.yticks(ytck, x=0.03)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.035, fontsize=9)
    plt.xlabel('deleterious fitness effect', labelpad=1)
    plt.legend(loc='upper left', ncol=ncol)

    # plot udel
    ax2 = plt.subplot(1, 9, 5)
    format_panels(ax2)
    plt.title('(ii)', loc='center', y=0.985)
    plt.ylabel(r'mutation rate (units of $\mu_0$)', labelpad=1)
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
    plt.ylim(0, 7)

    # plot pi0
    ax3 = plt.subplot(1,9,6)
    format_panels(ax3)
    plt.title('(iii)', loc='center', y=0.985)
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=1)
    plt.yticks(x=0.3)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    # plot CLH
    ax4 = plt.subplot(1, 9, 7)
    format_panels(ax4)
    plt.title('(iv)', loc='center', y=0.985)
    plt.ylabel(r'$- \Delta CLLH\ (x10^5)$', labelpad=3)
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

    # R^2 PLOT
    ax5 = plt.subplot(1,9, (8,9))
    format_panels(ax5)
    plt.title('B', loc='center', y=0.985)
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
    # plt.legend(prop=dict(size=9), loc='upper left')

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


#%% COMPARE CONS TO CHROMHMM HUVEC
flist = ['cadd94_gmask_mnb_378', 'genic_plus_ccre']
clist = ['darkorange', 'dodgerblue']
llist = ['CADD', 'coding/cCRE']
sllist = ['coding', 'ELS', 'CTCF']
sname = 'cCRE'
sfldr = 'sfigs'
# ytck = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
ytck = np.arange(8)
multi_bsanno(flist, llist, sllist, clist, sname, sfldr, ytck)


#%%
# --bs_annos=('cds', 'cCRE_PLS_filtered', 'cCRE_ELS_filtered', 'cCRE_CTCF_filtered')
n_els = 2.0816e+08
n_pls = 8.6598e+06
n_ctcf = 1.3568e+07
n_cds = 30245587.0

nlist = [n_cds, n_pls, n_els, n_ctcf]
flist = [n/2.88e9 for n in nlist]

ulist = [3.929897151899468e-09, 9.06532437039484e-28,
         5.046678154823988e-09, 9.360494448885902e-08]

bigulist = [n*u for (n, u) in zip(nlist, ulist)]
frac_bigu = [u/sum(bigulist) for u in bigulist]

print 'fraction of sites: cds={:.3f} PLS={:.3f} ELS={:.3f} CTCF={:.3f}'.format(*flist)
print 'U per annotation: cds={:.3f} PLS={:.3f} ELS={:.3f} CTCF={:.3f}'.format(*bigulist)
print 'fraction of total U per mutation: cds={:.3f} PLS={:.3f} ELS={:.3f} CTCF={:.3f}'.format(*frac_bigu)


#%%