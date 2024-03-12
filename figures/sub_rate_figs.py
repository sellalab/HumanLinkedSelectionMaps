__author__ = 'davidmurphy'


from classes.runstruct import root_dir, human_autosomes, np
from phast.neutcons_subrates import parse_exptotsub
from collections import defaultdict
from itertools import product
pdir = root_dir + '/data/phast/'
import matplotlib.pyplot as plt
import seaborn
from classes.knowngene import complement_strand, complement_base


#%%
def neutcons_uratio_context(spec, mod='U3S', n_cell=64, nonbgc=False):
    # get dict of conserved rates
    # mod = 'U3S'
    # spec = 'primate'

    # set input file formats
    cfmt = pdir + 'hcgo_cons/{ch}.hcgo.{sp}.cons.95.{md}.exptotsub'
    nfmt = pdir + 'hcgo_neutral/{ch}.hcgo.neut.{md}.exptotsub'

    # initialize empty arrays for rates
    cmat = np.zeros(shape=(n_cell, n_cell))
    nmat = np.zeros(shape=(n_cell, n_cell))

    # combine data across chromosomes
    for chrom in human_autosomes:
        # if chrom == 'chr2' or chrom == 'chr7':
        #     continue
        f_cxts = cfmt.format(ch=chrom, sp=spec, md=mod)
        f_nxts = nfmt.format(ch=chrom, md=mod)
        cmat += parse_exptotsub(f_cxts, n_cell)['hg19']
        nmat += parse_exptotsub(f_nxts, n_cell)['hg19']

    # convert sub counts to rates
    n_contexts, c_contexts = [], []
    for i in range(n_cell):
        c_sum = np.sum(cmat[i, :])
        n_sum = np.sum(nmat[i, :])
        c_contexts.append(c_sum)
        n_contexts.append(n_sum)
        cmat[i, :] /= c_sum
        nmat[i, :] /= n_sum

    n_contexts = [1.0 * n / sum(n_contexts) for n in n_contexts]
    c_contexts = [1.0 * c / sum(c_contexts) for c in c_contexts]

    # create trip pair dicts for neutral and cons rates
    ccd = dict()
    nnd = dict()
    cdict = defaultdict(dict)
    ndict = defaultdict(dict)
    triplets = list(product('ACGT', repeat=3))
    for i in range(n_cell):
        for j in range(n_cell):
            trip_1 = ''.join(triplets[i])
            trip_2 = ''.join(triplets[j])
            cdict[trip_1][trip_2] = cmat[i,j]
            ndict[trip_1][trip_2] = nmat[i,j]
            ccd[trip_1] = c_contexts[i]
            nnd[trip_1] = n_contexts[i]

    # get ratios for rates for single substitution triplet pairs
    ratios = []
    seen_trips = []
    s_type = []
    c_pct = []
    n = 0
    for trip in triplets:
        trip = ''.join(trip)
        if trip in seen_trips:
            continue
        seen_trips.append(trip)
        # for i in range(3):
        bases = list(trip)
        for b in 'ACGT':
            if (not nonbgc) and (bases[1] != b):
                new_bases = list(trip)
                new_bases[1] = b
                new_trip = ''.join(new_bases)

                # get reverse complement of each as well
                rtrip = complement_strand(trip)
                # put complement in the seen triplets list so its not doubled
                seen_trips.append(rtrip)
                rnew_trip = complement_strand(new_trip)
                crate = cdict[trip][new_trip] + cdict[rtrip][rnew_trip]
                nrate = ndict[trip][new_trip] + ndict[rtrip][rnew_trip]

                # crate = cdict[trip][new_trip] #* ccd[trip]
                # nrate = ndict[trip][new_trip] #* nnd[trip]
                # if spec == 'fish':
                #     print '{} {}/{} {}/{} {}'.format(n, trip, rtrip, new_trip,
                #                                      rnew_trip, crate / nrate)
                # print n, trip, new_trip, crate/nrate
                n += 1
                ratios.append(crate / nrate)
                s_type.append('{}:{}'.format(trip, new_trip))
                c_pct.append(ccd[trip])

            if nonbgc and (bases[1] == complement_base(b)):
                new_bases = list(trip)
                new_bases[1] = b
                new_trip = ''.join(new_bases)

                # get reverse complement of each as well
                rtrip = complement_strand(trip)
                # put complement in the seen triplets list so its not doubled
                seen_trips.append(rtrip)
                rnew_trip = complement_strand(new_trip)
                crate = cdict[trip][new_trip] + cdict[rtrip][rnew_trip]
                nrate = ndict[trip][new_trip] + ndict[rtrip][rnew_trip]

                # crate = cdict[trip][new_trip] #* ccd[trip]
                # nrate = ndict[trip][new_trip] #* nnd[trip]
                # if spec == 'fish':
                #     print '{} {}/{} {}/{} {}'.format(n, trip, rtrip, new_trip,
                #                                      rnew_trip, crate/nrate)
                # print n, trip, new_trip, crate/nrate
                n += 1
                ratios.append(crate/nrate)
                s_type.append('{}:{}'.format(trip, new_trip))
                c_pct.append(ccd[trip])

    return np.array(ratios), s_type, c_pct


def neutcons_uratio_single(spec, mod='UNREST', n_cell=4, nonbgc=False, pc=95):
    # get dict of conserved rates
    # mod = 'U3S'
    # spec = 'primate'

    # set input file formats
    cfmt = pdir + 'hcgo_cons/{ch}.hcgo.{sp}.cons.{pc}.{md}.exptotsub'
    nfmt = pdir + 'hcgo_neutral/{ch}.hcgo.neut.{md}.exptotsub'

    # initialize empty arrays for rates
    cmat = np.zeros(shape=(n_cell, n_cell), dtype=float)
    nmat = np.zeros(shape=(n_cell, n_cell), dtype=float)

    # combine data across chromosomes
    for chrom in human_autosomes:
        # if chrom == 'chr2' or chrom == 'chr7':
        #     continue
        f_cxts = cfmt.format(ch=chrom, sp=spec, md=mod, pc=pc)
        f_nxts = nfmt.format(ch=chrom, md=mod)
        print(parse_exptotsub(f_cxts, n_cell)['hg19'])
        cmat += parse_exptotsub(f_cxts, n_cell)['hg19']
        nmat += parse_exptotsub(f_nxts, n_cell)['hg19']

    # convert sub counts to rates
    n_contexts, c_contexts = [], []
    for i in range(n_cell):
        c_sum = np.sum(cmat[i, :])
        n_sum = np.sum(nmat[i, :])
        c_contexts.append(c_sum)
        n_contexts.append(n_sum)
        cmat[i, :] /= c_sum
        nmat[i, :] /= n_sum

    n_contexts = [1.0 * n / sum(n_contexts) for n in n_contexts]
    c_contexts = [1.0 * c / sum(c_contexts) for c in c_contexts]

    # create trip pair dicts for neutral and cons rates
    ccd = dict()
    nnd = dict()
    cdict = defaultdict(dict)
    ndict = defaultdict(dict)
    bases = 'ACGT'
    for i in range(n_cell):
        for j in range(n_cell):
            b_1 = bases[i]
            b_2 = bases[j]
            cdict[b_1][b_2] = cmat[i,j]
            ndict[b_1][b_2] = nmat[i,j]
            ccd[b_1] = c_contexts[i]
            nnd[b_1] = n_contexts[i]

    # get ratios for rates for single substitution triplet pairs
    ratios = []
    s_type = []
    c_pct = []
    n = 0
    for b_1 in bases:
        for b_2 in bases:
            if b_1 != b_2:
                crate = cdict[b_1][b_2]
                nrate = ndict[b_1][b_2]
                # print n, b_1, b_2, crate / nrate
                n += 1
                ratios.append(crate / nrate)
                s_type.append('{}:{}'.format(b_1, b_2))
                c_pct.append(ccd[b_1])
                if nonbgc and b_1 == complement_base(b_2):
                    crate = cdict[b_1][b_2]
                    nrate = ndict[b_1][b_2]
                    # print n, b_1, b_2, crate/nrate
                    n += 1
                    ratios.append(crate / nrate)
                    s_type.append('{}:{}'.format(b_1, b_2))
                    c_pct.append(ccd[b_1])

    return np.array(ratios), s_type, c_pct


#%%

cols = ['salmon', 'indianred', 'firebrick', 'darkturquoise', 'teal',
        'goldenrod', 'royalblue', 'darkorange']
spec = ['ape', 'primate', 'prosimian', 'euarchontoglires', 'laurasiatheria',
        'mammal','fish','cadd']


#%%
def combined_ratios():
    all_r = []
    non_bgc = []
    s_r = []
    for i in range(7):
        sp = spec[i]
        r, _, cpct = neutcons_uratio_context(sp)
        nr, _, ncpct = neutcons_uratio_context(sp, nonbgc=True)
        sr, _, scpct = neutcons_uratio_single(sp)
        # all_r.append(np.average(r, weights=cpct))
        # non_bgc.append(np.average(nr, weights=ncpct))
        # s_r.append(np.average(sr, weights=scpct))
        s_r.append(np.average(sr))
        all_r.append(np.average(r))
        non_bgc.append(np.average(nr))
        print(sp, sr.mean(), r.mean(), nr.mean())

    xv = np.arange(7)
    plt.figure(figsize=(12,4))

    # plot all subs bar
    plt.bar(xv-0.375, all_r, width=0.25, label='all substitutions')
    for i in range(7):
        xi = xv[i]-0.375
        yi = all_r[i]
        plt.text(xi, yi*0.5, '{:.3f}'.format(yi), rotation=90, ha='center',
                 color='white')

    # plot nonBGC subs bar
    plt.bar(xv-0.125, non_bgc, width=0.25, label='non-BGC substitutions')
    for i in range(7):
        xi = xv[i]-0.125
        yi = non_bgc[i]
        plt.text(xi, yi*0.5, '{:.3f}'.format(yi), rotation=90, ha='center',
                 color='white')

    # plot singleton subs bar
    plt.bar(xv+0.125, s_r, width=0.25, label='singleton model')
    for i in range(7):
        xi = xv[i]+0.125
        yi = s_r[i]
        plt.text(xi, yi*0.5, '{:.3f}'.format(yi), rotation=90, ha='center',
                 color='white')

    plt.ylabel('conserved/neutral ratio')
    plt.xlabel('phylogeny')
    plt.xticks(xv-0.125, spec)
    plt.legend(ncol=2)
    plt.ylim(0,0.71)
    # plt.show()
    # fig_save = pdir + 'BGC_nonBGC_singleton_ratios_weighted.png'
    fig_save = pdir + 'BGC_nonBGC_singleton_ratios.png'
    plt.savefig(fig_save, dpi=256)
    plt.close()


#%%
def single_utypes():
    plt.figure(figsize=(12,4))
    plt.subplots_adjust(bottom=0.2)
    nr, stypes, cpct = neutcons_uratio_single('fish', nonbgc=True)
    # nr, stypes, cpct = neutcons_uratio_context('fish', nonbgc=True)
    plt.plot(nr, label='fish', color=cols[-1])
    plt.ylabel('conserved/neutral ratio')
    plt.xticks(np.arange(len(nr)), stypes, rotation=90)
    plt.xlabel('substitution')
    plt.legend()
    # fig_save = pdir + 'fish_nonBGC_triplet_ratios.png'
    fig_save = pdir + 'fish_nonBGC_single_ratios.png'
    plt.savefig(fig_save, dpi=256)
    plt.close()


#%%
def unrest_ratios():
    s_r = []
    for i in range(7):
        sp = spec[i]
        sr, _, scpct = neutcons_uratio_single(sp)
        s_r.append(np.average(sr))
        print(sp, sr.mean())

    xv = np.arange(7)
    plt.figure(figsize=(12,4))

    # plot singleton subs bar
    plt.bar(xv, s_r, width=0.8, label='singleton model', color='firebrick')
    for i in range(7):
        xi = xv[i]
        yi = s_r[i]
        plt.text(xi, yi*0.5, '{:.3f}'.format(yi), rotation=90, ha='center',
                 color='white', fontsize=14)

    plt.ylabel('conserved/neutral ratio')
    plt.xlabel('phylogeny')
    plt.xticks(xv, spec)
    plt.legend(ncol=2)
    plt.ylim(0,0.71)
    # plt.show()
    # fig_save = pdir + 'BGC_nonBGC_singleton_ratios_weighted.png'
    fig_save = pdir + 'singleton_ratios.png'
    plt.savefig(fig_save, dpi=256)
    plt.close()


#%%
def unrest_ratios_all_pct(plot_type=None):
    uhi = 1.51
    shift = -0.4
    width = 0.8 / 6
    xv = np.arange(8)
    plt.figure(figsize=(10, 3.3))
    plt.subplots_adjust(bottom=0.15, top=0.98, right=1, left=0.07)
    percents = [91, 93, 94, 95, 97, 99]
    for (j, pct) in enumerate(percents):
        s_r = []
        for i in range(8):
            sp = spec[i]
            sr, _, scpct = neutcons_uratio_single(sp, pc=pct)
            s_r.append(np.average(sr))
            print(pct, sp, sr.mean())

        # plot singleton subs bar
        lbl = '{}%'.format(100-pct)
        if plot_type == 'udel':
            ud = uhi * (1.0 - np.array(s_r))
            plt.bar(xv+shift, ud, width=width, label=lbl, color=cols[j])
        else:
            plt.bar(xv+shift, s_r, width=width, label=lbl, color=cols[j])

        for i in range(8):
            xi = xv[i]+shift
            if plot_type == 'udel':
                ytext = 0.25
                yi = ud[i]
            else:
                ytext = 0.03
                yi = s_r[i]
            st = '{:.3f}'.format(yi)
            plt.text(xi, ytext, st, rotation=90, ha='center', va='bottom',
                     color='white', fontsize=10)
        shift += width

    # plt.ylabel('conserved/neutral ratio')
    if plot_type == 'udel':
        loc = 'upper left'
        plt.ylabel(r'$\mu_{del}$' + ' estimate')
    else:
        loc = 'upper right'
        plt.ylabel(r'$S_{cons}/S_{neut}$')

    plt.xlabel('annotation')
    xlab = [s.upper() if s == 'cadd' else s for s in spec]
    plt.xticks(xv-(0.5*width), xlab)

    plt.legend(ncol=2, loc=loc, title='threshold', frameon=1,
               framealpha=1, facecolor='white')
    # plt.ylim(0, 0.69)
    if plot_type == 'udel':
        fig_save = pdir + 'singleton_estimates_across_pct_cons_plus_CADD.png'
    else:
        fig_save = pdir + 'singleton_ratios_across_pct_cons_plus_CADD.png'

    plt.savefig(fig_save, dpi=256)
    plt.close()


def main():
    # unrest_ratios_all_pct()
    unrest_ratios_all_pct('udel')


if __name__ == '__main__':
    main()
