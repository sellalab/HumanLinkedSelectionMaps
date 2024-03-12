__author__ = 'davidmurphy'

from classes.runstruct import RunStruct, np, root_dir

#%%
def make_initfiles():
    ns = 8
    wn = 7000
    # for wn in 4000, 6000, 7000, 8000, 9000, 11000:
    for sl in [10, 5, 2]:
        pct = int(100.0 * 1.0 / sl)
        tk = '{}spc.{}.{}pct'.format(ns, wn, pct)
        # if cp:
        #     tk += '.filter.CpG'
        a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk,
                      bs_annos=('primate_cons94_Segments_rm_neut',),
                      bdir='primate_cons94_Segments_rm_neut',
                      bdfe=(np.power(10.0, np.arange(-3.5, -1.5, 0.5)),))
        a.save()

#%%
wdir = root_dir + '/data/cons/wigz'
fpro = wdir + '/prosimian/prosimian.scores.dist.txt'
fpri = wdir + '/primate/primate.scores.dist.txt'


def dist_percentile(fdist, pct):
    """return cutoff value for distribution of scores at given percentile"""
    # load distribution of scores from the file
    p, n = np.loadtxt(fdist).T
    # rescale scores as ints from 0-1000
    p = (p*1000).astype('u2')
    # rescale counts to cumulative fraction
    n /= n.sum()
    n = np.cumsum(n)
    # get the index where cumulative fraction is >= pct
    i = np.searchsorted(n, pct)

    # return the p value at the select percentile
    return p[i]
#%%

print dist_percentile(fpro, 0.95)
print dist_percentile(fpri, 0.95)

#%%