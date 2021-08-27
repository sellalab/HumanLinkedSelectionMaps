__author__ = 'davidmurphy'


import os
import gzip
from classes.runstruct import root_dir, np, human_autosomes
mcv_dir = root_dir + '/lsm_matlab/data/mcvick/bkgd_data'
hg18 = dict(chr1=247249719, chr2=242951149, chr3=199501827,
                chr4=191273063, chr5=180857866, chr6=170899992,
                chr7=158821424, chr8=146274826, chr9=140273252,
                chr10=135374737, chr11=134452384, chr12=132349534,
                chr13=114142980, chr14=106368585, chr15=100338915,
                chr16=88827254, chr17=78774742, chr18=76117153,
                chr19=63811651, chr20=62435964, chr21=46944323,
                chr22=49691432, chrX=154913754, chrY=57772954)
hg18_aut = sum(hg18[ch] for ch in human_autosomes)


#%%

def get_exon(chrom):
    f_path = mcv_dir + '/exon/features/{}.coords.gz'.format(chrom)
    arr_1 = np.zeros(shape=hg18[chrom], dtype=bool)
    with gzip.open(f_path, 'r') as f:
        for line in f:
            start, end = map(int, line.split()[1:])
            arr_1[start-1:end] = True

    f_path = mcv_dir + '/exon/features/{}_random.coords.gz'.format(chrom)
    arr_2 = np.zeros(shape=hg18[chrom], dtype=bool)
    if os.path.isfile(f_path):
        with gzip.open(f_path, 'r') as f:
            for line in f:
                start, end = map(int, line.split()[1:])
                arr_2[start-1:end] = True

    return arr_1, arr_2


def get_gcon_ex(chrom):
    f_path = mcv_dir + '/gcons_ex/features/{}.coords.gz'.format(chrom)
    arr_1 = np.zeros(shape=hg18[chrom], dtype=bool)
    with gzip.open(f_path, 'r') as f:
        for line in f:
            start, end = map(int, line.split()[1:])
            arr_1[start-1:end] = True

    f_path = mcv_dir + '/gcons_ex/features/{}_random.coords.gz'.format(chrom)
    arr_2 = np.zeros(shape=hg18[chrom], dtype=bool)
    if os.path.isfile(f_path):
        with gzip.open(f_path, 'r') as f:
            for line in f:
                start, end = map(int, line.split()[1:])
                arr_2[start-1:end] = True

    return arr_1, arr_2


def get_gcon_nex(chrom):
    f_path = mcv_dir + '/gcons_nex/features/{}.coords.gz'.format(chrom)
    arr = np.zeros(shape=hg18[chrom], dtype=bool)
    with gzip.open(f_path, 'r') as f:
        for line in f:
            start, end = map(int, line.split()[1:])
            arr[start-1:end] = True
    sum_1 = arr.sum()

    f_path = mcv_dir + '/gcons_nex/features/{}_random.coords.gz'.format(chrom)
    if os.path.isfile(f_path):
        arr = np.zeros(shape=hg18[chrom], dtype=bool)
        with gzip.open(f_path, 'r') as f:
            for line in f:
                start, end = map(int, line.split()[1:])
                arr[start-1:end] = True
        sum_2 = arr.sum()
    else:
        sum_2 = 0

    return sum_1+sum_2


def get_intersect(chrom):
    x1, x2 = get_exon(chrom)
    cx1, cx2 = get_gcon_ex(chrom)
    return np.sum(x1&cx1)+np.sum(x2&cx2), np.sum(cx1)+np.sum(cx2)

#%%
tot_int, tot_gcon = 0.0, 0.0
tot_gcon_nex = 0.0
for ch in human_autosomes:
    n, m = get_intersect(ch)
    tot_int += n
    tot_gcon += m
    tot_gcon_nex += get_gcon_nex(ch)

p_gnex = tot_gcon_nex / hg18_aut
p_int, p_gcon = tot_int/hg18_aut, tot_gcon/hg18_aut
print '{} {} {}'.format(p_int, p_gcon, p_gnex)
#%%