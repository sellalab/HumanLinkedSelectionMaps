__author__ = 'davidmurphy'


import os
from sys import argv, stdin


# [n=anno, ti=2.0, tj=4.5, ci=1, cj=22, pi=-1, pj=-1, plen=1e7, hr=24, fx=fexe, pts=0/1, mem=mem, ivl=ivl]
cfmt = 'sh/bkgd_batch.sh {fi} n={an} ti={tx} tj={tx} ci={cx} cj={cx} pi={px} pj={px} plen={pl} hr=48 mem=2G >> {an}.bjobs.log'
cfmt2 = 'sh/bkgd_batch.sh {fi} n={an} ti={tx} tj={tx} ci={cx} cj={cx} pi=0 plen={pl} hr=48 mem=2G >> {an}.bjobs.log'

ffmt = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files/YRI.{an}.BS1.6.CS0.0.NOT_STARTED.initial.txt'
# chr19 cadd95_gmask 00003162 p0 1e07
t_grid = [4.5, 4, 3.5, 3, 2.5, 2]
rev_t = dict(('{:.8f}'.format(10**-t)[2:], str(t)) for t in t_grid)
def main_1():
    for line in stdin:
        args = {}
        # error message for missing files
        if 'Errno 2' in line:
            fname = line.split()[-1]
        # # reads output where I've scanned for empty files in command line
        # elif 'dam2214' in line:
        #     fname = line.split()[-1]
        # reads error messages as I've set them to print out in bkgdmap.py
        else:
            fname = line.split()[1].split('=')[1].split('/')[-1]
        # print fname
        l = fname.split('.')
        ch, an, t, pi, plen = l[1], l[2], l[4], l[5], l[6]
        args['cx'] = ch[3:]
        args['an'] = an
        args['tx'] = rev_t[t]
        args['px'] = pi[1:]
        args['pl'] = plen
        args['fi'] = ffmt.format(**args)
        print cfmt.format(**args)


def main_2():
    if len(argv) < 2:
        print 'usage: resend_jobs <anno>'
        exit()
    anno = argv[1]
    fdir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/precalc/{}'.format(anno)
    ffmt1 = fdir + '/AA_Map_10kb.{}.{}.t{}.merge1e07.bkgd'
    ffmt2 = fdir + '/AA_Map_10kb.{}.{}.t{}.merge1e07.bkgd'
    args = {}
    for tx in [-4.5, -4, -3.5, -3, -2.5, -2]:
        tstr = '{:.8f}'.format(10**tx)
        for c in range(1,23):
            ch = 'chr{}'.format(c)
            f1 = ffmt1.format(ch, anno, tstr)
            f2 = ffmt2.format(ch, anno, tstr)
            if not (os.path.isfile(f1) and os.path.isfile(f2)):
                args['fi'] = ffmt.format(an=anno)
                args['cx'] = c
                args['tx'] = -1*tx
                args['an'] = anno
                args['pl'] = '1e07'
                print cfmt2.format(**args)
                # print f1

if __name__ == '__main__':
    # main_1()
    main_2()
