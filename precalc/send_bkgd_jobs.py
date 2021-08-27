__author__ = 'davidmurphy'


import re
import os
import sys
import glob
import time
import shutil
import subprocess
from likelihood.send_jackknife_jobs import jobcount, errmsg
from classes.runstruct import ChromStruct, root_dir, chromosome_length


run_path = '/ifs/data/c2b2/gs_lab/dam2214/run'
bsprog = run_path + '/sh/calc_bkgd.sh'
fexe = root_dir + '/debug_bkgd/src/calc_bkgd'
bpart_prog = run_path + '/sh/prcs_bparts.sh'


# complete command string for bmap command
cfmt1 = 'qsub -l mem=2G,time={h}:: -cwd -j y -o bkgd-batch-call.log {prog} \
{finit} {ch} 10**-{t} {an} {fx} {p} {plen} 0 >> {an}.bjobs.log'

# complete command string for processing b partitions
cfmt2 = 'qsub -l mem=3G,time=:20: -cwd -j y -o {log} {prog} {ch} {gm} {bdir} \
 {an} {t} {plen}'

lfmt = '{ch}.{gm}.{bdir}.{an}.{t}.{plen}.processbmap.log'


def send_command(command):
    """send command unless there are too many jobs - then wait 5min"""
    njb = jobcount()
    while njb > 900:
        msg = 'There are {} jobs. Pausing for 5min.'.format(njb)
        errmsg(msg)
        time.sleep(300)
        njb = jobcount()
    subprocess.call(command, shell=True)


def main():
    if len(sys.argv) == 3:
        args = {}
        args['prog'] = bsprog
        args['fx'] = fexe
        args['finit'] = sys.argv[1]
        args['an'] = sys.argv[2]
        print 'making maps for {}'.format(args['an'])
        t_vals = [4.5, 4, 3.5, 3, 2.5, 2]
        p_vals = map(int, [1e07, 2e07, 3e07, 3e07, 3e07, 3e07])
        h_vals = [24, 24, 24, 24, 24, 24]
        for i in range(6):
            args['t'] = t_vals[i]
            args['plen'] = p_vals[i]
            args['h'] = h_vals[i]
            for c in range(1, 23):
                chrom = 'chr{}'.format(c)
                args['ch'] = chrom
                max_p = 1 + int(chromosome_length(chrom) / p_vals[i])
                for p in range(0, max_p):
                    args['p'] = p
                    cmd = cfmt1.format(**args)
                    send_command(cmd)
                    # print cmd

    elif len(sys.argv) == 2:
        args = {}
        args['prog'] = bpart_prog
        args['gm'] = 'AA_Map_10kb'
        args['bdir'] = sys.argv[1]
        args['an'] = sys.argv[1]
        print 'consolidating maps for {}'.format(args['an'])
        t_vals = [4.5, 4, 3.5, 3, 2.5, 2]
        long_segs = ['ape_cons94_new', 'primate_cons94_new',
                     'prosimian_cons94_new']
        if args['an'] in long_segs:
            p_vals = '2e07 2e07 3e07 3e07 3e07 3e07'.split()
        else:
            p_vals = '1e07 1e07 3e07 3e07 3e07 3e07'.split()

        for i in range(6):
            args['t'] = '10**-{}'.format(t_vals[i])
            args['plen'] = p_vals[i]
            for c in range(1, 23):
                chrom = 'chr{}'.format(c)
                args['ch'] = chrom
                args['log'] = lfmt.format(**args)
                cmd = cfmt2.format(**args)
                # print cmd
                subprocess.call(cmd, shell=True)
                # send_command(cmd)
    else:
        print 'usage_1: send_bkgd_jobs <f_init> <anno>'
        print 'usage_2: process_partitions <anno>'
        exit(1)


if __name__ == '__main__':
    main()
