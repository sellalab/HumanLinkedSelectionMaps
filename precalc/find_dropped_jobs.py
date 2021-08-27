__author__ = 'davidmurphy'


import os
import subprocess
from sys import argv, stdin
import classes.runstruct as rst


run_path = '/ifs/data/c2b2/gs_lab/dam2214/run'
bsprog = run_path + '/sh/calc_bkgd.sh'
fexe = rst.root_dir + '/debug_bkgd/src/calc_bkgd'
bpart_prog = run_path + '/sh/prcs_bparts.sh'


# complete command string for bmap command
cfmt1 = 'qsub -l mem=2G,time={h}:: -cwd -j y -o bkgd-batch-call.log {prog} \
{finit} {ch} {t} {an} {fx} {p} {plen} 0 >> {an}.bjobs.log'

# complete command string for processing b partitions
cfmt2 = 'qsub -l mem=3G,time=:20: -cwd -j y -o {log} {prog} {ch} {gm} {bdir} \
 {an} {t} {plen}'

lfmt = '{ch}.{gm}.{bdir}.{an}.{t}.{plen}.processbmap.log'


def recover_error_jobs(err_log):
    # initialize dict for command string
    cdict = {}
    # set the path to sh scripts to be called by subprocess
    cdict['prg_1'] = '/ifs/data/c2b2/gs_lab/dam2214/run/sh/bkgd_rerun.sh'
    cdict['prg_2'] = '/ifs/data/c2b2/gs_lab/dam2214/run/sh/cprcs_bparts.sh'
    # parse error log to create rerun command string for each error
    with open(err_log, 'r') as f:
        for line in f:
            if 'IOError' in line:
                l = line.split('.')
                l1 = l[1].split('_')
                cdict['ch'] = l[0]
                cdict['sp'] = l1[0]
                cdict['pc'] = l1[1][4:]
                cdict['pn'] = l1[3][4:]
                if l[4] == 'processbmap':
                    cdict['ti'] = '{}.{}'.format(*l[2:4])
                    cdict['pi'] = l[10][1:]
                else:
                    cdict['ti'] = l[2]
                    cdict['pi'] = l[9][1:]
                # create calc_bkgd command string using parsed data
                cmd_1 = '{prg_1} {ch} {sp} {pc} {pn} {ti} {pi}'.format(**cdict)
                # create process bmap command with the parsed data
                cmd_2 = '{prg_2} {ch} {sp} {pc} {pn} {ti}'.format(**cdict)

                # print cmd_1

                # execute the command in the shell
                subprocess.call(cmd_2, shell=True)


def get_imax(chrom, block):
    """get indices for given chrom cut into x-length blocks"""
    return int(rst.chromosome_length(chrom) / block)


def main_1():
    if len(argv) != 6:
        print 'usage: find_dropped_jobs <init_file> <job_limit> <ti> <tj> <td>'
        exit(1)

    # path to rerun script
    f_sh = '/ifs/data/c2b2/gs_lab/dam2214/run/sh/bkgd_rerun.sh'

    # path to init file
    f_init = argv[1]
    assert os.path.isfile(f_init)
    cst = rst.ChromStruct(chrom='chr1', init=f_init)

    # set annotation used in bmaps and bdir
    an = cst.bdir

    # get list of files that DO exist in the given dir
    bdir = '{}/precalc/{}/'.format(rst.root_dir, an)
    b_files = [f for f in os.listdir(bdir) if f.endswith('.bkgd')]

    # template for filename: ch=chrom; t=coeff; p=block_idx; l=block_length
    f_fmt = 'AA_Map_10kb.{c}.{a}.t{t}.p{i}.{l}.bkgd'
    # template for shell command to rerun bmap
    c_fmt = '{sh} {f} {c} {a} {t} {i}'

    # negative exponents for t values
    # tex = [4.5, 4, 3.5, 3, 2.5, 2]
    ti, tj, td = map(float, argv[3:6])
    tex = rst.np.arange(ti, tj, td)

    # create list of expected filenames
    bl = 1e7
    l_str = '{:.0e}'.format(bl).replace('+', '')
    job_count = 0
    job_limit = int(argv[2])
    for ch in rst.human_autosomes:
        for tx in tex:
            t = 10 ** -tx
            t_str = '{:.8f}'.format(t)
            for ix in xrange(get_imax(ch, bl)+1):
                f = f_fmt.format(c=ch, a=an, t=t_str, i=ix, l=l_str)
                cmd = c_fmt.format(sh=f_sh, f=f_init, c=ch, a=an, t=tx, i=ix)
                under_limit = (job_count <= job_limit)
                if f not in b_files and under_limit:
                    subprocess.call(cmd, shell=True)
                    job_count += 1
                if not under_limit:
                    break

def main_2():
    with open(stdin, 'r') as pipe:
        for line in pipe:
            print 'piped: {}'.format(line)


def main_3():
    # AA_Map_10kb.chr11.cadd94_gmask_v1.6_without_bstat_nonexonic.t0.00010000.p9.1e07.bkgd
    bdir = argv[2]
    args = {}
    args['finit'] = argv[1]
    # args['an'] = argv[3]
    args['prog'] = bsprog
    args['fx'] = fexe
    for f in os.listdir(bdir):
        if f.startswith('AA_Map_10kb') and ('merge' not in f):
            fl = f.split('.')
            args['an'] = '.'.join(fl[2:-5])
            ch = fl[1]  # get chrom
            splen = fl[-2]  # get plen string
            plen = float(splen)
            spidx = fl[-3][1:]  # get partition index string
            pidx = int(spidx)
            maxp = get_imax(ch, plen)  # get max partition
            stv = '0.'+fl[-4]
            # if the partition is not the last one, expected size is plen
            if pidx < maxp:
                expected = plen
            # otherwise calculate the length of the final partition
            else:
                expected = rst.chromosome_length(ch) - (plen * maxp)

            # check the length of the file
            file_len = 0
            with open(bdir + '/' + f, 'r') as fread:
                for line in fread:
                    if not line.startswith('#'):
                        line = line.split()
                        if len(line) ==2:
                            file_len += int(line[1])

            if file_len != expected:
                msg = '{} only {} bp, expected {}. resending partition.'
                # print msg.format(f, file_len, expected)
                args['t'] = stv
                args['plen'] = int(plen)
                args['h'] = 72
                args['ch'] = ch
                args['p'] = pidx
                cmd = cfmt1.format(**args)
                subprocess.call(cmd, shell=True)
                # print cmd


if __name__ == '__main__':
    # if len(argv) != 2:
    #     print 'usage: recover_error_jobs <err_log>'
    #     exit(1)
    # err_log = argv[1]
    # recover_error_jobs(err_log)
    main_3()