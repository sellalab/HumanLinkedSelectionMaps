from classes.runstruct import ChromStruct
from datetime import datetime as dt
from sys import argv
import subprocess
import os

__author__ = 'davidmurphy'


def calc_bkgd(cst, coef, anno, udel, fexe, pidx=-1, plen=0, pts=0):
    """
    Write configuration file to McVicker's C code and then send shell command
    to run it
    :param cst: ChromStruct object
    :param coef: selection coefficient for mutations
    :param anno: annotation for file lookup
    :param udel: u deleterious initial for bmap
    :param fexe: path to calc_bkgd executable file
    :param csave: flag to save the cfg file, deleted by default
    """
    # hard-code sum approx flag to 0
    aprx = 0
    # hard-code config file delete at end of run
#    cfg_save = False
    
    # file names and string templates
    fcnfg = '{root}/result/configs/{chrom}.{a}.{gmap}.t{t:.8f}.conf'
    fpath = '{root}/precalc/{bdir}/'.format(**cst.dict)
    # fgmap = '{root}/data/maps/{gmap}/genmap_{gmap}_{chrom}.txt'
    fgmap = '{root}/data/maps/{gmap}/{chrom}_{gmap}.txt'
    # fanno = '{root}/data/bsanno/{a}/longest_transcript_533_{a}_{chrom}.coords'
    fanno = '{root}/data/bsanno/{a}/{chrom}.{a}.bed'
    fsave = '{gmap}.{chrom}.{a}.t{t:.8f}'.format(a=anno, t=coef, **cst.dict)
    # fsave = '{gmap}_{chrom}_{a}_t{t:.8f}'.format(a=anno, t=coef, **cst.dict)
    
    # create config file header
    h = ['CHROMOSOME_NAME={}'.format(cst.chrom),
         'CHROMOSOME_LENGTH={}'.format(cst.chlen),
         'CHR_PARTITION_INDEX={}'.format(pidx),
         'CHR_PARTITION_LENGTH={:.0e}'.format(plen),
         'OUTPUT_DIR={}'.format(fpath),
         'OUTPUT_TOKEN={}'.format(fsave),
         'RECOMB_RATE_TABLE=' + fgmap.format(**cst.dict),
         'RECOMB_RATE_SCALE={:g}'.format(1e-8),
         'CONS_TABLE=' + fanno.format(a=anno, **cst.dict),
         'PARAM_T={:g}'.format(coef),
         'PARAM_U={:g}'.format(udel),
         'PARAM_T_DIST_TYPE=POINT',
         'BKGD_SCALE={:.6f}'.format(cst.bscl),
         'BKGD_CHANGE_THRESH={}'.format(2.0 / cst.bscl),
         'BKGD_OFFSET={}'.format(1.0 / cst.bscl),
         'BKGD_INTERP_MAX_DIST=0.0001',
         'ENFORCE_BGI_PTS={}'.format(pts),
         'USE_SUM_APPROXIMATION={}'.format(aprx),
         'BKGD_PARAM_MAX_SUM_THRESH=0.001',
         'CALC_BKGD_EXECUTABLE={}'.format(fexe)]
    header = '\n'.join(h)

    # add index suffices to partitioned map files 
    cfg = fcnfg.format(a=anno, t=coef, **cst.dict)
    if (pidx >= 0) and (plen > 0):
        # cfg file
        suf = '.p{}.conf'.format(pidx)
        cfg = cfg.replace('.conf', suf)
        # map file
        suf = '.p{}.{:.0e}'.format(pidx, plen).replace('+', '')
        fsave += suf

    # write header to config file
    with open(cfg, 'w') as f:
        f.write(header + '\n')

    # create a path for new map folders if needed
    if not os.path.isdir(fpath):
        os.mkdir(fpath)

    # path for log files (create if it does not exist already)
    log_path = '{}/bkgd_run_logs/'.format(os.getcwd())
    if not os.path.isdir(log_path):
        os.mkdir(log_path)
    # log file name for current run
    log = log_path + '{}.{}.log'.format(cst.tkn, fsave)
    
#    print 'fsave={}'.format(fsave)
#    print 'log={}'.format(log)
#    print 'cfg={}'.format(cfg)

    # record the time just before calling calc_bkgd
    t_0 = dt.now()

    # format shell command and call calc_bkgd as a subprocess in the shell
    if cst.root.startswith('/Users/davidmurphy'):
        cmd_fmt = '{x} {c} {f}'
        command = cmd_fmt.format(x=fexe, c=cfg, f=fsave)    
        subprocess.call(command, shell=True)
        
    else:
        cmd_fmt = '{x} {c} {f} &> {g}'
        command = cmd_fmt.format(x=fexe, c=cfg, f=fsave, g=log)
        subprocess.call(command, shell=True)

        # append the final runtime to the log
        t_tot = dt.now() - t_0
        with open(log, 'a') as f:
            f.write('runtime={}\n'.format(t_tot))

    # # clean-up the configuration file, which is no longer needed
#    if not cfg_save:
#        os.remove(cfg)


def main_local():
#    ch = 'chr2'
    bd = 'test_mod'
    sc = 100.0
    tk = 'called.in.py'
    an = 'primate_cons95_Segments'
    udel = 7.4e-08   

    args = ['', 'chr10', '10**-4.5', '0', '1e06', '1']
# =============================================================================
#     if len(argv) != 6:
#         print 'usage: calc_bkgd <chrom> <coef> <pidx> <plen> <pts>'
#         exit(1)
# =============================================================================
        
    cst = ChromStruct(chrom=args[1], bdir=bd, bscl=sc, tkn=tk)
    fexe = '{}/debug_bkgd/src/calc_bkgd'.format(cst.root)
    coef, pidx, plen, pts = map(eval, args[2:])
    
    calc_bkgd(cst, coef, an, udel, fexe, pidx, plen, pts)
    
# =============================================================================
#     chrom = 'chr14'    
#     cst = ChromStruct(chrom=chrom, bmap_dir=bd, bscl=sc, tkn=tk)
#     coef = 0.01
#     pidx = 0
#     plen = 1e7
#     pts = 0
#  
#     if pts:
#         cst.tkn += '_pts'
#         cst.bmap_dir += '_pts'
#         cst.reset()
#  
#     for pidx in xrange(int(cst.chlen / plen)):
#         calc_bkgd(cst, coef, an, udel, fexe, pidx, plen, pts)
# =============================================================================
 

def main_remote():
    if len(argv) != 9:
        print 'usage: calc_bkgd <init> <chrom> <coef> <anno> <fexe>' \
              ' <pidx> <plen> <pts>'
        exit(1)
        
    cst = ChromStruct(init=argv[1], chrom=argv[2])
    coef = eval(argv[3])
    anno = argv[4]
    udel = cst.fixed.u_fix
    # udel = 1.0  # TODO: don't leave this alone...
    fexe = argv[5]
    pidx = eval(argv[6])
    plen = eval(argv[7])
    pts = eval(argv[8])
    if pts:
        cst.tkn += '_pts'
        cst.bdir += '_pts'
        cst.reset()
    calc_bkgd(cst, coef, anno, udel, fexe, pidx, plen, pts)


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy'):
        main_local()
    else:
        main_remote()
