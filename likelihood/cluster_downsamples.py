import numpy as np
from sys import argv
from shutil import move
from runinf import run_inference, evaluate
from classes.runstruct import ChromStruct, root_dir
from simulations import cluster_sim_config as simcfg

__author__ = 'davidmurphy'


"""
Usage: set initial parameters manually and run parallelized optimization on the
cluster.
"""


def main():
    if root_dir.startswith('/Users/davidmurphy'):
        # msk_id = '0000'
        iprm_idx = 0
        meth = simcfg.optim_meth

    else:
        if len(argv) != 3:
            print 'usage: cluster_downsamples <iprm_idx> <method>'
            exit(1)

        iprm_idx = int(argv[1])
        meth = (argv[2],)

    # create new chrom struct from config params
    cst = ChromStruct(chrom='chr1', bdir=simcfg.bdr, tkn=simcfg.tkn,
                      methods=meth)

    # set bounds and inits
    cst.fixed.min_bs = simcfg.min_bs
    cst.fixed.min_red = simcfg.min_red
    # cst.init_params = simcfg.initial_params
    cst.init_params = np.array(simcfg.iprm[iprm_idx])
    cst.params = cst.init_params
    # cst.msk_id = msk_id
    cst.vars.down_sample = simcfg.downsample

    # set nt files to sim names
    cst.files.fnt = cst.nt_simfiles(simcfg.newtkn)

    # run
    # TODO: find a more permanent solution for switching grad on/off
    grad = True if meth[0] in simcfg.grad_methods else False
    run_inference(cst, parallel=True, grad=grad)
    # evaluate(cst, True)

    # run ID for final file
    run_id = 'iprm_{:02}'.format(iprm_idx)
    fnstr = '{rslt}/final_files/{neut}.{lbl}.{rid}.{timestamp}.final.txt'
    new_fnm = fnstr.format(rid=run_id, **cst.dict)

    # replace default final file with final file including run ID
    old_fnm = cst.final_file
    move(old_fnm, new_fnm)


if __name__ == '__main__':
    main()
