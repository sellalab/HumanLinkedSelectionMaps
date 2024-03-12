__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from classes.runstruct import ChromStruct, root_dir


init_dir = root_dir + '/result/init_files'
final_dir = root_dir + '/result/final_files'
ffmt = init_dir + '/YRI.{an}.BS1.6.CS0.0.NOT_STARTED.initial.txt'


def jack_params(anno):
    # set path for saving jackknife results
    save_dir = final_dir + '/{an}_jackknife_results/'.format(an=anno)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # create list of filenames for saving or loading presaved data
    ftokens = ['pmf.npy', 'udl.npy', 'pi0.npy', 'clh.npy']

    # save data arrays generated from original files if flagged
    pmf = []
    udl = []
    pi0 = []
    clh = []
    for jkidx in range(1441):
        ji = '{:04}'.format(jkidx)
        # set foldrer path to current jackknife index
        fldr = '{an}_jkidx_{ji}'.format(an=anno, ji=ji)
        fpath = final_dir + '/' + fldr
        # skip jkidx folder paths that dont exist
        if (not os.path.isdir(fpath)) or (len(os.listdir(fpath)) == 0):
            continue
        # get list of files in the folder and find the "composite" file
        f_list = os.listdir(fpath)
        c = [f for f in f_list if 'composite' in f]
        # should only have one composite file
        if len(c) > 1:
            print "MULTIPLE COMPOSITES! {}".format(fpath)

        # initialize RunStruct with composite file
        f_jk = '{}/{}'.format(fpath, c[0])
        cst = ChromStruct('chr1', init=f_jk)

        # get dfe, udel, pi0 and CLLH from composite run
        pmf.append(cst.uvec[0] * 1e8)
        udl.append(sum(cst.uvec[0]) * 1e8)
        pi0.append(cst.params[-1] / cst.fixed.tau_init)
        clh.append(cst.stat.best_lh)

    # convert to arrays and save
    dtlists = [pmf, udl, pi0, clh]
    for ft, dt in zip(ftokens, dtlists):
        fsave = save_dir + ft
        np.save(fsave, np.array(dt))


def main():
    if len(argv) != 2:
        print 'usage: jackknife_params <anno>'
        exit(1)
    jack_params(argv[1])


if __name__ == '__main__':
    main()