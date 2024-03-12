__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from likelihood.cllh_functions import predicted_pi
from classes.runstruct import ChromStruct, root_dir, cst_from_fldr
from precalc.lh_inputs import load_saved, adjust_arrays, get_jackknife_mask


init_dir = root_dir + '/result/init_files'
final_dir = root_dir + '/result/final_files'
ffmt = init_dir + '/YRI.{an}.BS1.6.CS0.0.NOT_STARTED.initial.txt'


def jackknife_prediction_vector(fl_name):
    """
    assemble prediction vector from missing jackknife regions using params
    from jackknife sample where region was removed
    """
    # set path to main folder of each jackknife result
    fdir = '{}/result/final_files/{}/'.format(root_dir, fl_name)

    # use the general purpose init file for the jackknife folder to load arrays
    f_init = fdir + 'jk_init.txt'.format(fl_name)
    cst = ChromStruct(chrom='chr1', init=f_init)

    # load complete array data ONCE
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # save predictions for each jackknife sample in dict with index as keys
    pred_dict = {}
    for fldr in os.listdir(fdir):
        f_path = fdir + fldr
        if os.path.isdir(f_path):
            f_list = os.listdir(f_path)
            # check for contents in folder
            if len(f_list):
                c = [f for f in f_list if 'composite' in f]
                # should only have one composite file
                if len(c) > 1:
                    print("MULTIPLE COMPOSITES! {}".format(f_path))

                # initialize RunStruct with composite file
                f_init = '{}/{}'.format(f_path, c[0])
                cst = ChromStruct('chr1', init=f_init)

                # # switch off jackknife flag to load data normally
                # cst.vars.use_jackknife = False
                # # load complete array data
                # sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

                # get jackknife mask COMPLEMENT for sample
                jidx = cst.vars.jackknife_index
                jwin = cst.vars.jackknife_window
                m = ~get_jackknife_mask(sg, jidx, jwin)

                # print('upload?')
                # get JUST jackknife region using mask
                jsg, jnu, jnt, jdv, jpl = sg[m], nu[m], nt[m], dv[m], pl[m]
                if bs is not None:
                    jbs = bs[m]
                else:
                    jbs = None
                if cs is not None:
                    jcs = cs[m]
                else:
                    jcs = None
                args = (cst, jbs, jcs, jnu, jnt, jdv, jpl)

                # mask and rescale maps from jackknife region
                jbs, jcs, jnu = adjust_arrays(*args)[:3]

                # calculate predicted pi in jackknife region
                pred_dict[jidx] = predicted_pi(cst.params, cst, jnu, jbs, jcs)

    # create list of arrays in index order
    pred_list = []
    for k in sorted(pred_dict.keys()):
        pred_list.append(pred_dict[k])

    # concatenate predictions and save as single array
    fpred = fdir + 'jkpred.npy'
    np.save(fpred, np.concatenate(pred_list))


def jackknife_prediction_vector_2(anno):
    """
    assemble prediction vector from missing jackknife regions using params
    from jackknife sample where region was removed
    """
    # use init file for the jackknife folder to load arrays ONCE
    # f_init = ffmt.format(an=anno)
    # cst = ChromStruct(chrom='chr1', init=f_init)
    cst = cst_from_fldr(anno)
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # save predictions for each jackknife sample in dict with index as keys
    pred_list = []
    # make a separate list for chunks of uniform prediction for sorted pred analysis
    uniform_pred_list = []
    for jkidx in range(1441):
        ji = '{:04}'.format(jkidx)
        # set foldrer path to current jackknife index
        # fldr = '{an}_jkidx_{ji}'.format(an=anno, ji=ji)
        fldr = '{an}_jackknife_results/{an}_jkidx_{ji}'.format(an=anno, ji=ji)
        fpath = final_dir + '/' + fldr
        # skip jkidx folder paths that dont exist
        if (not os.path.isdir(fpath)) or (len(os.listdir(fpath)) == 0):
            continue
        # get list of files in the folder and find the "composite" file
        f_list = os.listdir(fpath)
        c = [f for f in f_list if 'composite' in f]
        # should only have one composite file
        if len(c) > 1:
            print("MULTIPLE COMPOSITES! {}".format(fpath))

        # initialize RunStruct with composite file
        f_jk = '{}/{}'.format(fpath, c[0])
        jk_cst = ChromStruct('chr1', init=f_jk)

        # get jackknife mask COMPLEMENT for sample
        jidx = jk_cst.vars.jackknife_index
        jwin = jk_cst.vars.jackknife_window
        m = ~get_jackknife_mask(sg, jidx, jwin)

        # get JUST jackknife region using mask
        jsg, jnu, jnt, jdv, jpl = sg[m], nu[m], nt[m], dv[m], pl[m]

        if bs is not None:
            jbs = bs[m]
        else:
            jbs = None
        if cs is not None:
            jcs = cs[m]
        else:
            jcs = None
        args = (jk_cst, jbs, jcs, jnu, jnt, jdv, jpl)

        # mask and rescale maps from jackknife region
        jbs, jcs, jnu = adjust_arrays(*args)[:3]
        jnu_const = np.full(len(jnu), cst.stat.meandiv)

        # calculate predicted pi in jackknife region and add to list of blocks
        jk_pred = predicted_pi(jk_cst.params, jk_cst, jnu, jbs, jcs)
        pred_list.append(jk_pred)
        # get predicted pi without in jackknife region without scaling by local substitution rate
        jk_pred_const = predicted_pi(jk_cst.params, jk_cst, jnu_const, jbs, jcs)
        uniform_pred_list.append(jk_pred_const)

    # concatenate predictions and save as single array
    save_dir = final_dir + '/{an}_jackknife_results/'.format(an=anno)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    fpred = save_dir + '{an}.jkpred.npy'.format(an=anno)
    np.save(fpred, np.concatenate(pred_list))
    fpred_const = save_dir + '{an}.jkpred_const.npy'.format(an=anno)
    np.save(fpred_const, np.concatenate(uniform_pred_list))


def main():
    # if len(argv) != 2:
    #     print 'usage: jackknife_prediction <folder_name>'
    #     exit(1)
    # jackknife_prediction_vector(argv[1])
    if len(argv) != 2:
        print('usage: jackknife_prediction <anno>')
        exit(1)
    jackknife_prediction_vector_2(argv[1])


if __name__ == '__main__':
    main()
