__author__ = 'davidmurphy'


import os
import re
import sys
import shutil
import subprocess
from time import sleep
from classes.runstruct import ChromStruct, root_dir, np


init_dir = root_dir + '/result/init_files'
final_dir = root_dir + '/result/final_files'
run_dir  = '/ifs/data/c2b2/gs_lab/dam2214/run'
ffmt = init_dir + '/YRI.{an}.BS1.6.CS0.0.NOT_STARTED.initial.txt'


def find_missing_indices(directory):
    """find missing indices in results folder"""
    # set ape/CADD init file
    if 'ape' in directory:
        tkn = 'ape'
        f_init = init_dir + '/YRI.ape_cons94_clean_extel.filter.gmap.edge.0.1.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    else:
        tkn = 'cadd'
        f_init = init_dir + '/YRI.cadd93_extel.filter.gmap.edge.0.1.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    # regex to get the iprm string from each file
    re_idx = re.compile('iprm_\d\d')
    # list of existing files in the directory
    f_list = [f for f in os.listdir(directory) if f.endswith('.txt')]
    # record all of the iprm strings from existing files
    existing_iprm = []
    for f_name in f_list:
        iprm_match = re_idx.search(f_name)
        if iprm_match:
            existing_iprm.append(iprm_match.group())

    # get job params from one of the existing files
    cst = ChromStruct('chr1', init=directory + '/' + f_list[0])
    if cst.fixed.min_bs is not None:
        min_b = np.exp(cst.fixed.min_bs)
    else:
        min_b = None
    min_bsx = cst.fixed.min_bsx
    # command to send missing job:
    prg = '/ifs/data/c2b2/gs_lab/dam2214/run/sh/inf1.sh'
    log = '{}.minbsx_{}.minb_{}.idx{{i}}'.format(tkn, min_bsx, min_b)
    # usage: cluster_runinf <init> <idx> <min_bsx> <min_b
    job = '{} {} {{i}} {} {}'.format(prg, f_init, min_bsx, min_b)
    cmd = 'qsub -l mem=24G,time=16:: -cwd -j y -o {} {}'.format(log, job)

    # check against a list of all 15 iprm strings that SHOULD exist
    all_iprm = ['iprm_{:02}'.format(i) for i in xrange(15)]
    for idx in all_iprm:
        if idx not in existing_iprm:
            command = cmd.format(i=int(idx.split("_")[1]))
            subprocess.call(command, shell=True)


def move_to_folder(filepath):
    f_name = filepath.split('/')[-1]
    if f_name.startswith('YRI.') and f_name.endswith('.final.txt'):
        # if '_nval' in f_name.split('.')[1]:
        #     fldr = f_name.split('.')[1]
        # else:
        #     # fldr = '_'.join(f_name.split('.')[1:3])
        #     fldr = f_name.split('.')[1] + '_fixed'
        fldr = f_name.split('.')[1]
        fldrpath = filepath.replace(f_name, fldr)
        if not os.path.isdir(fldrpath):
            os.mkdir(fldrpath)
        shutil.move(filepath, fldrpath)
        # print fldrpath


def move_to_folder_2():
    for f in os.listdir(final_dir):
        if f.endswith('final.txt'):
            f_list = f.split('.')
            an = f_list[1]
            fldr = an
            filepath = '{}/{}'.format(final_dir, f)
            fldrpath = '{}/{}'.format(final_dir, fldr)
            if not os.path.isdir(fldrpath):
                os.mkdir(fldrpath)
            shutil.move(filepath, fldrpath)


def move_to_folder_3():
    # YRI.fish_cons94_gmask.BS1.6.CS0.0.iprm_11.jkidx_0087.200908214213.final
    for f in os.listdir(final_dir):
        if f.endswith('final.txt'):
            if '.jkidx_' in f:
                f_list = f.split('.')
                an = f_list[1]
                jk = f_list[7]
                fldr = '{}_{}'.format(an, jk)
                fldrpath = '{}/{}'.format(final_dir, fldr)
                filepath = '{}/{}'.format(final_dir, f)
                if not os.path.isdir(fldrpath):
                    os.mkdir(fldrpath)
                shutil.move(filepath, fldrpath)
            else:
                f_list = f.split('.')
                pop, an = f_list[0], f_list[1]
                fldr = '{}_{}'.format(pop, an)
                fldrpath = '{}/{}'.format(final_dir, fldr)
                filepath = '{}/{}'.format(final_dir, f)
                # print fldrpath + '\n' + filepath
                if not os.path.isdir(fldrpath):
                    os.mkdir(fldrpath)
                shutil.move(filepath, fldrpath)


def get_missing_iprms(flist):
    """get list of missing indices from a folder of final inf run files"""
    ilist = [int(f.split('.')[6].split('_')[1]) for f in flist]
    imiss = []
    for i in xrange(15):
        if i not in ilist:
            imiss.append(i)

    return imiss


def find_missing_jackknife_jobs(anno):
    """scan the jackknife jobs, find missing indices and rerun them"""
    # args dict for string formats
    args = {}
    # set annotation from function args
    args['an'] = anno
    # prefill the program path and initialization file path
    args['fi'] = ffmt.format(**args)

    # parameters for qsub command
    args['prm'] = 'mem=24G,time=10::'

    # string template for qsub commands
    scmd1 = 'qsub -l {prm} -cwd -j y -o {log} {prg} {fi} {i} {ji} > /dev/null'
    scmd2 = 'qsub -l {prm} -cwd -j y -o {log} {prg} {fldr} > /dev/null'

    # string template for logs
    slog1 = run_dir + '/{an}.idx{i}.jkidx{ji}.inf1.log'
    slog2 = run_dir + '/{fldr}.inf2.log'

    # scan each jackknife index
    for jkidx in xrange(1441):
        # set the current arg value for jkidx
        args['ji'] = '{:04}'.format(jkidx)
        # set foldrer path to current jackknife index
        fldr = '{an}_jkidx_{ji}'.format(**args)
        args['fldr'] = fldr
        fpath = final_dir + '/' + fldr
        # skip jkidx folder paths that dont exist
        if not os.path.isdir(fpath):
            continue
        # get the list of files in paths that do exist and count final files
        flist = [f for f in os.listdir(fpath) if f.endswith('.txt')]
        num = len(flist)

        # if final files total < 15, find the missing indices
        if num < 15:
            # continue
            # set program to inf step 1 if rerunning step 1 jobs
            args['prg'] = run_dir + '/sh/jkinf1.sh'
            imiss = get_missing_iprms(flist)
            for i in imiss:
                # set arg value for current missing index
                args['i'] = i
                # set current arg value for log based on idx, jkidx
                args['log'] = slog1.format(**args)
                # set the command string using args dict
                cmd = scmd1.format(**args)
                print cmd
                # sys.stderr.write('{} {}\n'.format(fldr, i))
                # sys.stdout.flush()
                # subprocess.call(cmd, shell=True)
                # sleep(240)

        # if the folder contains 15 files, run inf step 2
        elif num == 15:
            # set program to inf step 2 if running step 2
            args['prg'] = run_dir + '/sh/jkinf2.sh'
            args['log'] = slog2.format(**args)
            cmd = scmd2.format(**args)
            print cmd
            # sys.stderr.write('{}\n'.format(fldr))
            # sys.stdout.flush()
            # subprocess.call(cmd, shell=True)
        else:
            n_composite = len([f for f in flist if 'composite' in f])
            if n_composite == 0:
                print fldr
            # continue


def compare_fish_cadd():
    an1 = 'fish_cons94_gmask'
    an2 = 'cadd94_gmask'

    # scan each jackknife index
    for jkidx in xrange(1441):
        # set the current arg value for jkidx
        ji = '{:04}'.format(jkidx)
        # set fish/cadd folder paths to current jackknife index
        fldr1 = '{}_jkidx_{}'.format(an1, ji)
        fldr2 = '{}_jkidx_{}'.format(an2, ji)
        fpath1 = final_dir + '/' + fldr1
        fpath2 = final_dir + '/' + fldr2

        # skip jkidx fish folder paths that dont exist
        if not os.path.isdir(fpath1):
            continue
        # print names of jkidx where fish folder exists but cadd doesn't
        else:
            if not os.path.isdir(fpath2):
                print fldr2


def main():
    move_to_folder_2()
    # if len(sys.argv) == 2:
    #     find_missing_jackknife_jobs(sys.argv[1])
    # else:
        # move_to_folder_3()
        # compare_fish_cadd()
    # if len(sys.argv) != 2:
    #     print 'usage: find_missing_inf1_runs <filepath>'
    #     exit(1)
    # filepath = sys.argv[1]
    # move_to_folder_2(filepath)


if __name__ == '__main__':
    main()

