from sys import argv
from classes.datasampler import RunStruct, DataSampler, os

__author__ = 'davidmurphy'


def init_files_batch(sampler, new_sample=False):
    """
    Generate and save a batch of bootstrap samples using data and bootstrap params encoded by the RunStruct.
    :param sampler: a DataSampler class for making bootstrap/jackknife samples
    :param new_sample: flag indicates that a new directory should be created and a new set of samples drawn
    :type sampler: DataSampler
    """
    assert sampler.sample_type in 'boot jack'.split()
    # add bt/jk numbers to this template for each file produced
    file_template = sampler.rst.txt_file
    # get the index once in order to update the jackknife
    sampler.get_idx()
    if sampler.sample_type == 'boot':
        for bi in range(sampler.rst.vars.bootstrap_samples):
            new_file = file_template.replace('.txt', '-bt{}.txt'.format(bi + 1))
            sampler.next_boot()
            sampler.get_boot()
            sampler.rst.save(txt_file=new_file)
    # jackknife files
    else:
        for ji in range(sampler.rst.vars.jackknife_samples):
            new_file = file_template.replace('.txt', '-jk{}.txt'.format(ji + 1))
            sampler.next_jack()
            sampler.rst.save(txt_file=new_file)

    return None


def main_local():
    r = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/init_files/'
    f = r + 'lamprey-cons-95-segs_BS1_CS0_initial.txt'
    t = 'jack'
    rst = RunStruct(init=f)
    rst.vars.sample_type = t
    init_files_batch(sampler=DataSampler(rst=rst), new_sample=False)


def main_remote():
    if len(argv) == 3 and argv[2] in 'boot jack'.split():
        rst = RunStruct(init=argv[1])
        rst.vars.sample_type = argv[2]
        init_files_batch(sampler=DataSampler(rst=rst), new_sample=False)
    else:
        print 'usage: bootstrap <init_file> <sample_type> [boot, jack]'


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
