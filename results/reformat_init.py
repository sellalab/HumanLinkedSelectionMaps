import re
from classes.runstruct import root_dir


__author__ = 'davidmurphy'


def reformat_init(fname):
    """
    A function to reformat older files by removing obsolete attributes,
    including new attributes and updating attributes whose names have changed
    :param fname: init_file/final_file for RunStruct based classes
    """
    # read the complete file as one string
    with open(fname, 'r') as f:
        text = f.read()

    # remove this, it is no longer being used
    text = text.replace("--naid='song'\n", "")
    # insert bkgd scale param
    text = text.replace("--bs_annos", "--bkgd_scale=100\n--bs_annos")
    # insert bmap/cmap dirs
    text = text.replace("--cs_annos=()",
                        "--bmap_dir='bmaps'\n--cs_annos=()\n--cmap_dir='cmaps'")
    # replace tvals/svals with bs_dfe, cs_dfe
    text = text.replace('tvals', 'bs_dfe')
    text = text.replace('svals', 'cs_dfe')
    text = text.replace("--vars.percentile='0.95'", "--vars.percentile=95")
    # remove obsolete vars attributes
    text = re.sub('--vars\.cutoff=.+\n', '', text)
    # remove obsolete keys from savekeys list
    text = text.replace("'naid', ", "")
    text = text.replace("'bs_annos', ",
                        "'bkgd_scale', 'bs_annos', 'bmap_dir', ")
    text = text.replace("'cs_annos', ", "'cs_annos', 'cmap_dir', ")

    with open(fname + '.tmp', 'w') as f:
        f.write(text)

fdir = root_dir + '/result/final_files/pr90to99clean'
fname = fdir + '/YRI.pr95.cleanrun.BS1.6.CS0.0.170714065758.final.txt'
reformat_init(fname)
