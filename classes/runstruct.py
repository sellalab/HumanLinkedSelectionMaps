# todo: create much more extensive log writers for data processing steps
# todo: improve "filter" annotations processing - just dump all in one folder?

from data_processing.data_tools import chromosome_length, str2time
from optimizerparams import OptimizerParams
from fixedparams import FixedParams
from collections import namedtuple
from statistics import Statistics
from filestruct import FileStruct
from variables import Variables
from dictlike import DictLike
from itertools import izip
from time import strftime
import numpy as np
import re
import os


__author__ = 'davidmurphy'


"""
GLOBAL VARIABLES FOR DEFAULT SETTINGS
"""
# GOOGLE DRIVE DIRECTORY
goog_dir = '/Users/davidmurphy/GoogleDrive'

# ROOT DIRECTORIES FOR MY COMPUTER AND CLUSTER
if os.getcwd().startswith('/Users/davidmurphy'):
    root_dir = '{}/linked_selection'.format(goog_dir)
else:
    root_dir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection'

# HUMAN AUTOSOMES
human_autosomes = tuple('chr{}'.format(c) for c in xrange(1, 23))

# DROSOPHILA AUTOSOMES
dmel_autosomes = ('2L', '2R', '3L', '3R')

# DEFAULT FITNESS EFFECTS GRID
default_fitness_effects = np.power(10.0, np.arange(-4.5, -1.5, 0.5))

# DEFAULT FILTERED REGIONS
default_filter_file = root_dir + '/data/coords/neutral_filter_files.txt'

# DEFAULT PATH TO COMPILED CALC_BKGD EXECUTABLE FILE
calc_bkgd_executable = '{}/bkgd/src/calc_bkgd'.format(root_dir)

# FIGURES DIRECTORY
figs_dir = '{}/murphy_et_al/figures'.format(root_dir)


def cst_from_fldr(fldr, ch='chr1'):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    f_init = fdir + flist[0]

    return ChromStruct(ch, init=f_init)


def chromosome_length(chrom):
    """return the integer chromosome length in bp for chrom"""
    hg18 = dict(chr1=247249719, chr2=242951149, chr3=199501827,
                chr4=191273063, chr5=180857866, chr6=170899992,
                chr7=158821424, chr8=146274826, chr9=140273252,
                chr10=135374737, chr11=134452384, chr12=132349534,
                chr13=114142980, chr14=106368585, chr15=100338915,
                chr16=88827254, chr17=78774742, chr18=76117153,
                chr19=63811651, chr20=62435964, chr21=46944323,
                chr22=49691432, chrX=154913754, chrY=57772954)
    chr_length = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
                  'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
                  'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431,
                  'chr10': 135534747, 'chr11': 135006516,
                  'chr12': 133851895, 'chr13': 115169878,
                  'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
                  'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983,
                  'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566,
                  'chrY': 59373566, 'chrX': 155270560,
                  # dmel chroms
                  '2L': 23011544, '2R': 21146708, '3L': 24543557,
                  '3R': 27905053}

    return chr_length[chrom]


class Chromosome:
    """simple class to represent chromosome properties"""
    def __init__(self, chr_id):
        self.chr_id = chr_id
        self.chr_len = chromosome_length(self.chr_id)


class RunStruct(DictLike):
    """
    The RunStruct is a central structure that configures all of the steps in the
    pipeline required to generate maps of linked selection from prepared data
    files.
    """
    @staticmethod
    def _readable(attribute_key, attribute_value):
        """
        Return a file readable string "value_string" for attribute_value by
        modifying repr(attribute_value) such that when "value_string" is read
        from a file, eval("value_string") yields a Python object identical to
        the original attribute_value. Use regex to fix problematic cases.
        """
        # re to change numpy repr = 'array' to 'np.array'
        re_nparray = re.compile('array')
        # re to change numpy repr = 'inf' to 'np.inf'
        # re_inf = re.compile('\s-?inf[,\s]')
        re_inf = re.compile('[^A-Za-z]?inf[^A-Za-z]?')
        # re to change numpy 'nan' into 'np.nan'
        re_nan = re.compile('[^A-Za-z]?nan[^A-Za-z]?')
        # re to remove for list, dict, array line continuations
        re_whitespace = re.compile('[\s\n]+')
        # re to formation function names
        re_func = re.compile('<function (\w+) at [\d\w]+>')
        # re to match class instances
        re_class = re.compile('<(\w+\.)*(\w+) object at [\d\w]+>')
        # initial string representation
        value_string = repr(attribute_value)
        # substitute correct format for special strings
        if re_nparray.search(value_string):
            # format array value strings with high precision
            np.set_printoptions(precision=18)
            value_string = re_nparray.sub('np.array', value_string)
            np.set_printoptions()
        if re_inf.search(value_string):
            # save numpy inifinity with the np. prefix
            value_string = re.sub('inf', 'np.inf', value_string)
        if re_nan.search(value_string):
            # save numpy inifinity with the np. prefix
            value_string = re.sub('nan', 'np.nan', value_string)
        if re_whitespace.search(value_string):
            value_string = re_whitespace.sub(' ', value_string)
        # skip eval test on func & class: just return ID string
        if re_func.search(value_string):
            func_name = re_func.search(value_string).group(1)
            value_string = "'{}'".format(func_name)
            return value_string
        if re_class.search(value_string):
            class_name = re_class.search(value_string).group(2)
            value_string = "'{}'".format(class_name)
            return value_string

        # test the string encoding of each value
        try:
            assert eval(value_string) == attribute_value
        except ValueError:
            cmp_arr = izip(eval(value_string), attribute_value)
            try:
                assert all(np.allclose(v1, v2) for (v1, v2) in cmp_arr)
            except TypeError:
                msg = 'I don\'t know'
                print msg
        except AssertionError:
            msg = 'error saving {}\nnew: {}\nold: {}'
            print msg.format(attribute_key, value_string, str(attribute_value))
        except SyntaxError:
            msg = 'error saving {}\nnew: {}\nold: {}'
            print msg.format(attribute_key, value_string, str(attribute_value))
        except NameError:
            msg = 'error saving {}\nnew: {}\nold: {}'
            print msg.format(attribute_key, value_string, str(attribute_value))

        # return the readable string
        return value_string

    def __init__(self, init=None, clone=None, make=False,
                 root=root_dir,
                 gmap='AA_Map_10kb',
                 refg='hg19',
                 mask='strictMask',
                 ncon='euarchontoglires',
                 npct=0.35,
                 nval=0,
                 gmsk=0.0,
                 cons='primate',
                 outg='rheMac3',
                 nspc=4,
                 wind=5000,
                 slid=True,
                 neut='YRI',
                 dpop=None,
                 phs='phase3',
                 tkn='MISSING.TOKEN',
                 fcl='nonsyn',
                 bscl=100,
                 bs_annos=('primate_cons95_Segments',),
                 bdir='bmaps',
                 cs_annos=(),
                 cdir='cmaps',
                 bdfe=(default_fitness_effects,),
                 cdfe=(),
                 nff='neutral_filter_files_2',
                 chroms=human_autosomes,
                 methods=('Nelder-Mead',)):

        """
        OPTIONAL ARGUMENT FLAGS TO BYPASS DEFAULT INITIALIZATION
        --------------------------------------------------------
        :param init:optional initialization file for the run (bypasses remaining
        args and sets them from file)

        :param clone: existing RunStruct instance used as a template to set all
        attributes of new instance

        :param make: bypass normal __init__, go straight to _make function


        EXPLICIT INITIALIZATION FOR ARG-SPECIFIED ATTRIBUTES
        ----------------------------------------------------
        :param root: root directory for file paths
        :type root: str

        :param gmap: genetic map file name
        :type gmap: str

        :param refg: the id of the reference genome
        :type refg: str

        :param mask: the name of the call mask to use for polymorphism data
        :type mask: str

        :param ncon: the name of the alignment used with phastCons to call
        neutral sites where score=0
        :type ncon: str

        :param cons: the name of the alignment used with phastCons to call
        conserved regions
        :type cons: str

        :param outg: mutation rate estimate outgroup genome (e.g. 'rheMac3')
        :type outg: str

        :param neut: neutral polymorphism file name
        :type neut: str

        :param phs: the phase of the 1000 genomes data used
        :type phs: str

        :param tkn: token for the run and processed files
        :type tkn: str

        :param fcl: substitutions for collated plots & selective sweeps
        :type fcl: str

        :param bscl: the scale for background selection maps
        :type bscl: int

        :param bs_annos: tuple of annotations for bs files
        :type bs_annos: tuple

        :param bdir: path to write/read bmap files from
        :type: str

        :param cs_annos: tuple of annotations for cs files
        :type cs_annos: tuple

        :param cdir: path to write/read cmap files from
        :type: str

        :param bdfe: tuple of selection coefficients for bs annos
        :type bdfe: tuple

        :param cdfe: tuple of selection coefficients for bs annos
        :type cdfe: tuple

        :param nff: string name of neutral filter file used for neutmask
        :type nff: str

        :param chroms: tuple of chromosomes
        :type chroms: tuple

        :param methods: set of algorithms to use for optimization, defaults to
        all of them
        :type methods: tuple
        """
        super(RunStruct, self).__init__()

        """
        INITIALIZATION BYPASSES
        """
        # INIT FROM FILE IF FILE PATH PROVIDED
        if init:
            assert isinstance(init, str)
            self._init(init_file=init)
            return

        # INIT BY CLONING EXISTING INSTANCE
        if clone:
            assert isinstance(clone, RunStruct)
            self.setitems(items=clone.items, safe=False)
            return

        # INIT WITH PRE-FILLED ARGS
        if make:
            self._make()
            return

        """
        INITIALIZATION OF ARG-SPECIFIED ATTRIBUTES
        """
        self.root = root
        self.gmap = gmap
        self.refg = refg
        self.mask = mask
        self.ncon = ncon
        self.npct = npct
        self.nval = nval
        self.gmsk = gmsk
        self.cons = cons
        self.outg = outg
        self.nspc = nspc
        self.wind = wind
        self.slid = slid
        self.neut = neut
        self.dpop = dpop
        self.phs = phs
        self.tkn = tkn
        self.fcl = fcl
        self.bscl = bscl
        self.bs_annos = bs_annos
        self.bdir = bdir
        self.cs_annos = cs_annos
        self.cdir = cdir
        self.bdfe = bdfe
        self.cdfe = cdfe
        self.nff = nff
        self.chroms = chroms
        self.methods = methods
        self.msk_id = '0000'

        """
        HIGH-LEVEL FILE DIRECTORIES
        """
        self.data = root + '/data'
        self.prcl = root + '/precalc'
        self.cmpr = root + '/compress'
        self.rslt = root + '/result'

        """
        NUMERIC ATTRIBUTES
        """
        self.bnum = len(self.bs_annos)
        self.bsgrid = len(self.bdfe[0]) if self.bdfe else 0
        self.bsparams = sum(len(t) for t in self.bdfe)
        self.cnum = len(self.cs_annos)
        self.csgrid = len(self.cdfe[0]) if self.cdfe else 0
        self.csparams = sum(len(t) for t in self.cdfe)
        self.nchroms = len(self.chroms)

        """
        INTERNALLY SET ATTRIBUTES
        """
        # time stamp to affix to saved file before/after optimization, etc.
        self.timestamp = 'NOT_STARTED'
        # label for files including info about LS maps
        l_str = '{tkn}.BS{bnum}.{bsgrid}.CS{cnum}.{csgrid}'
        self.lbl = l_str.format(**self.dict)
        # file used to initialize the RunStruct -- only set by _init
        self.init_file = None
        # save file ptah
        t_str = '{rslt}/init_files/{neut}.{lbl}.{timestamp}.initial.txt'
        self.txt_file = t_str.format(**self.dict)
        # calc_ckgd executable
        self.calc_bkgd_exe = calc_bkgd_executable

        """
        INITIALIZATION OF RUNSTRUCT SUB-CLASSES
        """
        # file path templates for complete project file structure
        self.files = FileStruct(kwd=self.dict)

        # bounds, initial parameter values and fixed parameters
        self.fixed = FixedParams(bn=self.bsparams, cn=self.csparams)

        # stats calculated on data and inference results
        self.stat = Statistics(initial_params=self.fixed.initial_parameters)

        # misc. variables and flags
        self.vars = Variables()

        # minimization algorithms
        self.optimizers = []
        # note - more than 3 can be used, currently only using 2
        self.op1 = OptimizerParams()
        self.op2 = OptimizerParams()
        self.op3 = OptimizerParams()
        self._init_optimizers()

        """
        INITIALIZATION OF DERIVED ATTRIBUTES
        """
        # free parameter vector for bs, cs, tau: default initialization
        self.params = self.fixed.initial_parameters
        self.stat.init_params = self.params
        # only save attrs specified by __init__, ignore inherited attrs
        self.savekeys = self.keys
        # container used to store params during optimization
        self.cached_params = []

    """
    PRIVATE CLASS METHODS
    """
    def _init(self, init_file):
        """
        Load an init file and set attributes in safe mode skipping
        unrecognized keys
        """
        # create a default RunStruct to store data from init file
        initer = RunStruct()

        # read each line of the file and update attributes
        with open(init_file, 'r') as f:
            for line in f:
                # split (key, value) pairs on '='
                key, value = line.strip('-\n').split('=')
                # if key contains a period it is an embedded subclass
                if '.' in key:
                    cls, k = key.split('.')
                    # use safe attribute settings to build initializer
                    initer[cls].setitem(key=k, value=eval(value), safe=True)
                else:
                    initer.setitem(key=key, value=eval(value), safe=True)

        # set init_file and txt_file to the file used for initialization
        initer.init_file = initer.txt_file = init_file

        # copy attributes from initializer to self
        self.setitems(initer.items, safe=False)

    def _make(self):
        """make a new RunStruct with the defaults in this function"""
        self.__init__(tkn='pr95.default.settings')

    def _init_optimizers(self):
        """add optimizer class for methods; save attribute names in list"""
        bs_params = map(len, self.bdfe)
        cs_params = map(len, self.cdfe)
        umax = self.fixed.u_max
        for (i, method) in enumerate(self.methods, start=1):
            op = OptimizerParams(method=method, bs_params=bs_params,
                                 cs_params=cs_params, umax=umax)
            optimizer_id = 'op{}'.format(i)
            self[optimizer_id] = op
            self.optimizers.append(optimizer_id)

    """
    PUBLIC CLASS METHODS
    """
    def reset(self):
        """Re-initialize the RunStruct using current values for init args"""
        self.__init__(root=self.root, gmap=self.gmap, refg=self.refg,
                      mask=self.mask, ncon=self.ncon, cons=self.cons,
                      outg=self.outg, neut=self.neut, phs=self.phs,
                      tkn=self.tkn, fcl=self.fcl, bscl=self.bscl,
                      bs_annos=self.bs_annos, bdir=self.bdir,
                      cs_annos=self.cs_annos, cdir=self.cdir,
                      bdfe=self.bdfe, cdfe=self.cdfe, nff=self.nff,
                      chroms=self.chroms, methods=self.methods)

    def optimizer_inputs(self, func, args, jac=False):
        for oid in self.optimizers:
            self[oid].target_func = func
            self[oid].args = args
            # self[oid].jac = jac
            # self[oid].jac = True
            # self[oid] .hess = '3-point'

    def record_timestamp(self):
        """create a timestamp for the current time, reset txt_file suffix"""
        # generate YYMMDDHHMMSS-formatted timestamp from current time
        current_time = strftime('%y%m%d%H%M%S')
        # update timestamp on txt_file
        self.txt_file = self.txt_file.replace(self.timestamp, current_time)
        # update timestamp attribute
        self.timestamp = current_time
        # update the suffix of txt_file with timestamp
        # txt = '{root}/result/init_files/{neut}.{lbl}.{timestamp}.initial.txt'
        # self.txt_file = txt.format(**self.dict)

    def save(self, txt_file=None):
        """
        Save text file formatted to preserve all class attributes.

        Convert each member variable (key, value) pair to a string formatted
        --key=repr(value) and write to either internally defined or externally
        supplied txt_file. Modify the string format of certain special cases
        such as numpy objects, functions and classes.
        :param txt_file: optional new filename to replace default constructed
        name
        """
        # if new txt file name given, update internal name before saving
        if txt_file is not None:
            self.txt_file = txt_file

        f = open(self.txt_file, 'w')
        # save attributes in 'save_keys' only
        for key in self.savekeys:
            value = self[key]
            # expand subclasses and write out their attributes individually
            if isinstance(value, DictLike):
                for (k, v) in value.items:
                    k = '{}.{}'.format(key, k)
                    value_string = self._readable(k, v)
                    f.write('--{}={}\n'.format(k, value_string))
            # treat others as single objects
            else:
                value_string = self._readable(key, value)
                f.write('--{}={}\n'.format(key, value_string))

        f.close()

    def makeraw(self, make_init=False):
        """
        make_struct but filling each field manually
        :param make_init: optional init file closer to the final, new file
        to be created
        """
        # init from file or _make defaults
        if make_init:
            self.__init__(make=True)

        # use current root for saving by default
        savedir = self.root

        # possible root directories
        loc = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS'
        rem = '/ifs/data/c2b2/gs_lab/dam2214/pyLS'

        # help message prints on start
        print 'local_root={}\nremote_root={}\n' \
              '"finish" + ' \
              'enter completes init_file with defaults in remaining fields\n' \
              '"exit" + ' \
              'enter kills the program without saving init_file\n' \
              'note: string variables require quotes'.format(loc, rem)

        # confirm preset values for each key or take new value
        # keys = 'root gmap mask ncon cons outg neut phase token ' \
        #        'rm_annos bs_annos cs_annos bdfe cdfe methods'.split()
        keys = ['neut', 'tkn', 'bscl', 'bs_annos', 'bdir',
                'bdfe', 'methods']

        for k in keys:

            # enter the input for the key
            next_prompt = 'press enter if {}={} or enter an alternative: '
            nextinput = raw_input(next_prompt.format(k, self[k]))

            # up to 5 tries to rename the values for each key in keys
            tries = 0
            while True:

                # don't allow more than 5 tries
                if tries > 5:
                    max_tries = 'too many input errors -- {} set to {}'
                    print max_tries.format(k, self[k])
                    break

                # accept the default value and continue
                elif nextinput == '':
                    break

                # finish with defaults for the remaining fields
                elif nextinput == 'finish':
                    self.reset()
                    savetxt = self.txt_file.replace(self.root, savedir)
                    fin_cmmd = 'defaults for remaining fields; saving as:\n{}'
                    print fin_cmmd.format(savetxt)
                    self.save(txt_file=savetxt)
                    exit(1)

                # kill the program without saving
                elif nextinput == 'exit':
                    print 'exiting without saving an init file'
                    exit(1)

                # try the user input
                else:
                    try:
                        self[k] = eval(nextinput)  # try setting the value
                        break
                    except SyntaxError:
                        tries += 1  # count up the failed attempts
                        err = 'press enter if {}={} or enter an alternative: '
                        nextinput = raw_input(err.format(k, self[k]))
                    except NameError:
                        tries += 1
                        err = 'press enter if {}={} or enter an alternative: '
                        nextinput = raw_input(err.format(k, self[k]))

        # reset the RunStruct and save an init file to the current environment
        self.reset()
        savetxt = self.txt_file.replace(self.root, savedir)
        self.save(txt_file=savetxt)

        print 'init file saved as:\n{}'.format(savetxt)

    def cur_params(self, nlp):
        """record the current batch of optimization results"""
        # TODO: expand for cs annos
        # create list of params and mutation rates
        # mu_rates = [sum(10 ** self.params[i:j]) for (i, j) in self.bsidx]
        # px = [p for p in self.params] + [nlp] + mu_rates
        px = [p for p in self.params] + [nlp]

        # create a list of formatting
        idx = xrange(1, len(px) + 1)
        fmts = ['{:>17.10e} ' if i % 4 else '{:>17.10e}\n' for i in idx]

        # join char joins optimization cycle entries
        join_char = '\n\n' if len(fmts) % 4 else '\n'

        # create complete string-record
        current_params = ''.join(fmts).format(*px) + join_char

        # store current params as a string
        self.cached_params.append(current_params)

        # check whether or not there are enough records to write to file
        if self.stat.function_calls % 50 == 0:
            self.log_params()
        
        return current_params

    def log_params(self):
        """write current optimization progress to file, clear cache"""
        # write results to log file every 50 function calls
        with open(self.optimization_log, 'a') as f:
            # f.write(''.join(self.cached_params) + '\n')
            f.write(''.join(self.cached_params))            
        self.cached_params = []

    """
    DERIVED PARAM PROPERTIES
    """
    @property
    def bsidx(self):
        bidx = []
        for i in xrange(0, self.bsparams, self.bsgrid):
            bidx.append((i, i + self.bsgrid))
        return bidx

    @property
    def csidx(self):
        cidx = []
        for i in xrange(0, self.csparams, self.csgrid):
            i += self.bsparams
            cidx.append((i, i + self.csgrid))
        return cidx

    """
    RUN & LOG FILES
    """
    @property
    def final_file(self):
        """return path to final file with timestamp affixed"""
        new_dir = self.txt_file.replace('/init_files/', '/final_files/')
        new_suf = new_dir.replace('.initial.txt', '.final.txt')
        return new_suf

    @property
    def optimization_log(self):
        """path to optimization params log file"""
        fname = '{root}/result/logs/{neut}.{lbl}.{timestamp}.log'
        return fname.format(**self.dict)

    @property
    def data_indices(self):
        """segment indices corresponding to resampling windows"""
        return '{root}/arrays/boot/{lbl}.indices.npy'.format(**self.dict)

    @property
    def oid_mthd(self):
        return dict((o, m) for o, m in izip(self.optimizers, self.methods))

    @property
    def mthd_oid(self):
        return dict((m, o) for o, m in izip(self.optimizers, self.methods))

    """
    STATISTICS AND SUMMARIES
    """
    @property
    def uvec(self):
        """vector of deleterious mutation rates across DFE of each bs anno"""
        i, j = 0, self.bsparams
        dims = self.bnum, self.bsgrid
        return 10**self.params[i:j].reshape(dims)

    @property
    def utot(self):
        """total deleterious mutation rate for each bs anno"""
        return np.sum(self.uvec, axis=1)

    @property
    def upmf(self):
        """fraction of utot at per coefficient for each bs anno"""
        return self.uvec / self.utot

    @property
    def tmean(self):
        """mean selection coefficient per bs anno"""
        return np.average(self.bdfe, axis=0, weights=self.upmf)

    @property
    def avec(self):
        """vector of beneficial substitution rates across DFE of each cs anno"""
        i, j = self.bsparams, len(self.params)-1
        dims = self.cnum, self.csgrid
        return 10**self.params[i:j].reshape(dims)

    @property
    def atot(self):
        """total fraction beneficial substitution in each cs anno"""
        return np.sum(self.avec, axis=1)

    @property
    def apmf(self):
        """fraction of adaptive subs per coefficient for each cs anno"""
        return self.avec / self.atot

    @property
    def smean(self):
        """mean selection coefficient per cs anno"""
        return np.average(self.cdfe, axis=0, weights=self.upmf)


class ChromStruct(RunStruct):
    """RunStruct connected to a particular chromosome"""
    def __init__(self, chrom, **kwargs):
        super(ChromStruct, self).__init__(**kwargs)
        self.chrom = chrom
        if 'clone' in kwargs:
            return

    def reset(self):
        """Re-initialize the ChromStruct using current values for init args"""
        self.__init__(root=self.root, gmap=self.gmap, refg=self.refg,
                      mask=self.mask, ncon=self.ncon, cons=self.cons,
                      outg=self.outg, neut=self.neut, phs=self.phs,
                      tkn=self.tkn, fcl=self.fcl, bscl=self.bscl,
                      bs_annos=self.bs_annos, bdir=self.bdir,
                      cs_annos=self.cs_annos, cdir=self.cdir,
                      bdfe=self.bdfe, cdfe=self.cdfe, nff=self.nff,
                      chroms=self.chroms, methods=self.methods,
                      chrom=self.chrom)

    @property
    def chlen(self):
        """chrom length of current chrom"""
        return chromosome_length(self.chrom)

    @property
    def cidx(self):
        """0-based chrom index of current chrom"""
        return self.chroms.index(self.chrom)

    """
    FILE PATHS REQUIRING ADDITIONAL ARGS
    """
    def bs_target(self, anno):
        return self.files.fbtgt.format(ch=self.chrom, an=anno)

    def bkgd_file(self, anno, tval, pidx=-1, plen=0, merge=False):
        if (pidx >= 0) and (plen > 0):
            suf = 'p{}.{:.0e}.'.format(pidx, plen).replace('+', '')
        elif merge and (plen > 0):
            suf = 'merge{:.0e}.'.format(plen).replace('+', '')
        else:
            suf = ''
        return self.files.fbmap.format(ch=self.chrom, an=anno, t=tval, sf=suf)

    def bxct_file(self, anno, tval, suf):
        return self.files.fbxct.format(ch=self.chrom, an=anno, t=tval, sf=suf)

    def rand_file(self, n_sites):
        return self.files.fnrdm.format(ch=self.chrom, n=n_sites)

    def nt_simfiles(self, simulation_id):
        """simulated neutral [hom, het] pairs for a given set of params"""
        return self.files.fnsim.format(ch='{ch}', sim=simulation_id)

    def cs_target(self, anno):
        """CS target filepath for a given annotation and chrom"""
        return self.files.fctgt.format(ch=self.chrom, an=anno)

    def cmap_file(self, anno, sval, suf=''):
        return self.files.fcmap.format(ch=self.chrom, an=anno, s=sval, sf=suf)
    #
    # TODO: update this
    # def cxct_file(self, foc, sval, suf):
    #     f = self.prcl_dir + '/{cdir}/{chrom}.{a}.cxct.{u}.s{s:.8f}.npy'
    #     return f.format(a=foc, s=sval, u=suf, **self.dict)

    """
    GENOME DATA FILES
    """
    @property
    def gmap_files(self):
        """genetic map files for each chromosome"""
        return self.files.fgmap.format(ch=self.chrom)

    @property
    def refg_files(self):
        """reference genome FASTA files"""
        return self.files.frefg.format(ch=self.chrom)

    @property
    def ancs_files(self):
        """ancestor reference genome FASTA files"""
        return self.files.fancs.format(ch=self.chrom)

    @property
    def call_mask(self):
        """mask file for callability in polymorphism data set"""
        return self.files.fcall.format(ch=self.chrom)

    @property
    def axt_outgroup(self):
        """human:outgroup .axt-format alignment files for each chromosome"""
        return self.files.faxt.format(ch=self.chrom)

    @property
    def wigz_ncons(self):
        """phastCons scores files used to identify neutral segments"""
        return self.files.fnwig.format(ch=self.chrom)

    @property
    def wigz_cons(self):
        """wiggle format phastCons scores files for phylogenetic cons"""
        return self.files.fcwig.format(ch=self.chrom)

    @property
    def cons_dist(self):
        """file with distribution of PHASTCONS scores for cons alignment"""
        fmt = '{data}/cons/wigz/{cons}/{cons}.scores.dist.txt'
        return fmt.format(**self.dict)

    @property
    def snp_files(self):
        """SNP files for each chromosome (source of SNPs used in neut files)"""
        return self.files.fsnp.format(ch=self.chrom)

    @property
    def cpg_isld(self):
        """return CpG island file for chrom"""
        return self.files.fcpgi.format(ch=self.chrom)

    """
    PROCESSED DATA FILES
    """
    @property
    def neut_masks(self):
        """chromosome fasta mask. encoding: neutral=1; other=0"""
        return self.files.fnm.format(ch=self.chrom)

    @property
    def nr_exons(self):
        """non-redundant exons bed file"""
        return self.files.fnrex.format(ch=self.chrom)

    """
    COMPRESSED ARRAY FILES
    """
    @property
    def sg_files(self):
        """array of each segment length"""
        return self.files.fsg.format(ch=self.chrom)

    @property
    def bs_files(self):
        """array of bs predictions for each tval of each bs_anno"""
        return self.files.fbs.format(ch=self.chrom)

    @property
    def cs_files(self):
        """array of cs predictions for each sval of each cs_anno"""
        return self.files.fcs.format(ch=self.chrom)

    @property
    def nu_files(self):
        """array of neutral u estimates, with 'nan' segments interpolated"""
        return self.files.fnu.format(ch=self.chrom)

    @property
    def nt_files(self):
        """neutral [ref pairs, alt_pairs] summed across each segment"""
        return self.files.fnt.format(ch=self.chrom)

    @property
    def dv_files(self):
        """divergence stats for neutral data segments"""
        return self.files.fdv.format(ch=self.chrom)

    @property
    def pl_files(self):
        """array of polymorphic site counts per segment"""
        return self.files.fpl.format(ch=self.chrom)

    @property
    def dsmsk_file(self):
        """down sampling mask file"""
        return self.files.fmds.format(mskid=self.msk_id)

    @property
    def gc_files(self):
        """GC content for each segment"""
        return self.files.fgc.format(ch=self.chrom)

    @property
    def il_files(self):
        """GC content for each segment"""
        return self.files.fil.format(ch=self.chrom)

    @property
    def ai_files(self):
        """GC content for each segment"""
        return self.files.fai.format(ch=self.chrom)

    @property
    def cn_files(self):
        """GC content for each segment"""
        return self.files.fcn.format(ch=self.chrom)

    @property
    def cm_files(self):
        """GC content for each segment"""
        return self.files.fcm.format(ch=self.chrom)
