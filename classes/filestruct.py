__author__ = 'davidmurphy'


from classes.dictlike import DictLike


class FileStruct(DictLike):
    """A class to organizes the file system used by the program"""
    def __init__(self, kwd):
        """

        :param kwd:
        """
        # initialize the parent class (DictLike)
        super(FileStruct, self).__init__()

        # # dictionary taken from RunStruct provides field values for strings
        # self.kwd = kwd

        """
        REFERENCE DATA FILES
        """
        # gmap files
        self.fgmap = '{data}/maps/{gmap}/{{ch}}_{gmap}.txt'.format(**kwd)
        # reference genome for current build
        self.frefg = '{data}/fa/ref/{{ch}}.fa.gz'.format(**kwd)
        # reconstructed HC ancestor genome (EPO pipeline)
        self.fancs = '{data}/fa/anc/hc/{{ch}}.hc.fa.gz'.format(**kwd)
        # sequencing call mask for 1000G data used
        self.fcall = '{data}/fa/mask/{{ch}}.{mask}.fa.gz'.format(**kwd)
        # axt formatted HM alignment file (for neutral u estimate)
        f_str = '{data}/div/axt/{outg}/{{ch}}.{refg}.{outg}.net.axt.gz'
        self.faxt = f_str.format(**kwd)
        # 100-vertebrate PHASTCONS scores (used to call neutral sites)
        f_str = '{data}/cons/wigz/{ncon}/{{ch}}.{ncon}.wig.gz'
        self.fnwig = f_str.format(**kwd)
        # 8-primate PHASTCONS scores (used to call conserved sites)
        self.fcwig = self.fnwig.replace(kwd['ncon'], kwd['cons'])
        # SNP data from 108 YRI individuals, 1000G, phs 3
        f_str = '{data}/snps/frq/{neut}/{{ch}}.{neut}.all.{phs}.frq.count.gz'
        self.fsnp = f_str.format(**kwd)
        # CpG island coordinates in BED file
        f_str = '{data}/coords/cpg_islands/{{ch}}.cpg.islands.bed'
        self.fcpgi = f_str.format(**kwd)
        # UCSC knownGene file for hg19
        self.fgene = '{data}/ch_features/knownGene.txt'.format(**kwd)
        # UCSC knownGene non-redundant exon BED file
        nr_prfx = '{data}/coords/nr/exon/{{ch}}'.format(**kwd)
        token = 'knownGene_nonredundant_exonSegments'
        self.fnrex = '{}_{}.bed'.format(nr_prfx, token)
        # UCSC knownGene non-redundant cds BED file
        nr_prfx = '{data}/coords/nr/cds/{{ch}}'.format(**kwd)
        token = 'knownGene_nonredundant_cdsSegments'
        self.fnrcds = '{}_{}.bed'.format(nr_prfx, token)
        # UCSC knownGene non-redundant gene file
        nr_prfx = '{data}/coords/nr/gene/{{ch}}'.format(**kwd)
        token = 'knownGene_coding_nonredundant_genes'
        self.fnrgene = '{}_{}.txt'.format(nr_prfx, token)
        # TODO: nr intron, etc.

        """
        PRECALCULATED LINKED SELECTION MAP FILES
        """
        f_str = '{prcl}/{bdir}/{gmap}.{{ch}}.{{an}}.t{{t:.8f}}.{{sf}}bkgd'
        self.fbmap = f_str.format(**kwd)
        f_str = '{prcl}/{cdir}/{gmap}.{{ch}}.{{an}}.s{{s:.8f}}.{{sf}}npz'
        self.fcmap = f_str.format(**kwd)

        """
        CONSTRUCTED DATA FILES
        """
        # TODO: use compressed bed files
        # annotations for regions subject to background selection
        self.fbtgt = '{data}/bsanno/{{an}}/{{ch}}.{{an}}.bed'.format(**kwd)
        # annotations for substitutions subject to classic sweeps
        self.fctgt = '{data}/csanno/{{an}}/{{ch}}.{{an}}.npz'.format(**kwd)
        # text file with list of coordinate files to mask neutral sites
        self.fnff = '{data}/coords/{nff}.txt'.format(**kwd)
        # neutrality, divergence and polymorphism status for every site
        if 'CpG' in kwd['tkn']:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.filter.CpG.nmsk.npz'
        elif 'BGC' in kwd['tkn']:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.filter.BGC.nmsk.npz'
        elif 'aligned.8' in kwd['tkn']:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.aligned.8.nmsk.npz'
        elif 'CG.mut' in kwd['tkn']:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.filter.CG.mut.nmsk.npz'
        elif 'gmap.edge' in kwd['tkn']:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.gmsk{gmsk}.nmsk.npz'
        elif 'cutoff_nval' in kwd['tkn']:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.nval{nval}.nmsk.npz'
        else:
            f_str = '{data}/nmsk/{{ch}}.{ncon}.{npct}.nmsk.npz'
        # f_str = '{data}/nmsk/{{ch}}.{outg}.{neut}.{tkn}.nmsk.npz'
        self.fnm = f_str.format(**kwd)
        # path to substitution count files used for estimating neutral urate
        #     sdir = rdir + 'aligned_{}_win_{}'.format(cst.nspc, cst.wind)
        #         fsub = '{}/{}.subcount.filled.txt'.format(sdir, cst.chrom)
        self.fsub = '{data}/phast/aligned_{nspc}_win_{wind}/{{ch}}.subcount.filled.txt'.format(**kwd)

        """
        PRECALCULATED LS MAPS & NEUTRAL DATA COMPRESSED IN SEGMENT GRID
        """
        # segment lengths after merging maps
        self.fsg = '{cmpr}/sarr/{{ch}}.{lbl}.sarr.npz'.format(**kwd)
        # (little) b values for each anno & tval in each segment
        self.fbs = '{cmpr}/barr/{{ch}}.{lbl}.barr.npz'.format(**kwd)
        # fixation rates due to CS for each focal & sval in each segment
        self.fcs = '{cmpr}/carr/{{ch}}.{lbl}.carr.npz'.format(**kwd)
        # estimated neutral substitution rate in each segment
        self.fnu = '{cmpr}/uarr/{{ch}}.{lbl}.{outg}.uarr.npz'.format(**kwd)
        # hom & het pairs for each segment
        self.fnt = '{cmpr}/narr/{{ch}}.{lbl}.{neut}.narr.npz'.format(**kwd)
        # GC content for each segment
        self.fgc = '{cmpr}/gc_arr/{{ch}}.{tkn}.gc.npz'.format(**kwd)
        # CpG island content for each segment
        self.fil = '{cmpr}/il_arr/{{ch}}.{tkn}.il.npz'.format(**kwd)
        # archaic introgression for each segment
        self.fai = '{cmpr}/ai_arr/{{ch}}.{tkn}.ai.npz'.format(**kwd)
        # mean cons for each segment
        self.fcn = '{cmpr}/cn_arr/{{ch}}.{tkn}.cn.npz'.format(**kwd)
        # mean cM/Mb for each segment
        self.fcm = '{cmpr}/cm_arr/{{ch}}.{tkn}.cm.npz'.format(**kwd)
        # # segment data for each segment
        # self.fsa = '{cmpr}/sg_arr/{{ch}}.{tkn}.sg.npz'.format(**kwd)
        # *simulated* hom & het pairs data for each segment
        self.fnsim = '{cmpr}/narr/{{ch}}.{lbl}.{{sim}}.narr.npz'.format(**kwd)
        # total neutral & total diverged neutral sites in each segment
        self.fdv = '{cmpr}/darr/{{ch}}.{lbl}.{outg}.darr.npz'.format(**kwd)
        # neutral polymorphic sites in each segment
        self.fpl = '{cmpr}/parr/{{ch}}.{lbl}.{neut}.parr.npz'.format(**kwd)

        """
        CALIBRATION AND TESTING FILES
        """
        # exact b calculation for a sample of sites (edge & random)
        # f_str = '{prcl}/{bdir}/{{ch}}.{{an}}.t{{t:.8f}}.bxct.{{sf}}.npz'
        f_str = '{prcl}/{bdir}/{{ch}}.{{an}}.t{{t:.8f}}.bxct.{{sf}}.npy'
        self.fbxct = f_str.format(**kwd)
        # random sample of neutral sites from the neutral mask file
        f_str = '{data}/nmsk/rdm/{{ch}}.{outg}.{tkn}.nrdm.n{{n}}.npz'
        self.fnrdm = f_str.format(**kwd)
        # exact cs calculation for a sample of sites
        f_str = '{prcl}/{cdir}/{{ch}}.{{an}}.s{{s:.8f}}.cxct.{{sf}}.npz'
        self.fcxct = f_str.format(**kwd)

        """
        DOWNSAMPLING AND MASKING FILES
        """
        self.fmds = '{cmpr}/marr/dnsmpl.mask.{{mskid}}.npz'.format(**kwd)