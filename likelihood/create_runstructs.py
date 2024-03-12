__author__ = 'davidmurphy'


import numpy as np
from classes.runstruct import root_dir, RunStruct, human_autosomes

gm = 0.1
nv = 1
sl = False
wn = 6000
ns = 8
pct = 94
npct = 35
cons = 'ape'
ncons = 'euarchontoglires'
spec = 'ape fish primate prosimian euarchontoglires laurasiatheria mammal'
cell_lines = ['GM12878', 'H1hesc', 'Huvec', 'Helas3', 'Hepg2', 'K562']
dfe_6pt = np.power(10.0, np.arange(-4.5, -1.5, 0.5))
dfe_4pt = np.power(10.0, np.arange(-3.5, -1.5, 0.5))
dpop = None


#%% USE ALTERNATE GENETIC MAPS
alt_gmap = 'deCODE_2019'
# alt_gmap = 'YRI_LD'

bs_annos = ['fish_cons94_gmask', 'cadd94_gmask']
for b_an in bs_annos:
    b_an_1 = b_an
    bdir = b_an
    neut = 'YRI'
    tk = b_an + '_' + alt_gmap
    a = RunStruct(gmap=alt_gmap,
                  nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()


#%% BS ACROSS POPS
pfile = root_dir + '/data/snps/grouped_pops.txt'
# create list of pops ordered by group and pop:group dict
bs_annos = ['fish_cons94_gmask', 'cadd94_gmask']
with open(pfile, 'r') as f:
    for line in f:
        pp, gp = line[:-1].split()
        if pp == 'YRI':
            continue
        for b_an in bs_annos:
            b_an_1 = b_an
            bdir = b_an
            neut = pp
            tk = b_an
            a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                          npct=npct / 100.0, nval=nv, gmsk=gm,
                          cons=cons, neut=neut, dpop=dpop,
                          bs_annos=(b_an_1,),
                          bdir=bdir,
                          bdfe=(dfe_6pt,),
                          cs_annos=(),
                          cdfe=(),
                          )
            a.fixed.cth = None
            a.save()


#%% BS EXONIC/NONEXONIC
b_an_1 = 'fish_cons94_gmask_exonic'
b_an_2 = 'fish_cons94_gmask_nonexonic'
bdir = 'fish_cons94_gmask_exnex'
tk = bdir
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% BS EXONIC/NONEXONIC INDIVIDUALLY (FISH)
b_an_1 = 'fish_cons94_gmask_exonic'
b_an_2 = 'fish_cons94_gmask_nonexonic'
for b_an in [b_an_1, b_an_2]:
    bdir = 'fish_cons94_gmask_exnex'
    tk = b_an
    neut = 'YRI'
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  # cdir=cdir,
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()


#%% GENIC ANNOS
b_an_1 = 'cds'
b_an_2 = 'peri'
bdir = 'genic'
tk = bdir + '_gmask'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% CHROMHMM ANNOS
b_an_1 = 'huvec_txn'
b_an_2 = 'huvec_pro'
b_an_3 = 'huvec_enh'
bdir = 'huvec_chromhmm'
tk = bdir + '_gmask'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% BS + nonsyn CS
b_an_1 = 'fish_cons94_gmask'
bdir = b_an_1
c_an_1 = 'YRI_nonsyn_s1'
cdir = c_an_1
tk = b_an_1 + '_' + c_an_1
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, ),
              bdir=bdir,
              bdfe=(dfe_4pt, ),
              cs_annos=(c_an_1,),
              cdir=cdir,
              cdfe=(dfe_4pt,),
              )
a.fixed.cth = None
a.save()


#%% BS + NS/OTHER CONS
b_an_1 = 'fish_cons94_gmask'
bdir = b_an_1
for pc in range(91, 99):
    c_an_1 = 'nonsyn_fish_cons{}_gmask_hc_subs'.format(pc)
    c_an_2 = 'other_fish_cons{}_gmask_hc_subs'.format(pc)
    cdir = 'fish_cons{}_gmask_hc_subs'.format(pc)
    tk = b_an_1 + '_' + cdir
    neut = 'YRI'
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1, ),
                  bdir=bdir,
                  bdfe=(dfe_4pt, ),
                  cs_annos=(c_an_1, c_an_2),
                  cdir=cdir,
                  cdfe=(dfe_4pt, dfe_4pt),
                  )
    a.fixed.cth = None
    a.save()


#%% BS + NS/OTHER CADD
b_an_1 = 'cadd94_gmask'
bdir = b_an_1
for pc in range(91, 99):
    c_an_1 = 'nonsyn_cadd{}_gmask_hc_subs'.format(pc)
    c_an_2 = 'other_cadd{}_gmask_hc_subs'.format(pc)
    cdir = 'cadd{}_gmask_hc_subs'.format(pc)
    tk = b_an_1 + '_' + cdir
    neut = 'YRI'
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1, ),
                  bdir=bdir,
                  bdfe=(dfe_4pt, ),
                  cs_annos=(c_an_1, c_an_2),
                  cdir=cdir,
                  cdfe=(dfe_4pt, dfe_4pt),
                  )
    a.fixed.cth = None
    a.save()


#%% CS ONLY (NONSYN/OTHER CADD)
for pc in range(91, 99):
    c_an_1 = 'nonsyn_cadd{}_gmask_hc_subs'.format(pc)
    c_an_2 = 'other_cadd{}_gmask_hc_subs'.format(pc)
    cdir = 'cadd{}_gmask_hc_subs'.format(pc)
    tk = cdir
    neut = 'YRI'
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(),
                  bdir=bdir,
                  bdfe=(),
                  cs_annos=(c_an_1, c_an_2),
                  cdir=cdir,
                  cdfe=(dfe_4pt, dfe_4pt),
                  )
    a.fixed.cth = None
    a.save()


#%% CS ONLY (NONSYN)
c_an_1 = 'YRI_nonsyn_s1'
cdir = c_an_1
tk = c_an_1
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(),
              bdir=bdir,
              bdfe=(),
              cs_annos=(c_an_1,),
              cdir=cdir,
              cdfe=(dfe_4pt,),
              )
a.fixed.cth = None
a.save()


#%% CS ONLY (NONSYN/OTHER CONS)
for pc in range(91, 99):
    c_an_1 = 'nonsyn_fish_cons{}_gmask_hc_subs'.format(pc)
    c_an_2 = 'other_fish_cons{}_gmask_hc_subs'.format(pc)
    cdir = 'fish_cons{}_gmask_hc_subs'.format(pc)
    tk = cdir
    neut = 'YRI'
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(),
                  bdir=bdir,
                  bdfe=(),
                  cs_annos=(c_an_1, c_an_2),
                  cdir=cdir,
                  cdfe=(dfe_4pt, dfe_4pt),
                  )
    a.fixed.cth = None
    a.save()


#%% BS EXONIC/NONEXONIC (CADD)
b_an_1 = 'cadd94_gmask_exonic'
b_an_2 = 'cadd94_gmask_nonexonic'
bdir = 'cadd94_gmask_exnex'
tk = bdir
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% BS EXONIC/NONEXONIC INDIVIDUALLY (CADD)
b_an_1 = 'cadd94_gmask_exonic'
b_an_2 = 'cadd94_gmask_nonexonic'
for b_an in [b_an_1, b_an_2]:
    bdir = 'cadd94_gmask_exnex'
    tk = b_an
    neut = 'YRI'
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  # cdir=cdir,
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()


#%% CG MUTATIONS FILTERED
b_an_1 = 'cadd94_gmask'
bdir = b_an_1
tk = b_an_1 + '.filter.CG.mut'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1,),
              bdir=bdir,
              bdfe=(dfe_6pt,),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% CHROMHMM MERGED INCLUDING MORE STATES
b_an_1 = 'txn_merged'
b_an_2 = 'pro_merged'
b_an_3 = 'enh_merged'
b_an_4 = 'ins_merged'
b_an_5 = 'rep_merged'
bdir = 'chromhmm_merged'
tk = bdir
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3, b_an_4, b_an_5),
              bdir=bdir,
              bdfe=(dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% CHROMHMM MERGED WITH EXON SUBSTITUTED FOR TRANSCRIPTION STATE
b_an_1 = 'exon'
b_an_2 = 'pro_merged'
b_an_3 = 'enh_merged'
b_an_4 = 'ins_merged'
b_an_5 = 'rep_merged'
bdir = 'chromhmm_merged'
tk = 'chromhmm_plus_exon'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3, b_an_4, b_an_5),
              bdir=bdir,
              bdfe=(dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% cCREs MERGED WITH CDS ANNOTATIONS, 4PTS
b_an_1 = 'cds'
# b_an_2 = 'peri'
b_an_2 = 'cCRE_PLS_filtered'
b_an_3 = 'cCRE_ELS_filtered'
b_an_4 = 'cCRE_CTCF_filtered'
bdir = 'cCRE_filtered'
tk = 'genic_plus_cCRE'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3, b_an_4),
              bdir=bdir,
              bdfe=(dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% cCREs MERGED WITH CDS ANNOTATIONS, 6PTS
b_an_1 = 'cds'
# b_an_2 = 'peri'
b_an_2 = 'cCRE_PLS_filtered'
b_an_3 = 'cCRE_ELS_filtered'
b_an_4 = 'cCRE_CTCF_filtered'
b_an_5 = 'cCRE_H3K4me3_filtered'
bdir = 'cCRE_filtered'
tk = 'cds_plus_all_cCRE'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3, b_an_4, b_an_5),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt),
              # bdfe=(dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt),

              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()

#%% cCREs MERGED WITH CDS ANNOTATIONS, 4PTS
b_an_1 = 'cds'
# b_an_2 = 'peri'
b_an_2 = 'PLS_union'
b_an_3 = 'ELS_union'
b_an_4 = 'CTCF_union'
b_an_5 = 'H3K4me3_union'
bdir = 'cCRE_filtered'
tk = 'cds_plus_all_cCRE'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3, b_an_4, b_an_5),
              bdir=bdir,
              # bdfe=(dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt),
              bdfe=(dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% cCREs SPLIT INTO ABOVE/BELOW MEDIAN CELL OVERLAPS
b_an_1 = 'cds'
# b_an_2 = 'peri'
b_an_2 = 'PLS_above_median'
b_an_3 = 'PLS_below_median'
b_an_4 = 'ELS_above_median'
b_an_5 = 'ELS_below_median'
b_an_6 = 'H3K4me3_above_median'
b_an_7 = 'H3K4me3_below_median'
b_an_8 = 'CTCF_above_median'
b_an_9 = 'CTCF_below_median'
bdir = 'cCRE_split'
tk = 'cds_plus_split_cCRE'
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2, b_an_3, b_an_4, b_an_5, b_an_6, b_an_7,
                        b_an_8, b_an_9),
              bdir=bdir,
              # bdfe=(dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt, dfe_4pt,
              #       dfe_4pt, dfe_4pt, dfe_4pt),
              bdfe=(dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt, dfe_6pt,
                    dfe_6pt, dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% NEW CADD ANNOTATIONS WITH AND WITHOUT MCVICKER B (CADD v1.6)
bs_annos = ['cadd94_gmask_v1.6_without_bstat', 'cadd94_gmask_v1.6']

for b_an_1 in bs_annos:
    bdir = b_an_1
    neut = 'YRI'
    tk = b_an_1
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()


#%% NEW CADD ANNOTATIONS WITHOUT MCVICKER B (CADD v1.6) ALL PCT
bs_annos = ['cadd94_gmask_v1.6_without_bstat', 'cadd94_gmask_v1.6']
for pct in [91, 92, 93, 95, 96, 97, 98]:
    b_an_1 = 'cadd{}_gmask_v1.6_without_bstat'.format(pct)
    bdir = b_an_1
    neut = 'YRI'
    tk = b_an_1
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()


#%% NEW CADD ANNOTATIONS WITHOUT MCVICKER B (CADD v1.6) EXNEX
b_an_1 = 'cadd94_gmask_v1.6_without_bstat_exonic'.format(pct)
b_an_2 = 'cadd94_gmask_v1.6_without_bstat_nonexonic'.format(pct)
bdir = 'cadd94_gmask_v1.6_without_bstat_exnex'
neut = 'YRI'
tk = bdir
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt),
              cs_annos=(),
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% ALL REPLACEMENT RUNS WITH NEW CADD SCORES

# store file names in single text file for job processing
file_names = []

# 1. new CADD with all pops
pfile = root_dir + '/data/snps/grouped_pops.txt'
b_an_1 = 'cadd94_gmask_v1.6_without_bstat'
bdir = b_an_1
tk = b_an_1
with open(pfile, 'r') as f:
    for line in f:
        pp, gp = line[:-1].split()
        if pp == 'YRI':
            continue
        neut = pp
        a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                      npct=npct / 100.0, nval=nv, gmsk=gm,
                      cons=cons, neut=neut, dpop=dpop,
                      bs_annos=(b_an_1,),
                      bdir=bdir,
                      bdfe=(dfe_6pt,),
                      cs_annos=(),
                      cdfe=(),
                      )
        a.fixed.cth = None
        a.save()
        file_names.append(a.txt_file)

# 2. use different choices phylo depth for neutral sites
neut = 'YRI'
for sp in spec.split():
    if sp == 'euarchontoglires':
        continue
    tk = '{}_{}'.format(b_an_1, sp)
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=sp,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()
    file_names.append(a.txt_file)

# 3. use different cutoff nvals
for ncutoff in [0, 1, 2, 3, 4, 5, 10, 20, 30]:
    tk = '{}_cutoff_nval_{}'.format(b_an_1, ncutoff)
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=ncutoff, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()
    file_names.append(a.txt_file)

# 4. use different genetic mask cutoff
for gcutoff in [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1]:
    tk = '{}_cutoff_gmask_{}'.format(b_an_1, gcutoff)
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gcutoff,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()
    file_names.append(a.txt_file)

# 5. use different window size for neutral substitution rate variation
for wsize in range(4000, 12000, 1000):
    tk = '{}_neutsubwindow_{}'.format(b_an_1, wsize)
    a = RunStruct(nspc=ns, wind=wsize, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()
    file_names.append(a.txt_file)


#%% 6. filter CG
b_an_1 = 'cadd94_gmask_v1.6_without_bstat'
bdir = b_an_1
tk = b_an_1 + '_filter.CG.mut'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1,),
              bdir=bdir,
              bdfe=(dfe_6pt,),
              cs_annos=(),
              cdfe=(),
              )
a.fixed.cth = None
a.save()
file_names.append(a.txt_file)
# OTHERS THAT CAN BE DONE WITHOUT CHANGING INIT FILES:
# -- change B thresholds


#%% FIXED PHASTCONS SCORES
# 2. use different choices phylo depth for neutral sites
neut = 'YRI'
for sp in spec.split():
    b_an_1 = '{}_cons94_new'.format(sp)
    bdir = b_an_1
    tk = b_an_1
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()
    # file_names.append(a.txt_file)


#%% BS EXONIC/NONEXONIC
b_an_1 = 'fish_cons94_new_exonic'
b_an_2 = 'fish_cons94_new_nonexonic'
bdir = 'fish_cons94_new_exnex'
tk = bdir
neut = 'YRI'
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1, b_an_2),
              bdir=bdir,
              bdfe=(dfe_6pt, dfe_6pt),
              cs_annos=(),
              # cdir=cdir,
              cdfe=(),
              )
a.fixed.cth = None
a.save()


#%% BS EXONIC/NONEXONIC INDIVIDUAL
b_an_1 = 'fish_cons94_new_exonic'
b_an_2 = 'fish_cons94_new_nonexonic'
bdir = 'fish_cons94_new_exnex'
neut = 'YRI'
for b_an in [b_an_1, b_an_2]:
    tk = b_an
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  # cdir=cdir,
                  cdfe=(),
                  )
    a.fixed.cth = None
    a.save()


#%% DROP EACH CHROM FROM DATASET

for c in range(1, 23):
    drop_chrom = 'chr{}'.format(c)
    b_an_1 = 'cadd94_gmask_v1.6_without_bstat'
    bdir = b_an_1
    neut = 'YRI'
    tk = b_an_1
    a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
                  npct=npct / 100.0, nval=nv, gmsk=gm,
                  cons=cons, neut=neut, dpop=dpop,
                  bs_annos=(b_an_1,),
                  bdir=bdir,
                  bdfe=(dfe_6pt,),
                  cs_annos=(),
                  cdfe=(),
                  )
    a.vars.drop_chrom = drop_chrom
    a.fixed.cth = None
    # a.txt_files = a.txt_file.replace(tk, tk+'_drop_{}'.format(drop_chrom))
    a.save(txt_file=a.txt_file.replace(tk, tk+'_drop_{}'.format(drop_chrom)))

#%%
