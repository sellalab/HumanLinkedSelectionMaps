__author__ = 'davidmurphy'


import numpy as np
from classes.runstruct import root_dir, RunStruct, human_autosomes

# parameter for the percentage of genetic map to chop of from the end of chroms
gm = 0.1
# parameter for the maximum phastcons value to keep for neutral sites (0-1000 scale, 0=least conserved)
nv = 1
# use sliding windows (not used in any paper/appendix results)
sl = False
# size of the window for neutral substitution rate estimates
wn = 6000
# number of species used in the phylo tree for neutral substitution rate estimation
ns = 8
# percentile of conservation (e.g., 94 = "top 6% conserved")
pct = 94
# bottom percentile of sites to keep for neutral data (deprecated but still used in file names)
npct = 35
# another deprecated label for the phastcons tree used for calling conservation
cons = 'fish'
# the name of the phylo tree used for calling neutral sites
ncons = 'euarchontoglires'
# the distribution of selection effects with a six-point grid
dfe_6pt = np.power(10.0, np.arange(-4.5, -1.5, 0.5))
# used for population divergence (not used in any of the paper/appendix results)
dpop = None
# population to use for neutral polymorphism
neut = 'YRI'
# the background selection map to use
b_an_1 = 'fish_cons94_new'
# the directory where bmaps are stored (same as the token used for the map names)
bdir = b_an_1
# token used for saved run files (in this case we're using the ID of the bmap as token)
tk = b_an_1

# make a new RunStruct object using the settings above
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons,
              npct=npct / 100.0, nval=nv, gmsk=gm,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(b_an_1,),
              bdir=bdir,
              bdfe=(dfe_6pt,),
              cs_annos=(),
              cdfe=(),
              )

a.files.fnm = '/Users/MURPHYD/Dropbox (OMRF)/linked_selection/lsm_run/data/nmsk/{ch}.euarchontoglires.0.35.MaskLowRecomb.nmsk.npz'
# turn of one of the experimental threshold params (not used in any paper/appendix results)
a.fixed.cth = None
# save the RunStruct data structure as a text file that can be read to initalize a RunStruct later
a.save()
