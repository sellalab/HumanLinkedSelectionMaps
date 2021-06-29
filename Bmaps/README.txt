B maps
------

The following maps are formatted [B, length], 
after McVicker et al., 2009, where B is the 
B value (scaled from 0-1000) covering a given 
chromosome segment and length is the length 
of the segment in base pairs. Segments sum to
the length of each chromosome in the hg19 
build.

For example, the row:
 
	717 29307

assigns a B value of 0.717 (after rescaling)
to a segment of 29307 base pairs.

The two best-fitting models presented in our
paper are saved under phastCons_bestfit and
CADD_bestfit. These correspond to models 
using the 6% most constrained sites from 
phastCons conservation in 99-vertebrates
and constraint measured in CADD v1.6 (excluding
McVicker's B statistic). They both use a 
threshold applied in lookup tables of B=0.6.