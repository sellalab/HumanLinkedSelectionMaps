import os
import gzip


# for c in {73..934}; do cat wigs/n_${c}.laurasiatheria.chr4.wig >> zips/laurasiatheria.chr4.wig.gz

zipped = gzip.open("zips/laurasiatheria.chr4.wig.gz", "a")

for n in xrange(73, 935):
    cur = "wigs/n_{c}.laurasiatheria.chr4.wig".format(c=n)
    if os.path.isfile(cur):
        with open(cur, "r") as f:
            zipped.write(f.read())
#    os.remove(cur)

    
