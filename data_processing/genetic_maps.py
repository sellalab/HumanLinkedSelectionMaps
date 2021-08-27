__author__ = 'davidmurphy'

from classes.geneticmap import *
from collections import namedtuple

# paths for genetic map files
maps_dir = root_dir + '/data/maps'
anjali_template = maps_dir + '/anjali/maps_{ch}.txt'
# decode_file = maps_dir + '/aau1043_DataS3'
decode_folder = maps_dir + '/deCODE_2019'
# create folder if it doesn't yet exist
if not os.path.isdir(decode_folder):
    os.mkdir(decode_folder)

# make a named tuple for data lines in deCODE file
DataLine = namedtuple('DataLine', 'ch start end cmmb cm')


#%% split the lifted over deCODE map file into chromosomes
def split_decode_file():
    """read the deCODE file and split into individual files per chrom"""
    # split files by chromosome
    header = []
    current_chrom = 'chr1'
    # file_template = decode_folder + '/{}.deCODE_2019.GRCh38.txt'
    file_template = decode_folder + '/{}.deCODE_2019_hg19.txt'
    decode_file = decode_folder + '/aau1043_DataS3_hg19_liftOver.bed'
    w = open(file_template.format(current_chrom), 'a')
    print('NOTE: appending to map files, not overwriting. may cause duplicates')
    with open(decode_file, 'r') as f:
        for line in f:
            # save the header info
            if line.startswith('#'):
                header.append(line)
            # save the column labels
            elif line.startswith('Chr'):
                header.append('# ' +  line)
                # write header to first file now
                w.write(''.join(header))
            # the remaining lines are data
            else:
                # get the chromosome for the current line
                ch = line.split()[0]
                # if the chromosome matches the open file, write to it
                if ch == current_chrom:
                    w.write(line)
                # if a new chromosome arises, switch to a new writefile
                else:
                    w.close()
                    current_chrom = ch
                    w = open(file_template.format(current_chrom), 'a')
                    # write header to file
                    w.write(''.join(header))
                    w.write(line)

        # close the last open file
        w.close()


split_decode_file()


#%% set path to the liftover file into hg19 coords
import matplotlib.pyplot as plt
chrom = 'chr9'
lifted_file = decode_folder + '/{}.deCODE_2019_hg19.txt'.format(chrom)
# get coordinates and genetic map values
gmap = np.loadtxt(lifted_file, usecols=(1,2,3,4))

#%% sort data by start positions
sidx = np.argsort(gmap[:,3])

#%% iterate over rows of the liftover map and fix any wierd gaps
last_start = 1  # initial pos
last_gpos = 0  # initial gpos
skipped_row = []
kept_row = []
n_skipped = 0
n_kept = 0
for row in gmap[sidx]:
    start, end, oldcmmb, gpos = row
    # calculate distance based on the last position that was kept
    dist = end - last_start
    # skip segments where the distances get mixed up
    if dist <= 0:
        skipped_row.append(row)
        n_skipped += 1  # count number of skipped rows
        continue
    # calculate genetic distance from last gpos
    gdist = gpos - last_gpos
    assert gdist >= 0
    # calculate recombination rate
    cmmb = 1e6 * gdist / dist
    # make a new row including the old cmmb
    new_row = last_start, end, cmmb, oldcmmb, gpos
    kept_row.append(new_row)
    n_kept += 1  # count number of kept rows
    # update last start and last gpos
    last_start = end
    last_gpos = gpos

new_gmap = np.array(kept_row)
print('kept: {}, skipped: {}'.format(n_kept, n_skipped))

#%% plot and compare maps with and without editing
plt.figure(figsize=(15,15))
plt.title('genetic map', fontsize=24)
plt.plot(new_gmap[:,0], new_gmap[:,4])
plt.show()

plt.figure(figsize=(15,15))
plt.title('cMMb', fontsize=24)
plt.plot(new_gmap[:,0], new_gmap[:,2], label='new')
plt.plot(new_gmap[:,0], new_gmap[:,3], label='old')
plt.legend()
plt.show()

#%% plot histogram of errors to the map
plt.figure(figsize=(15,15))
# plt.hist(np.log10(np.maximum(0.01,abs(new_gmap[:,2]-new_gmap[:,3]))), bins=500)
plt.hist(abs(new_gmap[:,2]-new_gmap[:,3]), bins=500, log=1)

plt.show()




#%% make the YRI LD map
for c in xrange(1,23):
    chrom = 'chr{}'.format(c)
    f_in = anjali_template.format(ch=chrom)
    anjali2gmap(chrom, f_in, 'YRI_LD', 1e4)


#%% make deCODE 2019 maps lifted into hg19 and formatted the standard way
def create_decode_maps(chrom, savefile=False):
    # set file path for hg19 lifted map files
    lifted_file = decode_folder + '/{}.deCODE_2019_hg19.txt'.format(chrom)
    # get coordinates and genetic map values
    gmap = np.loadtxt(lifted_file, usecols=(1, 2, 3, 4))
    # get sorting index using GENETIC MAP POS
    sidx = np.argsort(gmap[:, 3])
    # # separate map into start, end, cmmb, gpos and sort by gpos
    # start, end, cmmb, gpos = gmap.T

    # create the new map
    last_start = 1  # initial pos
    last_gpos = 0  # initial gpos
    skipped_row = []  # record skipped rows
    kept_row = []  # record kept rows
    n_skipped = 0  # count number of skipped
    n_kept = 0  # count number of kept
    for row in gmap[sidx]:
        start, end, oldcmmb, gpos = row

        # calculate distance based on the last position that was kept
        dist = end - last_start
        # skip segments where the distances get mixed up
        if dist <= 0:
            skipped_row.append(row)
            n_skipped += 1  # count number of skipped rows
            continue

        # calculate genetic distance from last gpos
        gdist = gpos - last_gpos
        assert gdist >= 0

        # calculate recombination rate (cM/Mb)
        cmmb = 1e6 * gdist / dist

        # make a new row including the old cmmb
        new_row = [last_start, cmmb, gpos]
        kept_row.append(new_row)
        n_kept += 1  # count number of kept rows

        # update last start and last gpos
        last_start = end
        last_gpos = gpos

    # at the end of the map, add a new row filling to the end of the chrom
    chlen = chromosome_length(chrom)
    last_row = [chlen, 0, last_gpos]
    kept_row.append(last_row)

    # reformat map in a single array
    new_gmap = np.array(kept_row)
    n_total = 0.01 * (n_kept + n_skipped)
    msg_args = (chrom, n_kept, n_kept/n_total, n_skipped, n_skipped/n_total)
    print('{}: kept={} ({:.0f}%), skipped={} ({:.0f}%)'.format(*msg_args))

    # optional map save
    if savefile:
        header = 'deCODE 2019 lifted from GRCh38 to hg19\npos cM/Mb cM'
        f_save = decode_folder + '/{}_deCODE_2019.txt'.format(chrom)
        np.savetxt(f_save, new_gmap, fmt='%d %.10f %.6f', header=header)


for c in xrange(1, 23):
    ch = 'chr{}'.format(c)
    create_decode_maps(ch, True)
#%%