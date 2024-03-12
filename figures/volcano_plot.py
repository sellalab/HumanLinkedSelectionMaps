__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
from openpyxl import load_workbook
from figures.common_functions import format_panels
import matplotlib.pyplot as plt


#%% set paths -- use raw string for Dropbox path to deal with space and parens
db_path = r'/Users/davidmurphy/Dropbox (OMRF)/'  # r'string' = raw string format
# this assertion checks that the file path was set up correctly and exists
assert os.path.isdir(db_path)
# set the path to the excel file
xl_path = db_path + 'DESeq2_Lesion_vs_Normal_CancerAndNonCancerPatients.xlsx'
# this assertion checks that the file path is correct and file exists
assert os.path.isfile(xl_path)
# load workbook using and access "Sheet 1"
wb = load_workbook(xl_path)
ws = wb['Sheet 1']


#%% append relevant data for making the plot to lists

# FDR-adjusted pvalues
padj = []
for col in ws.columns[47]:
    padj.append(col.value)

# log2 fold enrichment
xfold = []
for col in ws.columns[43]:
    xfold.append(col.value)

# gene names
genes = []
for col in ws.columns[2]:
    genes.append(col.value)


#%% convert lists of numeric values to numpy arrays
padj = np.array(padj[1:])
xfold = np.array(xfold[1:])

# # rescale fold change to log10
# log10_xfold = np.log10(2.0**xfold)

# get -log10 of pvalues
log10_pval = -np.log10(padj)


#%% create array indices for 3 groups of points

# group 1: all values less than or equal to FDR-adjusted p-value of -log10(0.05)
# OR without log2-fold change above 1 or below -1
p_threshold = -np.log10(0.05)
idx_1 = (log10_pval <= p_threshold) | (abs(xfold) <= 1)

# group 2: all values with greater than FDR-adjusted p-value of -log10(0.05)
# AND with log2-fold expression less than -1
idx_2 = (log10_pval > p_threshold) & (xfold < -1)

# group 3: all values with greater than FDR-adjusted p-value of -log10(0.05)
# # AND with log2-fold expression greater than -1
idx_3 = (log10_pval > p_threshold) & (xfold > 1)

# group 4: all "significant" values for heatmap (i.e., idx_2 or idx_3)
sig_idx = (idx_2 | idx_3)

#%% create the volcano plot figure

plt.figure(figsize=(4,3))
# create a single set of axes on the figure
ax1 = plt.subplot(111)
# use my function for formatting the appearance of the axes
format_panels(ax1)
# ax1.grid(False)  # turn off grid lines
# adjust the space between the axes and the figure
plt.subplots_adjust(bottom=0.12, left=0.13, top=0.99, right=0.98)

# PLOTTING EACH SET OF POINTS:
# use idx_1 to select points corresponding to group 1 and scatter plot
plt.scatter(xfold[idx_1], log10_pval[idx_1], s=2, alpha=0.5, color='gray',
            label='not significant')
# use idx_2 for points in group 2, colored blue
plt.scatter(xfold[idx_2], log10_pval[idx_2], s=2, alpha=0.5, color='dodgerblue',
            label='down')
# use idx_3 for points in group 3, colored red
plt.scatter(xfold[idx_3], log10_pval[idx_3], s=2, alpha=0.5, color='red',
            label='up')

# PLOTTING THRESHOLD LINES:
# set values for threshold line width, style and color
thresh_lw = 0.75
thresh_ls = ':'
thresh_col = 'black'
# plot the log10 FDR-adjusted p-value threshold as a horizontal (axhline) line
plt.axhline(p_threshold, color=thresh_col, lw=thresh_lw, ls=thresh_ls)
# plot the two fold-change thresholds as verticle (axvline) lines
plt.axvline(-1, color=thresh_col, lw=thresh_lw, ls=thresh_ls)
plt.axvline(1, color=thresh_col, lw=thresh_lw, ls=thresh_ls)

# now label the axes
plt.ylabel('-log10 FDR-adjusted p-value', labelpad=2)
plt.yticks(x=0.01)  # (this adjusts how close the ticks are to the axis)
plt.xlabel('log2 fold-change', labelpad=2)
xtck = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
plt.xticks(xtck, y=0.02)

# add the plot legend for the 3 labeled groups of points
plt.legend(ncol=1, frameon=1, framealpha=1, facecolor='white',
           borderaxespad=0.5, borderpad=0.5, handletextpad=0.1,
           labelspacing=0.25, prop=dict(size=9))


# set name for the saved file and save
fig_name = 'DESeq2_Lesion_vs_Normal_CancerAndNonCancerPatients_VolcanoPlot.png'
fig_path = db_path + fig_name
plt.savefig(fig_path, dpi=512)  # (dpi = dots oer inch - 512 is high res)
plt.close()


#%% create heatmap figure

tab_file = xl_path.replace('.xlsx', '.txt')
data_arr = np.loadtxt(tab_file, skiprows=1, usecols=(range(6, 42)))
d_use = data_arr[sig_idx]
d_adj = np.log10(np.maximum(1, d_use))
# heatmap = seaborn.clustermap(d_adj)
# fig_name = 'DESeq2_Lesion_vs_Normal_CancerAndNonCancerPatients_Heatmap.png'
# fig_path = db_path + fig_name
# heatmap.savefig(fig_path, dpi=512)
#%%
plt.figure(figsize=(4, 3))
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.99)
ax = plt.subplot(111)
format_panels(ax)
plt.hist(d_adj.flatten(), bins=100, log=True)
plt.ylabel('log10 number of gene/sample pairs')
plt.xlabel('log10 read counts')
plt.savefig(db_path + 'ReadCountHist.png', dpi=512)
plt.close()
# heatmap.

#%%
sidx = np.argsort(d_adj.flatten())
