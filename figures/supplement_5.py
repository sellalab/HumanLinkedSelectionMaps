__author__ = 'davidmurphy'

import numpy as np
import matplotlib.pyplot as plt
from figures.common_functions import format_panels
from classes.runstruct import root_dir, cst_from_fldr
from figures.sub_rate_figs import neutcons_uratio_single


final_dir = root_dir + '/result/final_files'

"""UDEL SUMMARY PLOT SUPPLEMENT SECTION 5"""


#%% GET DATA FOR SUPPLEMENTAL FIGURE 32

# upper and lower bounds of u total
ulo = 1.29
uhi = 1.51

# data from runs and phylofit analysis
pct = [96, 95, 94, 93, 92]
an_all = ['ape', 'primate', 'prosimian', 'euarchontoglires',
          'laurasiatheria', 'mammal', 'fish']
lb_all = ['ape', 'prim', 'pros', 'supr', 'laur', 'mamm', 'vert']
# s_ape = [1 - neutcons_uratio_single('ape', pc=p)[0].mean() for p in pct]
s_fis = [1 - neutcons_uratio_single('fish', pc=p)[0].mean() for p in pct]
s_cad = [1 - neutcons_uratio_single('cadd', pc=p)[0].mean() for p in pct]
s_all = [1 - neutcons_uratio_single(a, pc=94)[0].mean() for a in an_all]

# get inferred udel for ape and CADD
u_fis = []
u_cad = []
for p in pct:
    if p == 94:
        # cadfldr = 'cadd94_gmask_mnb_378'
        cadfldr = 'cadd94_gmask_v1.6_without_bstat'
        # fisfldr = 'fish_cons94_gmask_mnb_378'
        fisfldr = 'fish_cons94_new'
    else:
        # cadfldr = 'cadd{}_gmask'.format(p)
        cadfldr = 'cadd{}_gmask_v1.6_without_bstat'.format(p)
        # fisfldr = 'fish_cons{}_gmask'.format(p)
        fisfldr = 'fish_cons{}_new'.format(p)
    u_cad.append(cst_from_fldr(cadfldr).stat.utot[0] * 1e8)
    u_fis.append(cst_from_fldr(fisfldr).stat.utot[0] * 1e8)

# for fldr in ['ape_cons{}_clean'.format(p) for p in pct]:
#     cst = cst_from_fldr(fldr)
#     udel = cst.stat.utot[0] * 1e8
#     u_fis.append(udel)
# for fldr in ['cadd{}'.format(p) for p in pct]:
#     cst = cst_from_fldr(fldr)
#     udel = cst.stat.utot[0] * 1e8
#     u_cad.append(udel)


# SUPPLEMENTAL FIGURE 30
fig = plt.figure(figsize=(6.5, 5.4))
plt.subplots_adjust(wspace=0.05, hspace=0.2, left=0.1, right=0.995, top=0.95,
                    bottom=0.07)
x1 = np.arange(len(pct))
txt_y = 0.35
row1ytxt = 0.75
# PANEL A: ape result
ax1 = plt.subplot(231)
format_panels(ax1)
plt.text(-0.4, row1ytxt, 'a', ha='center', fontweight='bold')
plt.bar(x1, s_fis, color='darkorange')
for (x, p) in zip(x1, pct):
    plt.text(x, txt_y, 'vert. {}%'.format(int(100 - p)), color='white',
             rotation=90, va='top', ha='center')
plt.xticks(x1, ['' for _ in x1])
ylab = 'proportion of deleterious\nmutations'
plt.ylabel('proportion of\ndeleterious mutations', labelpad=2)
plt.ylim(0, 0.8)

# PANEL B: CADD result
ax2 = plt.subplot(232)
format_panels(ax2)
plt.text(-0.4, row1ytxt, 'b', ha='center', fontweight='bold')
plt.bar(x1, s_cad, color='firebrick')
for (x, p) in zip(x1, pct):
    plt.text(x, txt_y, 'CADD {}%'.format(int(100 - p)), color='white',
             rotation=90, va='top', ha='center')
plt.xticks(x1, ['' for _ in x1])
# plt.xlabel('putatively selected sites', labelpad=0)
plt.yticks(color='none')
plt.ylim(0, 0.8)

# PANEL C: all annotations
x2 = np.arange(len(s_all))
ax3 = plt.subplot(233)
format_panels(ax3)
plt.text(-0.4, row1ytxt, 'c', ha='center', fontweight='bold')
plt.bar(x2, s_all, color='dodgerblue')
for (x, l) in zip(x2, lb_all):
    plt.text(x, txt_y, '{} 6%'.format(l), color='white', rotation=90,
             va='top', ha='center')
plt.xticks(x2, ['' for _ in x2])
plt.yticks(color='none')
plt.ylim(0, 0.8)

# common title upper panels
fig.add_subplot(211, frame_on=False)
# plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.yticks([])
plt.xticks([])
plt.title('estimates based on evolutionary rates')
ymax = 1.25
txt_x = -0.2
txt_y = 1.16

# PANEL D: inferred vs. divergence udel range
# ape section
ax4 = plt.subplot(234)
format_panels(ax4)
utop = np.array([u / ulo for u in u_fis])
ubot = np.array([u / uhi for u in u_fis])
uheight = utop-ubot
# print utop
# print ubot
lbl = 'background selection'
plt.bar(x1, uheight, bottom=ubot, color='gray', alpha=0.75, label=lbl)
lbl = 'evolutionary rates'
plt.plot(x1, s_fis, 'o', color='k', ms=6, label=lbl)
plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
xtck = ['{}%'.format(int(100 - p)) for p in pct]
plt.xticks(x1, xtck, rotation=45, y=0.05)
plt.xlabel('phastCons', labelpad=2)
plt.xlim(-0.45, len(pct)-0.55)
plt.ylabel('proportion of\ndeleterious mutations', labelpad=2)
plt.yticks(x=0.045)
plt.ylim(0, ymax)
plt.text(txt_x, txt_y, 'd', ha='center', fontweight='bold')


# plt.text(-0.4, 2.1, 'D', ha='center')
# x1 = np.arange(len(pct))
# ubot = [s * ulo for s in s_ape]
# utop = [s * uhi - s * ulo for s in s_ape]
# # lbl = r'estimated $\mu_{del}$ range'
# lbl = 'evolutionary rates'
# plt.bar(x1, utop, bottom=ubot, color='darkorange', alpha=0.75, label=lbl)
# # lbl = r'inferred $\mu_{del}$'
# lbl = 'background selection'
# plt.plot(x1, u_fis, 'o', color='k', ms=6, label=lbl)
# xtck = ['{}%'.format(int(100 - p)) for p in pct]
# plt.xticks(x1, xtck, rotation=90)
# plt.xlabel('ape')
# plt.ylabel(r'deleterious mutation rate ($\times 10^{-8}$)')
# plt.legend(loc='upper right', frameon=1, framealpha=0.5, facecolor='white')
# plt.ylim(0, 2.3)

# CADD section
ax5 = plt.subplot(235)
format_panels(ax5)
utop = np.array([u / ulo for u in u_cad])
ubot = np.array([u / uhi for u in u_cad])
uheight = utop-ubot
lbl = 'background selection'
plt.bar(x1, uheight, bottom=ubot, color='gray', alpha=0.75, label=lbl)
lbl = 'evolutionary rates'
plt.plot(x1, s_cad, 'o', color='k', ms=6, label=lbl)
plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
xtck = ['{}%'.format(int(100 - p)) for p in pct]
plt.xticks(x1, xtck, rotation=45, y=0.05)
plt.xlabel('CADD', labelpad=2)
plt.xlim(-0.45, len(pct)-0.55)
plt.yticks(color='none')
plt.ylim(0, ymax)
plt.text(txt_x, txt_y, 'e', ha='center', fontweight='bold')
leg = plt.legend(loc='lower center', frameon=1, framealpha=0.75,
                     facecolor='white', handlelength=1,
                     handletextpad=0.4)
leg.set_title('estimates from:', prop={'size': 10})
leg._legend_box.align = 'left'

# plt.text(-0.4, 2.1, 'E', ha='center')
# ubot = [s * ulo for s in s_cad]
# utop = [s * uhi - s * ulo for s in s_cad]
# # lbl = r'estimated $\mu_{del}$ range'
# lbl = 'evolutionary rates'
# plt.bar(x1, utop, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
# # lbl = r'inferred $\mu_{del}$'
# lbl = 'background selection'
# plt.plot(x1, u_cad, 'o', color='k', ms=6, label=lbl)
# xtck = ['{}%'.format(int(100 - p)) for p in pct]
# plt.xticks(x1, xtck, rotation=90)
# plt.xlabel('CADD')
# plt.yticks(color='none')
# plt.ylim(0, 2.3)
# plt.legend(loc='upper right', frameon=1, framealpha=0.5, facecolor='white')

# PANEL F: table
ax6 = plt.subplot(236, frame_on=False)
format_panels(ax6)
plt.title('f', fontweight='bold')
plt.xticks([])
plt.yticks([])
# columns = ['annotation', r'BS ($\times10^{-8}$)',
#            r'ER ($\times10^{-8}$)']
columns = ['annotation', 'ER', 'BS']
data = []
# add ape data
for i in xrange(len(pct)):
    lbl = 'vert. {}%'.format(int(100 - pct[i]))
    # u = '{:.2f}'.format(u_fis[i])
    # ur = '[{:.2f}, {:.2f}]'.format(s_ape[i] * ulo, s_ape[i] * uhi)
    u = '{:.2f}'.format(s_fis[i])
    ur = '[{:.2f}, {:.2f}]'.format(u_fis[i] / uhi, u_fis[i] / ulo)
    data.append((lbl, u, ur))
# add CADD data
for i in xrange(len(pct)):
    lbl = 'CADD {}%'.format(int(100 - pct[i]))
    # u = '{:.2f}'.format(u_cad[i])
    # ur = '[{:.2f}, {:.2f}]'.format(s_cad[i] * ulo, s_cad[i] * uhi)
    u = '{:.2f}'.format(s_cad[i])
    ur = ' [{:.2f}, {:.2f}]'.format(u_cad[i] / uhi, u_cad[i] / ulo)
    data.append((lbl, u, ur))

the_table = plt.table(cellText=data, colLabels=columns,cellLoc='center',
                      colColours=['lightgray'] * 3, bbox=[0, 0, 1, 1])
the_table.auto_set_font_size(False)
the_table.set_fontsize(8)

# common title lower panels
fig.add_subplot(2, 3, (4,5), frame_on=False)
# plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.yticks([])
plt.xticks([])
plt.title('comparing the two kinds of estimates')

# # common X-axis
# fig.add_subplot(211, frame_on=False)
# plt.tick_params(labelcolor="none", bottom=False, left=False)
# plt.yticks([])
# plt.xticks([])
# plt.xlabel('putatively selected sites', fontsize=16, labelpad=25)
sdir = root_dir + '/result/final_files/sfigs'
f_save = sdir + '/fig_S30.udel_estimates.png'
plt.savefig(f_save, dpi=512)
plt.close()


#%% CONSTANT PER GAMETE RATE SUPPLEMENT FIG X

# get ape data
ape_file = root_dir + '/data/bsanno/ape_site_counts.txt'
ape_cnt = {}
with open(ape_file, 'r') as f:
    for line in f:
        pct, cnt = map(int, line.split())
        ape_cnt[pct] = float(cnt)

ape_udel = {}
for pct in xrange(90, 99):
    fldr = 'ape_cons{}_clean'.format(pct)
    ape_udel[pct] = cst_from_fldr(fldr).stat.utot[0]

# get CADD data
cadd_file = root_dir + '/data/bsanno/cadd_site_counts.txt'
cadd_cnt = {}
with open(cadd_file, 'r') as f:
    for line in f:
        pct, cnt = map(int, line.split())
        cadd_cnt[pct] = float(cnt)

cadd_udel = {}
for pct in xrange(90, 99):
    fldr = 'cadd{}'.format(pct)
    cadd_udel[pct] = cst_from_fldr(fldr).stat.utot[0]
#
pct_range = range(98, 89, -1)
xi = np.arange(9)
y_ape = np.array([ape_cnt[p]*ape_udel[p] for p in pct_range])
y_cadd = np.array([cadd_cnt[p]*cadd_udel[p] for p in pct_range])
plt.figure(figsize=(3, 2.5))
plt.subplots_adjust(top=1, bottom=0.18, left=0.22, right=1)
plt.bar(xi-0.2, y_ape, width=0.4, color='darkorange', alpha=0.8,
        label='conserved in apes')
plt.bar(xi+0.2, y_cadd, width=0.4, color='dodgerblue', alpha=0.8, label='CADD')
plt.yticks(x=0.03)
plt.ylabel('deleterious mutations per gamete\nper generation')
plt.ylim(0, 2.2)
plt.xticks(xi, ['{}%'.format(100-p) for p in pct_range], rotation=45, y=0.04)
plt.xlabel('threshold', labelpad=2)
plt.legend(loc='upper left')
f_save = final_dir + '/mainfigs/fig_4.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%% FIGURE 4 IN MAIN GET DATA

#  upper and lower bounds of u total
ulo = 1.29
uhi = 1.51

# data from runs and phylofit analysis
# pct = [96, 95, 94, 93, 92, 91]
pct = [94, 93, 92]

s_ape = [1 - neutcons_uratio_single('ape', pc=p)[0].mean() for p in pct]
s_cad = [1 - neutcons_uratio_single('cadd', pc=p)[0].mean() for p in pct]

# get inferred udel for ape and CADD
u_ape = []
for fldr in ['ape_cons{}_clean'.format(p) for p in pct]:
    cst = cst_from_fldr(fldr)
    udel = cst.stat.utot[0] * 1e8
    u_ape.append(udel)
u_cad = []
for fldr in ['cadd{}'.format(p) for p in pct]:
    cst = cst_from_fldr(fldr)
    udel = cst.stat.utot[0] * 1e8
    u_cad.append(udel)

#%% FIGURE 4 IN MAIN (MODIFIED FROM SUPPLEMENTAL FIGURE)

fig = plt.figure(figsize=(5, 3))
plt.subplots_adjust(wspace=0.1, left=0.12, right=1, top=0.94, bottom=0.17)

# panel A: ape rates
plt.subplot(121)
x1 = np.arange(len(pct))
utop = np.array([u / ulo for u in u_ape])
ubot = np.array([u / uhi for u in u_ape])
uheight = utop-ubot
print utop
print ubot
lbl = 'background selection'
plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
lbl = 'evolutionary rates'
plt.plot(x1, s_ape, 'o', color='k', ms=6, label=lbl)
plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
xtck = ['{}%'.format(int(100 - p)) for p in pct]
plt.xticks(x1, xtck, rotation=45, y=0.03)
plt.xlim(-0.45, len(pct)-0.55)
plt.xlabel('conservation threshold')
plt.ylabel('proportion of deleterious mutations')
# plt.legend(loc='upper right', frameon=1, framealpha=0.5, facecolor='white')
plt.ylim(0, 1.1)

# panel B: CADD rates
plt.subplot(122)
# ubot = [s * ulo for s in s_cad]
# utop = [s * uhi - s * ulo for s in s_cad]
utop = np.array([u / ulo for u in u_cad])
ubot = np.array([u / uhi for u in u_cad])
uheight = utop-ubot
print utop
print ubot
lbl = 'background selection'
plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
lbl = 'evolutionary rates'
plt.plot(x1, s_cad, 'o', color='k', ms=6, label=lbl)
plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
xtck = ['{}%'.format(int(100 - p)) for p in pct]
plt.xticks(x1, xtck, rotation=45, y=0.03)
plt.xlim(-0.45, len(pct)-0.55)
plt.xlabel('CADD threshold')
plt.yticks(color='none')
plt.ylim(0, 1.1)
plt.legend(title='estimates from:', loc='lower left', frameon=1,
           framealpha=0.5, facecolor='white')
plt.text(0.32, 0.95, 'A', transform=plt.gcf().transFigure)
plt.text(0.78, 0.95, 'B', transform=plt.gcf().transFigure)

sdir = root_dir + '/result/final_files/mainfigs'
f_save = sdir + '/fig_4.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%%