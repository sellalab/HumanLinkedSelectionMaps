import numpy as np


def py2mat(parr):
    # rearrange python array in matlab order (TAU + CS + BS)
    pidx = np.array([-1] + range(20, 35) + range(20))
    parr = parr[pidx]
    parr[1:16] = 10 ** parr[1:16]
    return parr


def matfree(marr):
    # get the free params from matlab results
    ifree = np.array(
        [3, 11, 12, 13, 14, 15, 23, 24, 25, 26, 27, 35, 36, 37, 38, 39,
         59, 60, 61, 62, 63, 71, 72, 73, 74, 75, 83, 84, 85, 86, 87, 95,
         96, 97, 98, 99]) - 1
    return marr[ifree]


def mat2py(marr):
    # rearrange matlab array in python order (BS + CS + TAU)
    pidx = np.array(range(16, 36) + range(1, 16) + [0])
    marr = marr[pidx]
    marr[20:35] = np.log10(marr[20:35])
    return marr


def main():
    # bkgd layout: ('exon', 'longintron', 'intergenic', 'UTR')
    # cs layout: ('exonicNS', 'intronic', 'UTR')

    mat_fix_cmmb = np.array([
        0, 0, 3.725463e+00, 0, 0, 0, 0, 0, 0, 0, 0, 4.975797e-03, 4.122033e-02,
        7.291548e-05, 3.034992e-01, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        3.513778e-07, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 5.028729e-02, 5.488467e-01,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -7.762214e+00,
        -8.540214e+00, -8.473901e+00, -9.936860e+00, -1.000000e+01, -10, -10,
        -10, -10, -10, -10, -8.167491e+00, -8.688135e+00, -9.999512e+00,
        -1.000000e+01, -1.000000e+01, -1.000000e+01, -10, -10, -10, -10, -10,
        -10, -8.167491e+00, -9.999477e+00, -9.999000e+00, -1.000000e+01,
        -1.000000e+01, -1.000000e+01, -10, -10, -10, -10, -10, -10,
        -8.167491e+00, -7.756361e+00, -8.212550e+00, -9.423483e+00,
        -9.121251e+00, -9.999006e+00, -10, -10, -10, -10, -10, -10,
        -8.167491e+00, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
        -8.167491e+00
    ])
    # organized with python layout
    mat_supp_table = np.array([2.73, 0.01, 0.75, 0, 0,  # exon
                               0.55, 0, 0, 0.01, 0,  # longintron
                               0.02, 0, 0, 0, 0,  # intergenic
                               6.95, 0.04, 0, 0.06, 0.01,  # UTR
                               0, 0.011, 0.02, 0.038, 0.365,  # exonicNS
                               0, 0, 0, 0, 0,  # intronic
                               0, 0, 0, 0.068, 0.345,  # UTRs
                               3.03531287])

    b_exp = np.log10(np.maximum(1e-10, mat_supp_table[:20]*3.5e-09))
    s_exp = np.log10(np.maximum(1e-10, mat_supp_table[20:35]))
    # s_exp[s_exp < -10] = 0
    mat_supp_table[:20] = b_exp
    mat_supp_table[20:35] = s_exp

    # py1 = np.array([3.86436435140000, 0.000775047000000000, 0.00455653000000000,
    #                 0.0305542000000000, 6.02149000000000e-17,
    #                 1.28470000000000e-18,
    #                 2.72029000000000e-18, 1.41557000000000e-18,
    #                 4.66678000000000e-14,
    #                 1.40588000000000e-17, 4.23656000000000e-16,
    #                 0.000215419000000000,
    #                 3.79920000000000e-23, 1.39129000000000e-13,
    #                 3.54934000000000e-15,
    #                 3.88940000000000e-17, -7.80850718350000, -9.58600667660000,
    #                 -22.7246212090000, -8.94443180110000, -21.9800642870000,
    #                 -19.8254901170000, -23.5656288320000, -23.8849747300000,
    #                 -23.8437202190000, -23.8660578030000, -26.6275585560000,
    #                 -20.6634246360000, -25.0186784800000, -23.0279087190000,
    #                 -24.0326092530000, -7.20712689760000, -21.8974961300000,
    #                 -22.0187250150000, -23.3194728250000, -24.5473981820000])

    # py_fixed_theta = np.array([4.17214624019399416, 7.57826460e-04, 4.28538765e-03, 2.87628680e-02,
    #      1.01815891e-14, 6.85354636e-16, 3.41208778e-20,
    #      8.23355351e-21, 2.21970608e-17, 7.26977240e-18,
    #      1.67929455e-19, 1.18734735e-04, 4.06491031e-23,
    #      1.90021229e-27, 6.12412261e-17, 1.15660394e-14,
    #      -7.911629938059484779, -9.218737208946818029,
    #      -23.469121828937151264, -8.962722549199678923,
    #      -20.058197112827539144, -53.454585365834695665,
    #      -23.920582859109309481, -25.584043787664420222,
    #      -26.58702587652468452, -26.061694492070706985,
    #      -20.726744117483530516, -21.123716499452388717,
    #      -22.861690930214162876, -23.720013593921933648,
    #      -26.269870553540549452, -7.21003809480269009,
    #      -10.779476946536936666, -19.250774910267725915,
    #      -10.054618088421159783, -19.358931375245653328])
    py_fix_cmmb = np.array([-7.656161782616624834, -8.635290990242886267,
                            -8.380658364905242408, -9.781538089444927309,
                            -18.871307648531157497, -8.847801662039149306,
                            -20.698348121985986836, -21.033498154921883128,
                            -23.815817569466481984, -24.713269391343452241,
                            -19.44554513613423552, -18.967247854107284866,
                            -27.969886094236176888, -25.167564050534906528,
                            -24.084509890739351334, -7.681030311951144718,
                            -16.601433520889390394, -25.75209228693328356,
                            -23.924267736386255478, -10.536303580801250845,
                            -13.115111284915768053, -2.451847625290505039,
                            -1.396210071989311574, -2.086325151687896451,
                            -0.454834299564158606, -16.261811675656126397,
                            -16.408653176667080231, -17.574888743825855641,
                            -22.360956729160857748, -35.438055318296370899,
                            -18.072048081503300665, -64.255547513596013687,
                            -3.121105048609699395, -1.098682996967399239,
                            -0.387952922055858351, 3.605367037574128908])
    # mat = matfree(mat_fix_cmmb)
    # py = np.maximum(-10, py_fix_cmmb)
    # prm2 = py2mat(py)
    # py = py2mat(mat_supp_table)

    prms = [
        py2mat(mat_supp_table),
        matfree(mat_fix_cmmb),
        py2mat(np.maximum(-10, py_fix_cmmb))
    ]

    labs = [
        'Elyashiv Supplement',
        'Independent MATLAB run',
        'Independent Python run'
    ]

    cols = 'darkorange darkturquoise fuchsia'.split()
    num = len(prms)

    import matplotlib.pyplot as plt
    import seaborn

    # CS params

    fig = plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=1, top=0.88, hspace=0.5)
    plt.suptitle('CS parameters compared', fontsize=16)

    j = np.arange(5) + 1
    coefs = np.arange(-1.5, -6, -1)
    map_names = ['exonicNS', 'intronic', 'UTR']
    for i in xrange(3):
        ax = plt.subplot(2, 2, i + 1)
        x = np.arange(5)
        width = 0.8 / num
        shift = -0.4
        all_nan = True
        for (p, c, l) in zip(prms, cols, labs):
            px = p[j]
            px[px < 0.01] = np.nan  # censor tiny values for CS
            if np.isfinite(px).any():
                all_nan = False
            plt.bar(x+shift, px, width, color=c, label=l)
            shift += width

        if not all_nan:
            plt.title(map_names[i], fontsize=14)
            plt.xticks(x, [r'$\mathrm{10^{%.1f}}$' % s for s in coefs], fontsize=12)
            plt.xlabel('selection coefficient', fontsize=14)
            plt.yticks(fontsize=14)
            plt.ylabel('fraction beneficial', fontsize=14)
            # plt.legend(loc='best')
        else:
            plt.title(map_names[i] + ' (N/A)', fontsize=14)
            plt.xticks([])
            plt.yticks([])
        j += 5
        if i == 0:
            plt.legend(loc='best')
        # if all_nan:
        #     fig.delaxes(ax)

    plt.subplot(224)
    u_mean = 0.14761379280275985
    pi_mean = 0.012652782640468048
    for (i, (p, c, l)) in enumerate(zip(prms, cols, labs)):
        tau = p[0]
        pi0 = u_mean / tau
        reduction = 1 - pi_mean / pi0
        plt.bar(i, reduction, 0.8, color=c, label=l)

    plt.title('diversity reduction', fontsize=14)
    plt.yticks(fontsize=14)
    # plt.ylabel(r'$\mathrm{1-\vfrac{\bar{\pi}}{\pi_0}}$', fontsize=16)
    plt.ylabel(r'$\mathrm{1-(\bar{\pi}/\pi_0})$', fontsize=16)
    plt.xticks([])

    # BS params
    fig = plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=1, top=0.88, hspace=0.5)
    plt.suptitle('BS parameters compared', fontsize=16)

    j = np.arange(5) + 16
    coefs = np.arange(-1.5, -6, -1)
    map_names = ['exon', 'longintron', 'intergenic', 'UTR']
    for i in xrange(4):
        ax = plt.subplot(2, 2, i + 1)
        x = np.arange(5)
        width = 0.8 / num
        shift = -0.4
        all_nan = True
        for (p, c, l) in zip(prms, cols, labs):
            px = (10 ** p[j]) / 3.5e-09
            px[px < 0.1] = np.nan
            if np.isfinite(px).any():
                all_nan = False
            plt.bar(x + shift, px, width, color=c, label=l)
            shift += width

        if not all_nan:
            plt.title(map_names[i], fontsize=14)
            plt.xticks(x, [r'$\mathrm{10^{%.1f}}$' % s for s in coefs], fontsize=12)
            plt.xlabel('selection coefficient', fontsize=14)
            plt.yticks(fontsize=14)
            # plt.ylabel('deleterious u', fontsize=14)
            plt.ylabel(r'$\mathrm{u_{del}/\mu}$', fontsize=14)
        else:
            plt.title(map_names[i] + ' (N/A)', fontsize=14)
            plt.xticks([])
            plt.yticks([])

        j += 5
        if i == 0:
            plt.legend(loc='best')
        # if all_nan:
        #     fig.delaxes(ax)

    plt.show()


if __name__ == '__main__':
    main()


# idx = np.array(
#     [3, 11, 12, 13, 14, 15, 23, 24, 25, 26, 27, 35, 36, 37, 38, 39, 59, 60,
#      61,
#      62, 63,
#      71, 72, 73, 74, 75, 83, 84, 85, 86, 87, 95, 96, 97, 98, 99])

# params = np.array([
#     0, 0, 3.739717e+00, 0, 0, 0, 0, 0, 0, 0, 1.404616e-03, 2.603632e-03,
#      3.687379e-02, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.747256e-06, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 1, -7.895465e+00, -9.989408e+00, -9.999000e+00,
#      -9.017598e+00, -9.999000e+00, -10, -10, -10, -10, -10, -10,
#      -8.167491e+00, -9.577779e+00, -9.999414e+00, -9.999001e+00,
#      -9.999794e+00, -9.999041e+00, -10, -10, -10, -10, -10, -10,
#      -8.167491e+00, -9.999005e+00, -9.999753e+00, -9.999603e+00,
#      -9.999015e+00, -9.999049e+00, -10, -10, -10, -10, -10, -10,
#      -8.167491e+00, -7.132977e+00, -9.878768e+00, -9.999935e+00,
#      -9.762116e+00, -9.999904e+00, -10, -10, -10, -10, -10, -10,
#      -8.167491e+00, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
#      -8.167491e+00])

# params_fixed_theta = np.array([
#     0, 0, 4.005568e+00, 0, 0, 0, 0, 0, 0, 0, 9.871867e-04, 3.298028e-03,
#     3.416207e-02, 6.161611e-09, 3.703897e-09, 0, 0, 0, 0, 0, 0, 1,
#     1.882661e-07,
#     1.017504e-09, 1.075642e-09, 1.176205e-09, 1.241630e-09, 0, 0, 0, 0, 0,
#     0, 1,
#     1.049747e-09, 1.901493e-09, 2.295629e-09, 2.202899e-09, 2.146478e-09, 0,
#     0,
#     0,
#     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -7.968431e+00,
#     -9.234633e+00,
#     -9.992394e+00, -9.088640e+00, -9.999577e+00, -10, -10, -10, -10, -10,
#     -10,
#     -8.167491e+00, -9.270446e+00, -9.929881e+00, -9.999384e+00,
#     -9.999122e+00,
#     -9.999002e+00, -10, -10, -10, -10, -10, -10, -8.167491e+00,
#     -9.702740e+00,
#     -9.999582e+00, -9.999172e+00, -9.999000e+00, -9.999001e+00, -10, -10,
#     -10,
#     -10,
#     -10, -10, -8.167491e+00, -7.171749e+00, -9.221302e+00, -9.879141e+00,
#     -9.629096e+00, -9.993798e+00, -10, -10, -10, -10, -10, -10,
#     -8.167491e+00,
#     -10,
#     -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -8.167491e+00])

# params for individual mean div per chrom
# mat1 = np.array([4.15899265273172, 0.00146317887212695, 0.00152870414627810,
#                  0.0335230522719113, 0, 0, 4.29565499232058e-07, 0, 0, 0, 0,
#                  0,
#                  0.000119370444382320, 0, 0, 0, -8.79249953831005,
#                  -9.80856436915452, -9.99942080802024, -8.97063706489873,
#                  -9.99999999588612, -9.01027739675119, -9.99900066920070,
#                  -9.99923038604574, -9.99900026568013, -9.99941093529020,
#                  -9.98494785997175, -9.99999999588612, -9.99999999588612,
#                  -9.99911321551354, -9.99999687242159, -7.03850412517983,
#                  -9.99900205904364, -9.99913429495595, -9.99387797717417,
#                  -9.99941391731486])
# mat2 = np.array([3.66956390000000, 0.00151390880000000, 0.00217520680000000,
#                  0.0377893670000000, 0, 0, 1.81853430000000e-08, 0, 0, 0, 0,
#                  0,
#                  0, 0, 0, 0, -7.77244000000000, -9.99505760000000,
#                  -9.99922430000000, -9.00598420000000, -9.99959980000000,
#                  -9.94161000000000, -9.99914620000000, -9.99909540000000,
#                  -9.99936110000000, -9.99910930000000, -9.99926310000000,
#                  -9.99930410000000, -9.99933650000000, -9.99900610000000,
#                  -9.99900000000000, -7.17393240000000, -9.81973840000000,
#                  -9.99900080000000, -9.70812070000000, -9.99900150000000])