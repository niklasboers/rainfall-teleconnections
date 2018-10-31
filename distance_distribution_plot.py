# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.



import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
import scipy.optimize as so
import scipy.special as spp
from scipy.constants import pi
from mpmath import gammainc
from sys import float_info
import scipy.stats as st
from scipy.interpolate import interp1d
from sklearn.neighbors import KernelDensity

import matplotlib as mpl
mpl.rc('text', usetex=True)






seasons = ['DJF', 'MAM', 'JJA', 'SON']
x_min = .5 * 1e2
x_max = 2 * 1e3

reg = [""]

data = 'trmm7'

if data == 'trmm7':
    x_min = 1.1*1e2
    x_max = 2.5 * 1e3

elif data == 'gpcp':
    x_min = 3 * 1e2
    x_max = 2.5 * 1e3



for tm in [10]:
   for j in [2]:
      for i in xrange(len(reg)):
         for perc in np.arange(80, 100):
            print "#############################################################################"
            print regions[i]
            print "Season %s"%seasons[j]
            print "%dth Percentile Events"%perc
            lat = np.load('data/%s_lat_global.npy'%data)
            lon = np.load('data/%s_lon_global.npy'%data)
            la = len(lat)
            lo = len(lon)
            nodes = la * lo
            m = nodes * (nodes - 1) / 2
            loghist = np.load('data/ang_dist_loghist40_%s_global%s_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(data, reg[i], perc, tm, j))


            loghist_ref = np.load('data/ang_dist_noe_ref_loghist40_%s_global%s_wd_perc%d_tm%d_season%d_nb_sig005.npy'%(data, reg[i], perc, tm, j))
            kde = np.load('data/ang_dist_log_kde_%s_global%s_wd_perc%d_tm%d_season%d_cor_2_nb_sig005.npy'%(data, reg[i], perc, tm, j))


            logx = loghist_ref[1][:-1] + (loghist_ref[1][1:] - loghist_ref[1][:-1]) / 2.
            tlogx = logx[logx > x_min]

            tkde = kde[logx > x_min]

            tlogy = loghist[0][logx > x_min]
            tlogy_ref = loghist_ref[0][logx > x_min]

            x_cut = loghist[1][loghist[1] < x_max]
            x_cut = x_cut[loghist[1][loghist[1] < x_max] > x_min]
            y_cut = loghist[0][loghist[1][:-1] < x_max]
            y_cut = y_cut[loghist[1][loghist[1] < x_max] > x_min]
            plfit = lambda params: np.sum(np.abs(params[0] * x_cut**(-params[1]) - y_cut))
            pl_params = so.fmin(plfit, [1,1], disp = False)

            fig = plt.figure(figsize = (5,2))
            ax1 = fig.add_subplot(111)

            plt.loglog(logx[logx <= 2500], loghist[0][logx <= 2500], color = 'r', ls = 'None', marker = 'o', fillstyle = 'none', label = r'Distance histogram ($d\leq$ 2500km)')
            plt.loglog(logx[logx > 2500], loghist[0][logx > 2500], color = 'b', ls = 'None', marker = 'o', fillstyle = 'none', label = r'Distance histogram ($d>$ 2500km)')
            plt.loglog(tlogx, pl_params[0] * tlogx**(-pl_params[1]), 'k--', alpha = .7, label = r'Powerlaw fit, $\alpha = %.3f$'%(pl_params[1]))
            plt.loglog(tlogx, tkde, 'k-', alpha = .7, label = 'Great-circle KDE')

            plt.axvline(x = x_max, color = 'm', alpha = .6)

            plt.grid()
            plt.ylim((10**-6,10**-2))
            plt.xlim(10**1, 10**5)
            plt.xlabel('Distance [km]')
            plt.ylabel('PDF')

            plt.legend(loc = 1, ncol = 2, fontsize = 8.)
            plt.savefig('pics/global/%s/densities/%s_ang_dist_simple_pl_perc%d_tm%d_season%d_nb_sig005.pdf'%(data, data, perc, tm, j), bbox_inches = 'tight')
            plt.close()
