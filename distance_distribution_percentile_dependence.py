# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.


import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)

cut = 2500.
for tm in [10]:
    for j in [2]:
        c = 0
        linkfracs = np.zeros((3, 20))
        percs = np.zeros((3, 20))
        for perc in np.arange(80, 100):
            print perc
            cu_deg = np.load('cu_degree01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))
            ang_dist = np.load('/p/tmp/boers/ang_dist01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))
            print "ang dist shape:", ang_dist.shape[0]
            ang_dist = np.random.choice(ang_dist, int(1e8))
            print "ang dist shape sampled:", ang_dist.shape[0]
            linkfracs[0, c] = ang_dist.shape[0]
            linkfracs[1, c] = ang_dist[ang_dist <= cut].shape[0]
            linkfracs[2, c] = ang_dist[ang_dist > cut].shape[0]

            percs[0, c] = st.scoreatpercentile(ang_dist, 25)
            percs[1, c] = st.scoreatpercentile(ang_dist, 50)
            percs[2, c] = st.scoreatpercentile(ang_dist, 75)
            c += 1
        np.save('trmm7_linkfracs_global_wd_tm%d_season%d_2_nb_sig005.npy'%(tm, j), linkfracs)
        np.save('trmm7_dist_percs_global_wd_tm%d_season%d_2_nb_sig005.npy'%(tm, j), percs)



lat = np.load('data/trmm7_lat_global.npy')
lon = np.load('data/trmm7_lon_global.npy')
la = len(lat)
lo = len(lon)
n = la * lo


for tm in [10]:
    for j in [2]:
        linkfracs = np.load('data/trmm7_linkfracs_global_wd_tm%d_season%d_2_nb_sig005.npy'%(tm, j))
        percs = np.load('data/trmm7_dist_percs_global_wd_tm%d_season%d_2_nb_sig005.npy'%(tm, j))
        print linkfracs
        fig = plt.figure(figsize = (6,2))
        ax = fig.add_subplot(111)
        ax.plot(np.arange(80, 100), linkfracs[2] / linkfracs[0], marker = 'o', color = 'royalblue')
        ax.set_ylim((.3,1.))
        ax.set_xlabel(r'Rainfall event percentile')
        ax.set_ylabel(r'Fraction of links with $d>$ 2500km')
        ax.yaxis.label.set_color('royalblue')
        ax.tick_params(axis='y', colors='royalblue')
        plt.locator_params(nbins=11)
        ax2 = ax.twinx()
        ax2.plot(np.arange(80, 100), percs[1], marker = 'o', color = 'darkorange')
        ax2.yaxis.label.set_color('darkorange')
        ax2.tick_params(axis='y', colors='darkorange')
        ax2.set_ylabel(r'Median link distance [km]')
        ax2.set_ylim((0.,10000.))
        plt.savefig('pics/global/trmm7/linkfracs/linkfracs_tm%d_season%d_2_nb_sig005.pdf'%(tm, j), bbox_inches = 'tight')
