# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# BSD license.



import numpy as np

lat = np.load('data/trmm7_lat_global.npy')
lon = np.load('data/trmm7_lon_global.npy')
print lat
print lon

la = len(lat)
lo = len(lon)
nodes = la * lo
m = nodes * (nodes - 1) / 2

idx_all = np.arange(nodes).reshape((la, lo))

idx_tropical = idx_all[np.where(np.logical_and(lat > -23.5, lat < 23.5))[0]].flatten()
idx_north = idx_all[np.where(lat > 23.5)[0]].flatten()
idx_south = idx_all[np.where(lat < -23.5)[0]].flatten()
idx_extra = np.concatenate((idx_north, idx_south))

print idx_tropical

idx_tropical = np.random.choice(idx_tropical, int(.1 * idx_tropical.shape[0]))
idx_north = np.random.choice(idx_north, int(.1 * idx_north.shape[0]))
idx_south = np.random.choice(idx_south, int(.1 * idx_south.shape[0]))
idx_extra = np.concatenate((idx_north, idx_south))

# dat = np.zeros(nodes)
# dat[idx_tropical] = 3
# dat[idx_north] = 1
# dat[idx_south] = 2
# dat[idx_extra] = 1
#
# import matplotlib.pyplot as plt
# from basemap_fct import basemap
# plt.figure()
# basemap(dat.reshape((la, lo)), lat, lon)
# plt.show()

x_min = 1e2
x_max = 2.5 * 1e3


for tm in [10]:
    print tm
    for j in [0,2]:
        print j
        for perc in np.arange(95, 96):
            print perc
            deg = np.load('degree01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))
            cu_deg = np.load('cu_degree01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))
            nlist = np.load('/p/tmp/boers/nlist01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))
            ang_dist = np.load('/p/tmp/boers/ang_dist01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))

            geodists_tropical = np.array([])
            geodists_tropical_ns = np.array([])
            geodists_only_tropical_ns = np.array([])
            geodists_north = np.array([])
            geodists_south = np.array([])
            for i in idx_tropical:
                geodists_tropical =  np.concatenate((geodists_tropical, ang_dist[cu_deg[i] : cu_deg[i + 1]][np.where(np.in1d(nlist[cu_deg[i] : cu_deg[i + 1]], idx_tropical))[0]]))
                geodists_tropical_ns =  np.concatenate((geodists_tropical_ns, ang_dist[cu_deg[i] : cu_deg[i + 1]]))
                geodists_only_tropical_ns =  np.concatenate((geodists_only_tropical_ns, ang_dist[cu_deg[i] : cu_deg[i + 1]][np.where(np.in1d(nlist[cu_deg[i] : cu_deg[i + 1]], idx_extra))[0]]))
            print "idx_tropical done"
            for i in idx_north:
                geodists_north =  np.concatenate((geodists_north, ang_dist[cu_deg[i] : cu_deg[i + 1]][np.where(np.in1d(nlist[cu_deg[i] : cu_deg[i + 1]], idx_north))[0]]))
            print "idx_north done"
            for i in idx_south:
                geodists_south =  np.concatenate((geodists_south, ang_dist[cu_deg[i] : cu_deg[i + 1]][np.where(np.in1d(nlist[cu_deg[i] : cu_deg[i + 1]], idx_south))[0]]))
            print "idx_south done"
            np.save('trmm7_tropical_ns_link_distsperc%d_tm%d_season%d_nb_sig005_test2.npy'%(perc, tm, j), geodists_tropical_ns)
            np.save('trmm7_only_tropical_ns_link_distsperc%d_tm%d_season%d_nb_sig005_test2.npy'%(perc, tm, j), geodists_only_tropical_ns)
            np.save('trmm7_tropical_link_distsperc%d_tm%d_season%d_nb_sig005_test2.npy'%(perc, tm, j), geodists_tropical)
            np.save('trmm7_north_link_dists_perc%d_tm%d_season%d_nb_sig005_test2.npy'%(perc, tm, j), geodists_north)
            np.save('trmm7_south_link_dists_perc%d_tm%d_season%d_nb_sig005_test2.npy'%(perc, tm, j), geodists_south)

            geodists_only_tropical_ns = np.load('data/trmm7_only_tropical_ns_link_distsperc%d_tm%d_season%d_nb_sig005.npy'%(perc, tm, j))
            geodists_tropical_ns = np.load('data/trmm7_tropical_ns_link_distsperc%d_tm%d_season%d_nb_sig005.npy'%(perc, tm, j))
            geodists_tropical = np.load('data/trmm7_tropical_link_distsperc%d_tm%d_season%d_nb_sig005.npy'%(perc, tm, j))
            geodists_north = np.load('data/trmm7_north_link_dists_perc%d_tm%d_season%d_nb_sig005.npy'%(perc, tm, j))
            geodists_south = np.load('data/trmm7_south_link_dists_perc%d_tm%d_season%d_nb_sig005.npy'%(perc, tm, j))

            geodists_tropical_ns = np.concatenate((geodists_tropical_ns, geodists_only_tropical_ns))

            loghist = np.load('data/ang_dist_loghist40_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(perc, tm, j))
            loghist_ref = np.load('data/ang_dist_noe_ref_loghist40_trmm7_global_wd_perc%d_tm%d_season%d.npy'%(perc, tm, j))
            kde = np.load('data/ang_dist_log_kde_trmm7_global_wd_perc%d_tm%d_season%d_cor_2.npy'%(perc, tm, j))


            logx = loghist_ref[1][:-1] + (loghist_ref[1][1:] - loghist_ref[1][:-1]) / 2.
            tlogx = logx[logx > x_min]

            tkde = kde[logx > x_min]
            loghist_only_tropical_ns = np.histogram(geodists_only_tropical_ns, bins = loghist_ref[1], density = True)
            loghist_tropical_ns = np.histogram(geodists_tropical_ns, bins = loghist_ref[1], density = True)
            loghist_tropical = np.histogram(geodists_tropical, bins = loghist_ref[1], density = True)
            loghist_north = np.histogram(geodists_north, bins = loghist_ref[1], density = True)
            loghist_south = np.histogram(geodists_south, bins = loghist_ref[1], density = True)

            import matplotlib.pyplot as plt
            import matplotlib as mpl
            mpl.rc('text', usetex=True)

            fig = plt.figure()
            ax1 = fig.add_subplot(111, aspect = 1.)
            plt.loglog(logx, loghist[0], color = 'k', ls = 'None', marker = 'o', fillstyle = 'none', label = r'All', alpha = .8)
            plt.loglog(logx, loghist_tropical[0], color = 'g', ls = 'None', marker = 's', fillstyle = 'none', label = r'Tropics', alpha = .8)
            plt.loglog(logx, loghist_tropical_ns[0], color = 'm', ls = 'None', marker = 'D', fillstyle = 'none', label = r'Tropics - Tropics - Extratropics', alpha = .8)
            plt.loglog(logx, loghist_only_tropical_ns[0], color = 'c', ls = 'None', marker = '*', fillstyle = 'none', label = r'Tropics - Extratropics', alpha = .8)
            plt.loglog(logx, loghist_north[0], color = 'r', ls = 'None', marker = '^', fillstyle = 'none', label = r'NH', alpha = .7)
            plt.loglog(logx, loghist_south[0], color = 'b', ls = 'None', marker = 'v', fillstyle = 'none', label = r'SH', alpha = .7)

            plt.grid()

            plt.ylim((.5 * 10**-5,10**-2))
            plt.xlim(10**1, 10**5)
            plt.xlabel('Distance [km]')
            plt.ylabel('PDF')

            plt.legend(loc = 1, ncol = 2, fontsize = 8.)

            plt.savefig('pics/global/trmm7/densities/trmm7_ang_dist_tns_pl_perc%d_tm%d_season%d_nb_sig005.pdf'%(perc, tm, j), bbox_inches = 'tight')
            plt.close()
