# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.



import time
import numpy as np
import powerlaw
import weave
from scipy.stats import chi2
from scipy.stats import gaussian_kde
import scipy.optimize as so
import scipy.special as spp
from scipy.constants import pi
from mpmath import gammainc
from sys import float_info
import scipy.stats as st
from scipy.interpolate import interp1d

from mpl_toolkits.basemap import shiftgrid
from sklearn.neighbors import KernelDensity
# from fitting_functions import *
import matplotlib.pyplot as plt
from basemap_fct import basemap

import matplotlib as mpl
mpl.rc('text', usetex=True)

def coordinate_indices_from_ra(trmm_lat, trmm_lon, lat_max, lat_min, lon_max, lon_min):
	la = len(trmm_lat)
	lo = len(trmm_lon)
	n = la * lo
	indices = np.arange(n).reshape((la, lo))
	lat_max = trmm_lat[np.argmin(np.abs(trmm_lat - lat_max))]
	lat_min = trmm_lat[np.argmin(np.abs(trmm_lat - lat_min))]
	lon_max = trmm_lon[np.argmin(np.abs(trmm_lon - lon_max))]
	lon_min = trmm_lon[np.argmin(np.abs(trmm_lon - lon_min))]
	ra_indices = indices[np.where(trmm_lat == lat_max)[0] : np.where(trmm_lat == lat_min)[0] + 1 , np.where(trmm_lon == lon_min)[0] : np.where(trmm_lon == lon_max)[0] + 1 ]
	return np.unique(ra_indices.flatten())


def get_index_from_coord(lats, lons, lat_t, lon_t):
   from itertools import product
   lat = lats[np.argmin(np.abs(lats - lat_t))]
   lon = lons[np.argmin(np.abs(lons - lon_t))]
   coords = np.array(list(product(lats, lons)))
   idx = np.intersect1d(np.where(coords[:, 0] == lat), np.where(coords[:,1] == lon))
   return idx


def geographic_link_dist_ref_sn(lat, lon, lat_t, lon_t, noe, frac):
   idx = int(get_index_from_coord(lat, lon, lat_t, lon_t))
   la = len(lat)
   lo = len(lon)
   n = la * lo
   lat_seq = np.repeat(lat, lo)
   lon_seq = np.tile(lon, la)
   sin_lon = np.sin(lon_seq * np.pi / 180)
   cos_lon = np.cos(lon_seq * np.pi / 180)
   sin_lat = np.sin(lat_seq * np.pi / 180)
   cos_lat = np.cos(lat_seq * np.pi / 180)
   rn = np.where(noe > 2)[0].shape[0]
   ang_dist = np.zeros(int(2 * rn * frac), dtype = 'uint16')
   var = ['cos_lat', 'sin_lat', 'cos_lon', 'sin_lon', 'ang_dist', 'n', 'noe', 'idx', 'frac']
   code = r"""
   srand48((unsigned int) time(NULL));
   long long j, u;
   long long c = 0;
   float expr;
   float r;
 	for (j = 0; j < n; j++) {
      if(noe[idx] > 2 && noe[j] > 2){
         r = drand48();
         if(r < frac){
           	expr = sin_lat[idx]*sin_lat[j] + cos_lat[idx]*cos_lat[j] * (sin_lon[idx]*sin_lon[j] + cos_lon[idx]*cos_lon[j]);
   			if (expr > 1) {
   				expr = 1.;
   			}
   			else if (expr < -1) {
               	expr = -1.;
   			}
   			ang_dist[c] = (int) (acos(expr) * 6371);
            c++;
         }
      }
   }
   """
   weave.inline(code, var)
   return ang_dist[ang_dist != 0]


seasons = ['DJF', 'MAM', 'JJA', 'SON']
x_min = 100.
func_names = ['PL', 'TPL', 'SEXP', 'EXP', 'LN']
regions = ["Global", "Land", "Ocean"]
reg = [""]

data = 'trmm7'

lat = np.load('data/%s_lat_global.npy'%data)
lon = np.load('data/%s_lon_global.npy'%data)

la = len(lat)
lo = len(lon)
nodes = la * lo
m = nodes * (nodes - 1) / 2

perc = 96
i = 0

bws = ['02']
regions = ['nism']

for bw in bws:
    for perc in [80, 85, 90, 94, 95, 96]:
        print perc
        for tm in [3,10,30]:
          print '######################'
          print "tm ", tm
          for region in regions:

             noe = np.load('data/%s_global%s_wd_nob_cor_seasonal_rain_perc%d_season%d.npy'%(data, reg[i], perc, j))

             links = np.load('data/%s_geo_links_samples_perc%d_tm%d_season%d_%s_nb_sig005.npy'%(data, perc, tm, j, region))
             gdists = np.load('data/%s_geo_gdists_samples_perc%d_tm%d_season%d_%s_nb_sig005.npy'%(data, perc, tm, j, region))

             index = get_index_from_coord(lat, lon, lat_t, lon_t)
             ang_dist = gdists.copy()

             ang_dist_ref = geographic_link_dist_ref_sn(lat, lon, lat_t, lon_t, noe, .1)
             ang_dist_ref = ang_dist_ref[ang_dist_ref > x_min]
             logbins_reg = np.logspace(np.log10(ang_dist_ref.min()), np.log10(ang_dist_ref.max()), 30)
             loghist_ref_reg = np.histogram(ang_dist_ref, bins = logbins_reg, density = True)
             logx_reg = logbins_reg[:-1] + (logbins_reg[1:] - logbins_reg[:-1]) / 2.
             loghist_reg = np.histogram(ang_dist, bins = logbins_reg, density = True)

             start = time.time()
             kernel = gaussian_kde(ang_dist_ref)
             end =  time.time()

             start = time.time()
             kde_ref = kernel.evaluate(logx_reg)

             kde_reg = kernel.evaluate(logx_reg)

             tkde_reg = kde_reg[logx_reg > x_min]

             kernel = kernel.evaluate(ang_dist)
             end =  time.time()

             tlogx_reg = logx_reg
             tkde_reg = kde_reg[logx_reg > x_min]


             tlogy_reg = loghist_reg[0]
             tlogy_ref_reg = loghist_ref_reg[0][logx_reg > x_min]

             from itertools import product
             coords = np.array(list(product(lat, lon)))
             lat0 = coords[index, 0]
             lon0 = coords[index, 1]
             lat1 = coords[links, 0]
             lon1 = coords[links, 1]
             nol = links.shape[0]


             links_short = links[gdists < 2500.]
             links_long = links[gdists > 2500.]

             ang_dist_short = ang_dist[ang_dist < 2500.]
             ang_dist_long = ang_dist[ang_dist >= 2500.]


             lat1_short = coords[links_short, 0]
             lon1_short = coords[links_short, 1]

             lat1_long = coords[links_long, 0]
             lon1_long = coords[links_long, 1]

             dat = np.load('data/%s_geo_link_tc_sphkde_perc%d_tm%d_season%d_%s_bw%s_nb_sig005.npy'%(data, perc, tm, j, region, bw))
             print dat.shape
             mean, std, perc90, perc95, perc99, perc995, perc999 = np.load('data/%s_geo_link_tc_sphkde_stats_perc%d_tm%d_season%d_%s_bw%s_nb.npy'%(data, perc, 10, j, region, bw))


             sig = perc999

             fac = np.intersect1d(np.where(dat.flatten() > sig.flatten())[0], links)
             ix = np.in1d(links.ravel(), fac).reshape(links.shape)
             lat1 = coords[links[ix], 0]
             lon1 = coords[links[ix], 1]
             lat1_uc = coords[links, 0]
             lon1_uc = coords[links, 1]
             area = np.ones((len(lat), len(lon))) * -1
             area[dat>sig] = 10
             area[0] = -1
             area[-1] = -1
             area_noe = np.ones((len(lat), len(lon)))
             area_noe[noe.reshape((la, lo)) >= 3] = -1

             fac_short = np.intersect1d(np.where(dat.flatten() > sig.flatten())[0], links_short)
             ix_short = np.in1d(links_short.ravel(), fac_short).reshape(links_short.shape)
             lat1_short = coords[links_short[ix_short], 0]
             lon1_short = coords[links_short[ix_short], 1]

             fac_long = np.intersect1d(np.where(dat.flatten() > sig.flatten())[0], links_long)
             ix_long = np.in1d(links_long.ravel(), fac).reshape(links_long.shape)
             lat1_long = coords[links_long[ix_long], 0]
             lon1_long = coords[links_long[ix_long], 1]

             print "area fraction of significant links = ", dat.flatten()[dat.flatten() > sig.flatten()].shape[0] / float(dat.flatten().shape[0])
             probs = .01 * np.ones_like(links)
             probs = probs[ix]
             probs[0] = .06
             print "fraction of significant links = ", links[ix].shape[0] / float(links.shape[0])

             ### uncorrected:
             datnan = np.ones((len(lat), len(lon))) * -1
             fig = plt.figure(figsize = (5,2))
             fig.add_subplot(111)
             basemap(datnan, lat, lon, res = 'c', shift = lon_t, proj = 'moll', loc = (lon_t, lat_t), lsmask = True, hillshade = False, etopo = False, scale = .5, gc_lon0 = lon0, gc_lat0 = lat0, gc_lon1 = lon1_uc, gc_lat1 = lat1_uc, gdists = probs, gc_contours = np.linspace(0.05, probs.max(), 6), gc_color = None, gc_cbar = False, gc_lw = .2, color = plt.cm.Blues, alpha = .5, colorbar = False, cbar_title = False, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-50,75,25))
             plt.savefig('pics/reg_kde_estimates/maps/%s/%s_geo_link_tc_cor_perc%d_tm%d_season%d_moll002_%s_vs_links_uncor_nokde_perc99_sig005.png'%(data, data, perc, tm, j, region), dpi = 300, bbox_inches = 'tight')

             ### corrected:
             sigdat = dat.copy()
             sigdat[dat > mean + 5 * std] = 5.5
             sigdat[dat <= mean + 5 * std] = 4.5
             sigdat[dat <= mean + 4 * std] = 3.5
             sigdat[dat <= mean + 3 * std] = 2.5
             sigdat[dat <= mean + 2 * std] = 1.5
             sigdat[-1] = 1.5
             sigdat[0] = 1.5
             conts = [2., 3., 4., 5.]

             fig = plt.figure(figsize = (5,2))
             fig.add_subplot(111)
             basemap(sigdat, lat, lon, res = 'c', shift = lon_t, proj = 'moll', loc = (lon_t, lat_t), lsmask = True, hillshade = False, etopo =False, scale = .5, gc_lon0 = lon0, gc_lat0 = lat0, gc_lon1 = lon1_long, gc_lat1 = lat1_long, gc_lon2 = lon1_short, gc_lat2 = lat1_short, gdists = probs, gc_contours = np.linspace(0.05, probs.max(), 6), gc_color = None, gc_cbar = False, gc_lw = .2, contours = conts, alpha = .5, area1 = area, line_color1 = 'k', hatch = area_noe, color = plt.cm.Blues, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-50,75,25), cbar_title = r'[Standard deviations $\sigma$ above mean]')
             plt.savefig('pics/reg_kde_estimates/maps/%s/%s_geo_link_tc_cor_perc%d_tm%d_season%d_moll002_%s_vs_dragon_links_cor_bw%s_perc99_nb_sig005.png'%(data, data, perc, tm, j, region, bw), dpi = 300, bbox_inches = 'tight')
             mdats = mean
             fig = plt.figure()
             fig.add_subplot(111, aspect = .4)
             basemap(mdats, lat, lon, res = 'c', shift = lon_t, proj = 'moll', loc = (lon_t, lat_t), lsmask = True, hillshade = False, etopo = False,scale = .5, contours = np.linspace(mdats.min(), mdats.max(), 11), area1 = area, line_color1 = 'w', hatch = area_noe, alpha = .5, color = plt.cm.plasma_r, colorbar = True, extend = 'max', meridians = np.arange(0,360, 40), parallels = np.arange(-50,75,25))
             plt.savefig('pics/reg_kde_estimates/maps/%s/%s_geo_link_tc_cor_perc%d_tm%d_season%d_moll002_%s_vs_meansig_bw%s_perc99_nb_sig005.pdf'%(data, data, perc, tm, j, region, bw), bbox_inches = 'tight')

             mdats = std
             fig = plt.figure()
             fig.add_subplot(111, aspect = .4)
             basemap(mdats, lat, lon, res = 'c', shift = lon_t, proj = 'moll', loc = (lon_t, lat_t), lsmask = True, hillshade = False, etopo = False,scale = .5, contours = np.linspace(mdats.min(), mdats.max(), 11), area1 = area, line_color1 = 'w', hatch = area_noe, alpha = .5, color = plt.cm.plasma_r, colorbar = True, extend = 'max', meridians = np.arange(0,360, 40), parallels = np.arange(-50,75,25))
             plt.savefig('pics/reg_kde_estimates/maps/%s/%s_geo_link_tc_cor_perc%d_tm%d_season%d_moll002_%s_vs_cor2s_stdsig_bw%s_perc99_nb_sig005.pdf'%(data, data, perc, tm, j, region, bw), bbox_inches = 'tight')


             loghist_reg_cor = np.histogram(ang_dist_cor, bins = logbins_reg, density = True)

             x_cut = loghist_reg_cor[1][loghist_reg_cor[1] < x_max]
             y_cut = loghist_reg_cor[0][loghist_reg_cor[1][:-1] < x_max]

             plfit = lambda params: np.sum(np.abs(params[0] * x_cut**(-params[1]) - y_cut))
             pl_params = so.fmin(plfit, [1,1], disp = False)

             fig = plt.figure()
             ax1 = fig.add_subplot(111, aspect = .35)

             kde_test = gaussian_kde(ang_dist, bw_method = .1 * ang_dist.shape[0]**(-1. / 5))


             plt.loglog(tlogx_reg, tkde_reg, 'k-', alpha = .7, label = 'Great-circle KDE')
             plt.loglog(logx_reg, loghist_reg[0], color = 'r', ls = 'None', marker = 'o', fillstyle = 'none', label = 'Histogram (all)')
             plt.loglog(logx_reg, loghist_reg_cor[0], color = 'b', ls = 'None', marker = 'o', fillstyle = 'none', label = 'Histogram (corrected)')
             plt.loglog(tlogx_reg, pl_params[0] * tlogx_reg**(-pl_params[1]), 'k--', alpha = .7, label = r'Powerlaw fit, $\alpha = %.3f$'%(pl_params[1]))

             plt.grid()

             plt.ylim((10**-6,10**-2))
             plt.xlim(10**1, 10**5)
             plt.xlabel('Distance [km]')
             plt.ylabel('PDF')

             plt.legend(loc = 1, ncol = 2, fontsize = 8.)

             plt.savefig('pics/reg_kde_estimates/densities/%s/%s_ang_dist_pl_%s_mle_hist40_noe_ref_global_wd_perc%d_tm%d_season%d_1fig_perc99_nb_sig005.pdf'%(data, data, region, perc, tm, j), bbox_inches = 'tight')
             plt.close()
