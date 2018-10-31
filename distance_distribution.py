# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.




import time
import numpy as np
from scipy import weave
import scipy.optimize as so
from scipy.stats import gaussian_kde
import scipy.special as spp
from scipy.stats import chi2
from scipy.constants import pi
from sys import float_info



def geographic_link_dist_ref(lat, lon, noe):
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
   print rn
   ang_dist = np.zeros(.5 * rn * (rn - 1), dtype = 'uint16')
   var = ['cos_lat', 'sin_lat', 'cos_lon', 'sin_lon', 'ang_dist', 'n', 'noe']
   code = r"""
   long long i, j, u;
   long long c = 0;
   float expr;
   for (i = 0; i < n; i++) {
    	for (j = 0; j < i; j++) {
         if(noe[i] > 2 && noe[j] > 2){
           	expr = sin_lat[i]*sin_lat[j] + cos_lat[i]*cos_lat[j] * (sin_lon[i]*sin_lon[j] + cos_lon[i]*cos_lon[j]);
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
   return ang_dist

def geographic_link_dist_ref_sample(lat, lon, noe, frac):
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
   print rn
   ang_dist = np.zeros(int(2 * frac * .5 * rn * (rn - 1)), dtype = 'uint16')
   var = ['cos_lat', 'sin_lat', 'cos_lon', 'sin_lon', 'ang_dist', 'n', 'noe', 'frac']
   code = r"""
   srand48((unsigned int) time(NULL));
   long long i, j, u;
   long long c = 0;
   float expr;
   float r;
   for (i = 0; i < n; i++) {
    	for (j = 0; j < i; j++) {
         if(noe[i] > 2 && noe[j] > 2){
            r = drand48();
            if(r < frac){
              	expr = sin_lat[i]*sin_lat[j] + cos_lat[i]*cos_lat[j] * (sin_lon[i]*sin_lon[j] + cos_lon[i]*cos_lon[j]);
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
   }
   """
   weave.inline(code, var)
   return ang_dist[ang_dist != 0]


x_min = 100.

frac = 0.001

reg = [""]

for tm in [10]:
   for j in [0]:
      for i in xrange(len(reg)):
         for perc in xrange(95, 96):
            print '#############'
            print 'Season %d'%j
            print '%s'%reg[i]
            print '%dth Percentile'%perc
            lat = np.load('trmm7_lat_global.npy')
            lon = np.load('trmm7_lon_global.npy')
            la = len(lat)
            lo = len(lon)
            nodes = la * lo
            m = nodes * (nodes - 1) / 2
            noe = np.load('trmm7_global%s_wd_nob_cor_seasonal_rain_perc%d_season%d.npy'%(reg[i], perc, j))
            start = time.time()
            ang_dist_ref = geographic_link_dist_ref_sample(lat, lon, noe, frac)
            end = time.time()
            dt = end - start
            print "ang_dist_ref took ", dt
            logbins = np.logspace(np.log10(ang_dist_ref.min()), np.log10(ang_dist_ref.max()), 40)
            loghist_ref = np.histogram(ang_dist_ref, bins = logbins, density = True)
            np.save('ang_dist_noe_ref_loghist40_trmm7_global%s_wd_perc%d_tm%d_season%d_nb_sig005.npy'%(reg[i], perc, tm, j), loghist_ref)
            ang_dist_ref = ang_dist_ref[ang_dist_ref > x_min]
            print "original ang_dist_ref size:", ang_dist_ref.shape[0]

            samples1 = int(1e6)
            ang_dist_ref = np.random.choice(ang_dist_ref, samples1)

            kernel = gaussian_kde(ang_dist_ref)
            logx = logbins[:-1] + (logbins[1:] - logbins[:-1]) / 2.
            kde = kernel.evaluate(logx)
            np.save('ang_dist_log_kde_trmm7_global%s_wd_perc%d_tm%d_season%d_cor_2_nb_sig005.npy'%(reg[i], perc, tm, j), kde)
            print "kernel estimated"
            del ang_dist_ref

            ang_dist = np.load('/p/tmp/boers/ang_dist01_trmm7_global%s_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(reg[i], perc, tm, j))
            logbins = np.logspace(np.log10(ang_dist.min()), np.log10(ang_dist.max()), 40)
            ang_dist = np.random.choice(ang_dist, int(1e7))
            loghist = np.histogram(ang_dist, bins = logbins, density = True)
            np.save('ang_dist_loghist40_trmm7_global%s_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(reg[i], perc, tm, j), loghist)
