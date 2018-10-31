# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# BSD license.



import numpy as np
import scipy.stats as st
from netCDF4 import Dataset

trmm = Dataset('/p/tmp/boers/TRMM3B42_GLOBAL/TRMM3B42V7_daily_global.nc')
lat = trmm.variables['latitude'][:]
lon = trmm.variables['longitude'][:]
time = trmm.variables['time'][:]
pcp = trmm.variables['pcp']
t = time.shape[0]
la = lat.shape[0]
lo = lon.shape[0]
n = la * lo

y = 16
index_seasons = np.zeros((4, 16), dtype = 'object')
index_seasons[0, 0] = np.arange(0,59,1)
index_seasons[1, 0] = np.arange(59,151,1)
index_seasons[2, 0] = np.arange(151,243,1)
index_seasons[3, 0] = np.arange(243,334,1)
for i in range(0, y - 1):
   index_seasons[0, i + 1] = np.arange(334 + 365 * i, 424 + 365 * i, 1)
   index_seasons[1, i + 1] = np.arange(424 + 365 * i, 516 + 365 * i, 1)
   index_seasons[2, i + 1] = np.arange(516 + 365 * i, 608 + 365 * i, 1)
   index_seasons[3, i + 1] = np.arange(608 + 365 * i, 699 + 365 * i, 1)

index_season = np.zeros((4), dtype = 'object')
index_season[0] = np.arange(0,59,1)
index_season[1] = np.arange(59,151,1)
index_season[2] = np.arange(151,243,1)
index_season[3] = np.arange(243,334,1)

for i in range(0, y - 1):
   index_season[0] = np.concatenate((index_season[0], np.arange(334 + 365 * i, 424 + 365 * i, 1)))
   index_season[1] = np.concatenate((index_season[1], np.arange(424 + 365 * i, 516 + 365 * i, 1)))
   index_season[2] = np.concatenate((index_season[2], np.arange(516 + 365 * i, 608 + 365 * i, 1)))
   index_season[3] = np.concatenate((index_season[3], np.arange(608 + 365 * i, 699 + 365 * i, 1)))

def ec_wd(ts, perc):
   th = st.scoreatpercentile(ts[ts > 1], perc)
   ts = np.where(ts > th)[0]
   noe = ts.shape[0]
   if noe > 2 and th > 2:
      return ts, noe, th
   else:
      return np.zeros(1), 0, 0



for perc in xrange(80, 100):
   for j in [2]:
      mnoe = index_season[j].shape[0] * (1 - perc / 100.)
      print mnoe
      ev = np.zeros((n, int(mnoe) + 1), dtype = 'uint16')
      noe = np.zeros(n, dtype = 'uint16')
      th = np.zeros(n)
      for l in xrange(la):
         rain = np.zeros((t, lo))
         precip = pcp[index_season[j], l, :]
         if np.ma.is_masked(precip) is True:
            rain[index_season[j], :] = pcp[index_season[j], l, :].filled(0)
         else:
            rain[index_season[j], :] = pcp[index_season[j], l, :]
         for k in xrange(lo):
            events, noe[l * lo + k], th[l * lo + k] = ec_wd(rain[:,k], perc)
            ev[l * lo + k, :noe[l * lo + k]] = events
      np.save('trmm7_global_wd_score_cor_seasonal_rain_perc%d_season%d'%(perc, j), th)
      np.save('/p/tmp/boers/trmm7_global_wd_events_cor_seasonal_rain_perc%d_season%d'%(perc, j), ev)
      np.save('trmm7_global_wd_noe_cor_seasonal_rain_perc%d_season%d'%(perc, j), noe)

      ev_nb = np.zeros_like(ev)
      nob = np.zeros_like(noe)

      for i in xrange(ev.shape[0]):
          noce = np.where(np.diff(ev[i]) == 1)[0].shape[0]
          ev_nb[i] = np.concatenate((np.delete(ev[i], np.where(np.diff(ev[i]) == 1)[0]+1), np.zeros(noce)))
          nob[i] = np.where(ev_nb[i] != 0)[0].shape[0]
      np.save('/p/tmp/boers/trmm7_global_wd_bursts_cor_seasonal_rain_perc%d_season%d.npy'%(perc, j), ev_nb)
      np.save('trmm7_global_wd_nob_cor_seasonal_rain_perc%d_season%d'%(perc, j), nob)
