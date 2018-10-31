# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# BSD license.

import numpy as np
from netCDF4 import Dataset

dat = Dataset('data/ncep/gph_winds_ncep_1998_2014.nc')

precip = Dataset('data/ncep/pr_ncep_1998_2014.nc')
print precip.variables

pr1 = precip.variables['pr_wtr'][:, :, :]

trmm = Dataset('data/TRMM3B42V7_daily_ncep_grid_global.nc')
pr2 = trmm.variables['pcp']


# [ 1000.   925.   850.   700.   600.   500.   400.   300.   250.   200. 150.   100.    70.    50.    30.    20.    10.]
lat = dat.variables['lat'][:]
lon = dat.variables['lon'][:]
time = dat.variables['time'][:]
level = dat.variables['level'][:]

data = np.zeros((3, time.shape[0], 3, lat.shape[0], lon.shape[0]))
data[0] = dat.variables['hgt'][:, [1, 5, 8], :, :]
data[1] = dat.variables['uwnd'][:, [1, 5, 8], :, :]
data[2] = dat.variables['vwnd'][:, [1, 5, 8], :, :]

print "levels: ", dat.variables['level'][[1, 5, 8]]

y = 16

index_june = np.arange(151, 181, 1)
index_july = np.arange(181, 212, 1)
index_august = np.arange(212, 243, 1)

for i in range(0, y - 1):
    index_june = np.concatenate((index_june, np.arange(516 + 365 * i, 546 + 365 * i, 1)))
    index_july = np.concatenate((index_july, np.arange(546 + 365 * i, 577 + 365 * i, 1)))
    index_august = np.concatenate((index_august, np.arange(577 + 365 * i, 608 + 365 * i, 1)))



jja_indices = np.sort(np.concatenate((index_june, index_july, index_august)))

monsoons = ['ISM1', 'EUR']
syncs = ['ISM1-EUR']
sync_pairs = [[0, 1]]




clim = np.zeros((5, data.shape[2], data.shape[3], data.shape[4]))
clim[:3] = np.mean(data[:, jja_indices], axis = 1)
clim[3, 0] = np.mean(pr1[jja_indices], axis = 0)
clim[4, 0] = np.mean(pr2[jja_indices], axis = 0)


taumax = 10
shift = [-14, -12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14]

monsoon_comp = np.zeros((len(monsoons), len(shift), 5, data.shape[2], data.shape[3], data.shape[4]))
sync_comp = np.zeros((len(syncs), len(shift), 5, data.shape[2], data.shape[3], data.shape[4]))


monsoon_ano = np.zeros((len(monsoons), len(shift), 5, data.shape[2], data.shape[3], data.shape[4]))
sync_ano = np.zeros((len(syncs), len(shift), 5, data.shape[2], data.shape[3], data.shape[4]))

for direct in ['12', '21']:
    for sth in [90]:
        for cutoff in [8, 10, 12]:
            print "Cut = %d"%cutoff
            for perc in xrange(95, 96):

                sync_times = np.load('data/monsoon_sync_times_tm%d_lpw%d_perc%s_jja_sth%d_nb.npy'%(taumax, cutoff, perc, sth))
                for k in xrange(sync_times.shape[0]):
                    if direct == '12':
                        sync_tidx = sync_times[k, 0]
                    elif direct == '21':
                        sync_tidx = sync_times[k, 1]
                    if m == 1:
                        sync_tidx = np.intersect1d(sync_tidx, np.where(mjo12 > 0)[0])

                    for i in xrange(len(shift)):
                        sync_comp[k, i, :3] = np.mean(data[:, sync_tidx + shift[i]], axis = 1)
                        sync_comp[k, i, 3, 0] = np.mean(pr1[sync_tidx + shift[i]], axis = 0)
                        sync_comp[k, i, 4, 0] = np.mean(pr2[sync_tidx + shift[i]], axis = 0)
                        sync_ano[k, i] = sync_comp[k, i] - clim
                    print "sync times %s: %d"%(syncs[k], sync_tidx.shape[0])
                    np.save('data/sync_hgt_wind_anomalies_jja_perc%d_lp%d_sth%d_mjo%d_%s_nb'%(perc, cutoff, sth, m, direct), sync_ano)
