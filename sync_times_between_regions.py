# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# BSD license.


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.signal import filtfilt, cheby1, butter, lombscargle, argrelmax, welch, periodogram
from lagcorr import lagcorr
import cython_func


def cheby_lowpass(cutoff, fs, order, rp):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = cheby1(order, rp, normal_cutoff, btype='low', analog=False)
    return b, a

def cheby_lowpass_filter(x, cutoff, fs, order, rp):
    b, a = cheby_lowpass(cutoff, fs, order, rp)
    y = filtfilt(b, a, x)
    return y


def coordinate_indices_from_ra(lat, lon, lat_max, lat_min, lon_max, lon_min):
    la = len(lat)
    lo = len(lon)
    n = la * lo
    indices = np.arange(n).reshape((la, lo))
    lat_max = lat[np.argmin(np.abs(lat - lat_max))]
    lat_min = lat[np.argmin(np.abs(lat - lat_min))]
    lon_max = lon[np.argmin(np.abs(lon - lon_max))]
    lon_min = lon[np.argmin(np.abs(lon - lon_min))]
    ra_indices = indices[np.where(lat == lat_min)[0][0] : np.where(lat == lat_max)[0][0] + 1 , np.where(lon == lon_min)[0][0] : np.where(lon == lon_max)[0][0] + 1 ]
    return np.unique(ra_indices.flatten())


y = 16

index_june = np.arange(151, 181, 1)
index_july = np.arange(181, 212, 1)
index_august = np.arange(212, 243, 1)
index_september = np.arange(243, 273, 1)

for i in range(0, y - 1):
    index_june = np.concatenate((index_june, np.arange(516 + 365 * i, 546 + 365 * i, 1)))
    index_july = np.concatenate((index_july, np.arange(546 + 365 * i, 577 + 365 * i, 1)))
    index_august = np.concatenate((index_august, np.arange(577 + 365 * i, 608 + 365 * i, 1)))


jja_indices = np.sort(np.concatenate((index_june, index_july, index_august)))

tlen = 5844
data = 'trmm7'
lat = np.load('data/%s_lat_global.npy'%data)
lon = np.load('data/%s_lon_global.npy'%data)
n = lat.shape[0] * lon.shape[0]

monsoon_indices = np.zeros(2, dtype = 'object')
## ISM1
lat_max = 32.
lat_min = 25.
lon_max = 88.
lon_min = 71.
monsoon_indices[0] = coordinate_indices_from_ra(lat, lon, lat_max, lat_min, lon_max, lon_min)



## EUR
lat_max = 50.
lat_min = 42.
lon_max = 15.
lon_min = 3.
monsoon_indices[1] = coordinate_indices_from_ra(lat, lon, lat_max, lat_min, lon_max, lon_min)



monsoons = ['NISM','EUR']
syncs = ['ISM1-EUR']
sync_pairs = [[0, 1]]



taumax = 10

for perc in xrange(95, 96):
    ev = np.zeros((n, 1))
    for j in xrange(4):
        monsoon_ts = np.zeros((len(monsoons), tlen))
        evs = np.load('data/trmm7_global_wd_bursts_cor_seasonal_rain_perc%d_season%d.npy'%(perc, j))
        ev = np.concatenate((ev, evs), axis = 1)
    ev = np.sort(ev, axis = 1)
    print ev
    for i in xrange(7):
        for t in xrange(1, tlen):
            monsoon_ts[i, t] = np.where(ev[monsoon_indices[i]] == t)[0].shape[0]
    np.save('data/monsoon_timeseries_perc%d_nb'%(perc), monsoon_ts)
    i = 0
    t = np.zeros((len(syncs), tlen))
    t12 = np.zeros((len(syncs), tlen))
    t21 = np.zeros((len(syncs), tlen))

    for p in sync_pairs:
        e1 = np.array(ev[monsoon_indices[p[0]]], dtype = 'float')
        e2 = np.array(ev[monsoon_indices[p[1]]], dtype = 'float')
        t[i], t12[i], t21[i], delay12, delay21 = cython_func.ESreg(e1, e2, t[i], t12[i], t21[i], taumax, tlen)

        i += 1
    sync_ts = [t, t12, t21]
    np.save('data/monsoon_sync_timeseries_tm%d_perc%s_nb.npy'%(taumax, perc), sync_ts)




for sth in [90]:
    for cutoff in [8, 10, 12]:
        print "Cutoff at %d days "%cutoff
        for perc in xrange(95, 96):
            monsoon_ts = np.load('data/monsoon_timeseries_perc%d_nb.npy'%(perc))
            t, t12, t21 = np.load('data/monsoon_sync_timeseries_tm%d_perc%s_nb.npy'%(taumax, perc))


            monsoon_times = np.zeros(len(monsoons), dtype = 'object')
            for i in xrange(len(monsoons)):
                rmt = cheby_lowpass_filter(monsoon_ts[i], .95 * 1. / cutoff, 1, 8, .05)
                locmax = np.array(argrelmax(rmt)[0])
                monsoon_times[i] = np.intersect1d(locmax, np.where(rmt > st.scoreatpercentile(rmt[jja_indices], sth))[0])

                monsoon_times[i] = np.intersect1d(monsoon_times[i], jja_indices)
                print "monsoon times %s: "%monsoons[i], monsoon_times[i].shape[0]


            np.save('data/monsoon_times_tm%d_lp%d_perc%s_jja_sth%d_nb.npy'%(taumax, cutoff, perc, sth), monsoon_times)


            i = 0

            sync_times = np.zeros((len(sync_pairs), 2), dtype = 'object')
            for p in sync_pairs:

                rmt = cheby_lowpass_filter(t[i] + t12[i], .95 * 1. / cutoff, 1, 8, .05)
                locmax = np.array(argrelmax(rmt)[0])
                sync_times[i, 0] = np.intersect1d(locmax, np.where(rmt > st.scoreatpercentile(rmt[jja_indices], sth))[0])


                rmt = cheby_lowpass_filter(t[i] + t21[i], .95 * 1. / cutoff, 1, 8, .05)
                locmax = np.array(argrelmax(rmt)[0])
                sync_times[i, 1] = np.intersect1d(locmax, np.where(rmt > st.scoreatpercentile(rmt[jja_indices], sth))[0])


                sync_times[i, 0] = np.intersect1d(sync_times[i, 0], jja_indices)
                sync_times[i, 1] = np.intersect1d(sync_times[i, 1], jja_indices)
                print "sync times %s 12: "%syncs[i], sync_times[i, 0].shape[0]
                print "sync times %s 21: "%syncs[i], sync_times[i, 1].shape[0]


                i += 1

            np.save('data/monsoon_sync_times_tm%d_lpw%d_perc%s_jja_sth%d_nb.npy'%(taumax, cutoff, perc, sth), sync_times)

            for s in xrange(len(sync_pairs)):
                m1 = sync_pairs[s][0]
                m2 = sync_pairs[s][1]
                t1 = cheby_lowpass_filter(monsoon_ts[m1], .95 * 1. / cutoff, 1, 8, .05)
                t2 = cheby_lowpass_filter(monsoon_ts[m2], .95 * 1. / cutoff, 1, 8, .05)
                ts12 = cheby_lowpass_filter(t[s] + t12[s], .95 * 1. / cutoff, 1, 8, .05)
                ts21 = cheby_lowpass_filter(t[s] + t21[s], .95 * 1. / cutoff, 1, 8, .05)


                lag = np.arange(-40, 41)
                lcor = lagcorr(t1, t2, lag, cor = 'S')[:, 0]
                pv = lagcorr(t1, t2, lag, cor = 'S')[:, 1]
                maxcor_idx = argrelmax(lcor)[0]
                maxcor_idx = maxcor_idx[pv[maxcor_idx] < .05]
                maxcor_idx = maxcor_idx[lcor[maxcor_idx] > 0.]
                mci = maxcor_idx[np.argmin(np.abs(maxcor_idx - 40))]
                mlag = lag[mci]
                mcor = lcor[mci]

                fig = plt.figure(figsize = (10, 2))
                ax = fig.add_subplot(111)
                ax.plot(lag, lcor, 'k-', label = 'Cor(%s,%s)'%(monsoons[m1], monsoons[m2]))
                ax.axvline(x = mlag, color = 'k', lw = 2)
                ax.axvline(x = 0, color = 'k', ls = '--', lw = .5)
                ax.axhline(y = 0, color = 'k', ls = '--', lw = .5)
                ax.set_ylabel(r'Correlation')
                ax2 = ax.twinx()
                ax2.plot(lag, pv, label = r'$p$-value')
                ax2.axhline(y = 0.05, ls = '--', lw = .5)
                ax2.set_ylabel(r'$p$-value')
                ax.set_xlabel('Lag [days]')
                ax.set_ylabel('Correlation')
                ax.legend(loc = 2)
                ax2.legend(loc = 1)
                plt.savefig('pics/monsoon_sync_times/lags/monsoon_sync_lags_%s_tm%d_lpw%d_perc%s_jja_sth%d_nb.pdf'%(syncs[s], taumax, cutoff, perc, sth), bbox_inches = 'tight')
                plt.close('all')
