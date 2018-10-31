# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.


import numpy as np
import matplotlib.pyplot as plt
from basemap_fct import basemap
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy.signal import lombscargle, welch

def coordinate_indices_from_ra(lat, lon, lat_max, lat_min, lon_max, lon_min):
    la = len(lat)
    lo = len(lon)
    n = la * lo
    indices = np.arange(n).reshape((la, lo))
    lat_max = lat[np.argmin(np.abs(lat - lat_max))]
    lat_min = lat[np.argmin(np.abs(lat - lat_min))]
    lon_max = lon[np.argmin(np.abs(lon - lon_max))]
    lon_min = lon[np.argmin(np.abs(lon - lon_min))]
    ra_indices = indices[np.where(lat == lat_max)[0][0] : np.where(lat == lat_min)[0][0] + 1 , np.where(lon == lon_min)[0][0] : np.where(lon == lon_max)[0][0] + 1 ]
    return np.unique(ra_indices.flatten())


lat = np.load('data/ncep_lat.npy')

# lon = np.load('data/ncep_lon.npy') -180.
lon = np.load('data/ncep_lon.npy')

la = len(lat)
lo = len(lon)
nodes = la * lo

levels = ['925hPa', '500hPa', '250hPa']
seasons = ['DJF', 'MAM', 'JJA', 'SON']

monsoons = ['ISM1']
syncs = ['ISM1-EUR']
sync_pairs = [[0, 1]]

sync2 = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2]

sync_idx = [[0, 1]]


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
shift = [-14, -12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14]
contlevel1 = [5., 20., 30., 40.]
contlevel2 = [5., 2*1e6, 3*1e6, 4*1e6]
contlevel3 = [2., 3., 4.]
contlevel31 = [4., 5., 8.]
contlevel32 = [15., 20., 30.]
contlevel41 = [1., 2., 3.]
contlevel42 = [4., 5., 8.]


m = 0
for season in ['jja']:
    for sth in [90]:
        for cutoff in [8, 10, 12]:
            for perc in xrange(95, 96):


                for direct in ['21']:
                    sync_ano = np.load('data/sync_hgt_wind_anomalies_%s_perc%d_lp%d_sth%d_mjo%d_%s_nb.npy'%(season, perc, cutoff, sth, m, direct))
                    sync_ano[np.isnan(sync_ano)] = 0

                    for k in [2]:

                        areaA = np.ones(nodes)
                        areaA[monsoon_indices[sync_idx[k][0]]] = -1
                        areaA = areaA.reshape((lat.shape[0], lon.shape[0]))
                        areaB = np.ones(nodes)
                        areaB[monsoon_indices[sync_idx[k][1]]] = -1
                        areaB = areaB.reshape((lat.shape[0], lon.shape[0]))

                        for l in xrange(len(shift)):
                            fig = plt.figure(figsize = (7.5, 5 ))
                            ax = fig.add_subplot(4, 2, 1)
                            prwt = sync_ano[k, l, 3, 0]
                            cont = contlevel1[0]
                            basemap(prwt, lat, lon, lat_min = -10, lat_max = 70, lon_min = -80, lon_max = 180, res = 'c', area1 = areaA, line_color1 = 'm', area2 = areaB, line_color2 = 'm', lms = 2., proj = 'mill', shift = 70., contours = [-5, -3, -1, 1, 3, 5], color = plt.cm.RdBu_r, alpha = .7, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-80,80,40), cbar_title = 'PRWT [kg/m^2]')

                            ax = fig.add_subplot(4, 2, 2)
                            rain = sync_ano[k, l, 4, 0]
                            rain[lat > 50.] = 0.
                            rain[lat < -50.] = 0.
                            cont = contlevel2[0]
                            basemap(rain, lat, lon, lat_min = -10, lat_max = 70, lon_min = -80, lon_max = 180, res = 'c', area1 = areaA, line_color1 = 'm', area2 = areaB, line_color2 = 'm', lms = 2., proj = 'mill', shift = 70., contours = [-5, -3, -1, 1, 3, 5], color = plt.cm.RdBu_r, alpha = .7, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-80,80,40), cbar_title = 'Rainfall [mm/d]')


                            for i in xrange(3):
                                ax = fig.add_subplot(4, 2, 2 * i + 3)
                                gph = sync_ano[k, l, 0, 2 - i]
                                za_gph = np.mean(gph, axis = 1)
                                gph = gph - za_gph.reshape((gph.shape[0], 1))
                                gphref = sync_ano[k, 2, 0, i]
                                za_gphref = np.mean(gphref, axis = 1)
                                gphref = gphref - za_gphref.reshape((gphref.shape[0], 1))
                                uwind = sync_ano[k, l, 1, 2 - i]
                                vwind = sync_ano[k, l, 2, 2 - i]

                                cont = contlevel1[2 - i + 1]

                                basemap(gph, lat, lon, lat_min = -10, lat_max = 70, lon_min = -80, lon_max = 180, res = 'c', area1 = areaA, line_color1 = 'm', area2 = areaB, line_color2 = 'm', lms = 2., proj = 'mill', shift = 70., uwinds = uwind, vwinds = vwind, asp = 4, quiverscale = 200, contours = [-25, -15, -5, 5, 15, 25], color = plt.cm.RdBu_r, alpha = .7, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-80,80,40), cbar_title = 'GPH %s [m]'%levels[2 - i])

                                w = VectorWind(uwind, vwind)
                                sf, vp = w.sfvp()
                                za_sf, za_vp = np.mean(sf, axis = 1), np.mean(vp, axis = 1)
                                sf, vp = sf - za_sf.reshape((sf.shape[0], 1)), vp - za_vp.reshape((sf.shape[0], 1))

                                sf_lat = np.mean(sf[int(np.where(lat == 47.5)[0]) : int(np.where(lat == 37.5)[0]), :], axis = 0)
                                f, Pxx = welch(sf_lat, 1. / 2.5)
                                wno = f[np.argmax(Pxx)] * 360.

                                ax = fig.add_subplot(4, 2, 2 * i + 4)
                                cont = contlevel2[i + 1]
                                basemap(sf * 1e-5, lat, lon, lat_min = -10, lat_max = 70, lon_min = -80, lon_max = 180, res = 'c', area1 = areaA, line_color1 = 'm', area2 = areaB, line_color2 = 'm', lms = 2., proj = 'mill', shift = 70., uwinds = uwind, vwinds = vwind, asp = 4, quiverscale = 200, contours = [-25, -15, -5, 5, 15, 25], color = plt.cm.RdBu_r, alpha = .7, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-80,80,40), cbar_title = 'PSI %s [m^2/s]'%levels[2 - i])



                            fig.suptitle('Anomalies %s %s Day %d Cutoff = %d STH = %d | k = %d'%(syncs[k], season, shift[l], cutoff, sth, wno))
                            fig.savefig('pics/composites/%s/gph/gph_zar_day%d_monsoon_%s_ano_lp%d_%s_sth%d_mjo%d_%s_nb.pdf'%(season, shift[l], syncs[k], cutoff, season, sth, m, direct), bbox_inches = 'tight')
                            plt.close('all')
                        for l in xrange(len(shift)):
                            fig = plt.figure(figsize = (5, 3))
                            ax = fig.add_subplot(2, 1, 1)

                            rain = sync_ano[k, l, 4, 0]
                            rain[lat > 50.] = 0.
                            rain[lat < -50.] = 0.
                            cont = contlevel2[0]
                            basemap(rain, lat, lon, lat_min = -10, lat_max = 70, lon_min = -80, lon_max = 180, res = 'c', area1 = areaA, line_color1 = 'm', area2 = areaB, line_color2 = 'm', lms = 2., proj = 'mill', shift = 70., contours = [-5, -3, -1, 1, 3, 5], color = plt.cm.RdBu_r, alpha = .7, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-80,80,40), cbar_title = 'Rainfall [mm/day]')

                            i = 2
                            gph = sync_ano[k, l, 0, i]
                            uwind = sync_ano[k, l, 1, i]
                            vwind = sync_ano[k, l, 2, i]


                            ax = fig.add_subplot(2, 1, 2)

                            cont = contlevel3[i]
                            basemap(vwind, lat, lon, lat_min = -10, lat_max = 70, lon_min = -80, lon_max = 180, res = 'c', area1 = areaA, line_color1 = 'm', area2 = areaB, line_color2 = 'm', lms = 2., proj = 'mill', shift = 70., contours = [-4, -3, -2, -1, 1, 2, 3, 4], color = plt.cm.RdBu_r, alpha = .7, colorbar = True, extend = 'both', meridians = np.arange(0,360, 40), parallels = np.arange(-80,80,40), cbar_title = 'v %s [m/s]'%levels[i])

                            vw_lat = np.mean(vwind[int(np.where(lat == 47.5)[0]) : int(np.where(lat == 37.5)[0]), :], axis = 0)
                            f, Pxx = welch(vw_lat, 1. / 2.5)
                            wno = f[np.argmax(Pxx)] * 360.



                            fig.suptitle('Anomalies %s %s Day %d Cutoff = %d STH = %d | k = %d'%(season, syncs[k], shift[l], cutoff, sth, wno))
                            fig.savefig('pics/composites/%s/vwind/vwind2_zar_day%d_monsoon_%s_ano_lp%d_%s_sth%d_mjo%d_%s_nb.pdf'%(season, shift[l], syncs[k], cutoff, season, sth, m, direct), bbox_inches = 'tight')
                            plt.close('all')
