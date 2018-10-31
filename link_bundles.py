# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

def get_index_from_coord(lats, lons, lat_t, lon_t):
   from itertools import product
   lat = lats[np.argmin(np.abs(lats - lat_t))]
   lon = lons[np.argmin(np.abs(lons - lon_t))]
   coords = np.array(list(product(lats, lons)))
   idx = np.intersect1d(np.where(coords[:, 0] == lat), np.where(coords[:,1] == lon))
   return idx

def shperical_kde(values, xy, bw_opt):
   kde = KernelDensity(bandwidth=bw_opt, metric='haversine', kernel='gaussian', algorithm='ball_tree')
   kde.fit(values)
   datss = np.exp(kde.score_samples(xy))
   return datss





regions = ['nism']
prefix = ''
prefix2 = '/p/tmp/boers/'

data = 'trmm7'
print data
lat = np.load(prefix + '%s_lat_global.npy'%data)
lon = np.load(prefix + '%s_lon_global.npy'%data)

for perc in [80,85, 90, 94, 95, 96]:

   for tm in [3,10,30]:
      print perc
      print "tm = ", tm
      for region in regions:
         print "region: ", region

         if region == 'sam2':
            j = 0
         else:
            j = 2

         events = np.load(prefix + '%s_global_wd_nob_cor_seasonal_rain_perc%d_season%d.npy'%(data, perc, j))
         deg = np.load(prefix + 'degree01_%s_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(data, perc, tm, j))
         cu_deg = np.load(prefix + 'cu_degree01_%s_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(data, perc, tm, j))
         nlist = np.load(prefix2 + 'nlist01_%s_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(data, perc, tm, j))
         ang_dist_all = np.load('/p/tmp/boers/ang_dist01_%s_global_wd_perc%d_tm%d_season%d_2_nb_sig005.npy'%(data, perc, tm, j))

         links = np.array([], dtype = 'int')
         gdists = np.array([])
         index = np.array([], dtype = 'int')


         for lat_t in [27.5, 27.75, 28., 28.25, 28.5]:
            for lon_t in [78.5, 78.75, 79., 79.25, 79.5]:
               index = np.concatenate((index, get_index_from_coord(lat, lon, lat_t, lon_t)))

         lat_t = 28
         lon_t = 79

         index = np.unique(index)

         for l in xrange(index.shape[0]):
            links = np.concatenate((links, nlist[cu_deg[index[l]] : cu_deg[index[l] + 1]]))
            gdists = np.concatenate((gdists, ang_dist_all[cu_deg[index[l]] : cu_deg[index[l] + 1]]))
         print "nol = ", links.shape[0]
         selector = np.zeros(links.shape[0], dtype = 'int')
         selector[:10000] = 1
         selector = np.random.permutation(selector)
         links = links[selector == 1]
         gdists = gdists[selector == 1]

         np.save('%s_geo_links_samples_perc%d_tm%d_season%d_%s_nb_sig005.npy'%(data, perc, tm, j, region), links)
         np.save('%s_geo_gdists_samples_perc%d_tm%d_season%d_%s_nb_sig005.npy'%(data, perc, tm, j, region), gdists)

         print "nol = ", links.shape[0]
         index = get_index_from_coord(lat, lon, lat_t, lon_t)

         from itertools import product
         coords = np.array(list(product(lat, lon)))
         lat0 = coords[index, 0]
         lon0 = coords[index, 1]
         lat1 = coords[links, 0]
         lon1 = coords[links, 1]
         nol = links.shape[0]

         values = np.vstack([coords[links, 0], coords[links, 1]]).T
         values *= np.pi / 180.
         X, Y = np.meshgrid(lon, lat)
         xy = np.vstack([Y.ravel(), X.ravel()]).T
         xy *= np.pi / 180.

         bw_opt = .2 * values.shape[0]**(-1./(2+4))
         kde = KernelDensity(bandwidth=bw_opt, metric='haversine', kernel='gaussian', algorithm='ball_tree')
         kde.fit(values)
         print "kde fitted"
         dat = np.exp(kde.score_samples(xy))
         print "kde evaluated"
         dat = dat.reshape(X.shape)
         np.save(prefix + '%s_geo_link_tc_sphkde_perc%d_tm%d_season%d_%s_bw02_nb_sig005'%(data, perc, tm, j, region), dat)
