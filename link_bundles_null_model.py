# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.




import numpy as np
import scipy.stats as st
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
import mpi

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



def master():
   regions = ['nism']
   prefix = ''
   prefix2 = '/p/tmp/boers/'

   data = 'trmm7'

   lat = np.load(prefix + '%s_lat_global.npy'%data)
   lon = np.load(prefix + '%s_lon_global.npy'%data)

   c = 0
   for perc in [80, 85, 90, 94, 95, 96]:
      print perc
      for tm in [3,10,30]:
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

            links = np.array([], dtype = 'int')
            index = np.array([], dtype = 'int')

            for lat_t in [27.5, 27.75, 28., 28.25, 28.5]:
                for lon_t in [78.5, 78.75, 79., 79.25, 79.5]:
                    index = np.concatenate((index, get_index_from_coord(lat, lon, lat_t, lon_t)))

            lat_t = 28
            lon_t = 79

            index = np.unique(index)

            links = np.load('%s_geo_links_samples_perc%d_tm%d_season%d_%s_nb_sig005.npy'%(data, perc, tm, j, region))

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
            print ".2 * bw_scott = ", bw_opt

            nos = 15
            ser = 35
            dats = np.zeros((nos * ser, X.shape[0], X.shape[1]))
            for k in xrange(ser):
               for i in xrange(nos):
                  noel2 = np.where(events > 2)[0]
                  values = np.vstack([np.random.choice(coords[noel2, 0], nol), np.random.choice(coords[noel2, 1], nol)]).T
                  values *= np.pi / 180.
                  mpi.submit_call("shperical_kde", (values, xy, bw_opt), id = c * nos * ser + k * nos + i)
                  print "bash %d submitted"%(k * nos + i)
               for i in xrange(nos):
                  datss = mpi.get_result(id = c * nos * ser + k * nos + i)
                  dats[k * nos + i] = datss.reshape(X.shape)
                  print "bash %d picked up"%(k * nos + i)
            np.save(prefix2 + '%s_geo_link_tc_sphkde_nm_perc%d_tm%d_season%d_%s_bw02_nb'%(data, perc, tm, j, region), dats)
            mean = np.mean(dats, axis = 0)
            std = np.std(dats, axis = 0)
            perc90 = st.scoreatpercentile(dats, 90, axis = 0)
            perc95 = st.scoreatpercentile(dats, 95, axis = 0)
            perc99 = st.scoreatpercentile(dats, 99, axis = 0)
            perc995 = st.scoreatpercentile(dats, 99.5, axis = 0)
            perc999 = st.scoreatpercentile(dats, 99.9, axis = 0)

            stats = np.array([mean, std, perc90, perc95, perc99, perc995, perc999])
            np.save('%s_geo_link_tc_sphkde_stats_perc%d_tm%d_season%d_%s_bw02_nb.npy'%(data, perc, tm, j, region), stats)
            c += 1

mpi.run()
