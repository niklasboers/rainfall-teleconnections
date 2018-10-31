# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# BSD license.



import numpy as np
import matplotlib.pyplot as plt
from basemap_fct import basemap


data = 'trmm7'
lat = np.load('data/%s_lat_global.npy'%data)
lon = np.load('data/%s_lon_global.npy'%data)

if data == 'gpcp':
   origin = 180
   shift = 180
elif data == 'trmm7':
   origin = 0
   shift = 0


la = len(lat)
lo = len(lon)
nodes = la * lo
m = nodes * (nodes - 1) / 2
for tm in [10]:
   print "tm ", tm
   for j in [2]:
      print "season ", j
      for perc in np.arange(95, 96):
         print "perc ", perc
         th = np.load('data/%s_global_wd_score3_seasonal_rain_perc%d_season%d.npy'%(data, perc, j))
         noe = np.load('data/%s_global_wd_noe3_seasonal_rain_perc%d_season%d.npy'%(data, perc, j))
         nob = np.load('data/%s_global_wd_nob_cor_seasonal_rain_perc%d_season%d.npy'%(data, perc, j))
         area = np.where(nob.reshape((la, lo)) < 3, 1, -1)
         deg = np.load('data/degree01_%s_global_wd_perc%d_tm%d_season%d_2.npy'%(data, perc, tm, j))
         print np.mean(deg[noe < 3])

         fig = plt.figure(figsize = (8,5))
         fig.add_subplot(211)
         basemap(th.reshape((la, lo)), lat, lon, proj = 'mill', origin = origin, shift = shift, hatch = area, hatch_color = 'w', hillshade = False, alpha = .7, contours = np.linspace(0, 100, 11), extend = 'max', color = plt.cm.viridis_r, meridians = np.arange(0,360, 60), parallels = np.arange(-90,90,30), cbar_title = '95th percentile of rainfall [mm/day]')
         fig.add_subplot(212)
         basemap(nob.reshape((la, lo)), lat, lon, proj = 'mill', origin = origin, shift = shift, hatch = area, hatch_color = 'w', hillshade = False, alpha = .7, contours = np.linspace(0, 70, 8), extend = 'max', color = plt.cm.viridis_r, meridians = np.arange(0,360, 60), parallels = np.arange(-90,90,30), cbar_title = '# EREs')
         plt.savefig('pics/%s/clim/%s_01_global_wd_clim_perc%d_tm%d_season%d_2_cor_mill_nob'%(data, data, perc, tm, j), dpi = 400, bbox_inches = 'tight')
         plt.close()
