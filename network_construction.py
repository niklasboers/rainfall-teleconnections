# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.



import numpy as np
from scipy import weave
import scipy.stats as st


def geographic_link_dist_dist(lat, lon, deg, cu_deg, nlist):
   la = len(lat)
   lo = len(lon)
   n = la * lo
   lat_seq = np.repeat(lat, lo)
   lon_seq = np.tile(lon, la)
   sin_lon = np.sin(lon_seq * np.pi / 180)
   cos_lon = np.cos(lon_seq * np.pi / 180)
   sin_lat = np.sin(lat_seq * np.pi / 180)
   cos_lat = np.cos(lat_seq * np.pi / 180)
   ang_dist_sum = np.zeros(n)
   median_ang_dist = np.zeros(n)
   mean_ang_dist = np.zeros(n)
   ang_dist_5p = np.zeros(n)
   ang_dist_95p = np.zeros(n)
   ang_dist = np.zeros(len(nlist), dtype = 'uint16')
   var = ['cos_lat', 'sin_lat', 'cos_lon', 'sin_lon', 'ang_dist', 'ang_dist_sum', 'deg', 'cu_deg', 'nlist', 'n']
   code = r"""
   long long i, j, u;
   float expr;
   for (i = 0; i < n; i++) {
    	for (j = cu_deg[i]; j < cu_deg[i + 1]; j++) {
   		u = nlist[j];
        	expr = sin_lat[i]*sin_lat[u] + cos_lat[i]*cos_lat[u] * (sin_lon[i]*sin_lon[u] + cos_lon[i]*cos_lon[u]);
   		if (expr > 1) {
   			expr = 1.;
   		}
   		else if (expr < -1) {
            	expr = -1.;
   		}
   		ang_dist[j] = (int) (acos(expr) * 6371);
   	}
   }
   """
   weave.inline(code, var)
   return ang_dist


noc = 15


lat = np.load('trmm7_lat_global.npy')
lon = np.load('trmm7_lon_global.npy')
la = len(lat)
lo = len(lon)
nodes = la * lo
m = nodes * (nodes - 1) / 2

for tm in [10]:
   for j in [2]:
      for perc in [83]:
         deg = np.zeros(nodes, dtype = 'uint32')
         for core in xrange(noc):
            l = nodes * ((nodes / noc) - 1) / 2 + core  * (nodes / noc)
            print l
            e = 0
            Q_slice = np.load('/p/tmp/boers/Q_trmm7_global_wd_perc%d_tm%d_season%d_bool_nb_sig005_slice%d.npy'%(perc, tm, j, core))
            print "Q_slice loaded, shape:", Q_slice.shape[0]
            print nodes
            var = ['Q_slice', 'noc', 'core', 'nodes', 'deg']
            src = r"""
            long long i, k;
            long long e = 0;
            long long e1, e2;
            i = core;
            k = 0;
            for(i = core; i < nodes; i+=noc){
            	for(k = 0; k < i; k+=1){
                  if(Q_slice[e]){
                     /*printf("Q = %d\n", Q_slice[e]);
                     printf("P = %d\n", P[((e1 - 2) * (e1 - 3) / 2 + e2 - 3) * 3 + 2]);
                     printf("i, k = %d, %d\n", i, k);*/
                     deg[i]+=1;
                     deg[k]+=1;
                     /*printf("degi, degk = %d, %d\n", deg[i], deg[k]);*/
                  }
            		e++;
            	}
            }
            """
            weave.inline(src, var)
            del Q_slice
         print "deg done"
         np.save('degree01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005'%(perc, tm, j), deg)
         cu_deg = np.zeros(nodes + 1, dtype = 'uint64')
         var = ['deg', 'cu_deg', 'nodes']
         src = r"""
         long long i, s;
         s = 0;
         for(i = 0; i < nodes + 1; i++){
            cu_deg[i] = s;
            s += deg[i];
         }
         """
         weave.inline(src, var)
         np.save('cu_degree01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005'%(perc, tm, j), cu_deg)
         print "cu_deg done"
         positions = cu_deg.copy()
         nlist =  np.zeros(np.sum(deg), dtype = 'uint32')
         for core in xrange(noc):
            Q_slice = np.load('/p/tmp/boers/Q_trmm7_global_wd_perc%d_tm%d_season%d_bool_nb_sig005_slice%d.npy'%(perc, tm, j, core))
            var = ['Q_slice', 'core', 'noc', 'nodes', 'nlist', 'positions']
            src = r"""
            long long i, k, v, w;
            long long e = 0;
            long long e1, e2;
            for(i = core; i < nodes; i += noc){
               for(k = 0; k < i; k++){
                  v = positions[k];
                  w = positions[i];
                  if(Q_slice[e]){
                     nlist[v] = i;
                     nlist[w] = k;
                     positions[k]++;
                     positions[i]++;
                  }
                  e += 1;
               }
            }
            """
            weave.inline(src, var)
            del Q_slice
         print "nlist done"
         ld = np.sum(deg) / float(nodes * (nodes - 1))
         print "link density for global and perc %d for tm %d in season %d is %f"%(perc, tm, j, ld)
         print "nol = ", nlist.shape[0]
         print "saving nlist ..."
         np.save('/p/tmp/boers/nlist01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005'%(perc, tm, j), nlist)
         print "calculating geographical link distances ..."
         ang_dist = geographic_link_dist_dist(lat, lon, deg, cu_deg, nlist)
         print "saving ..."
         np.save('/p/tmp/boers/ang_dist01_trmm7_global_wd_perc%d_tm%d_season%d_2_nb_sig005'%(perc, tm, j), ang_dist)
