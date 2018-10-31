# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.




import sys
import mpi
import numpy as np
from scipy import weave
import scipy.stats as st
import time

def EvSync(e, noe, nodes, core, noc, perc, tm, season, P, noep):
   print "batch %d running ..."%core
   l = nodes * ((nodes / noc) - 1) / 2 + core  * (nodes / noc)
   Q = np.zeros(l, dtype='bool')
   taumax = tm
   var = ['e', 'Q', 'nodes', 'noe', 'taumax', 'core', 'noc', 'P', 'noep']
   src = r"""
	int tau, tmp, count, dst;
	long i, k, m, n, t, c, li, lk;
   long long e1, e2;
	int q;
	c = 0;
	for(i = core; i < nodes; i += noc){
	    for(k = 0; k < i; k++){
	        count = 0;
           li = 0;
           lk = 0;
	        for(m = 1; m < noep[i] - 1; m++) {
				if(e[i * noe +  m] > 0){
					li++;
	            for(n = 1; n < noep[k] - 1; n++) {
						if(e[k * noe +  n] > 0){
							lk++;
							dst = e[i * noe + m] - e[k * noe +  n];
							if(dst > taumax)
								continue;
		                	tmp = e[i * noe + m + 1] - e[i * noe + m];
		                	if(tmp > e[i * noe + m] - e[i * noe + m - 1])
		                    	tmp = e[i * noe + m] - e[i * noe + m - 1];
		                	tau = e[k * noe +  n + 1] - e[k * noe +  n];
		                	if(tau > e[k * noe +  n] - e[k * noe +  n - 1])
		                    	tau = e[k * noe +  n] - e[k * noe +  n - 1];
		                	if(tau > tmp)
		                    	tau = tmp;
		                	tau /= 2;
	                		if(abs(e[i * noe + m] - e[k * noe + n]) <= taumax && abs(e[i * noe + m] - e[k * noe + n]) < tau)
	                 			count++;
							if(dst < -taumax)
							break;
						}
					}
	            }
	        }
			if(li<3)
			{
				q = 0;
			}
			else if(lk<3)
			{
				q = 0;
			}
			else
			{
	        	q = count;
			}
         e1 = noep[i];
         e2 = noep[k];
         if(e1 < e2){
            e1 = noep[k];
            e2 = noep[i];
         }
         if(e1 > 2 && e2 > 2 && q > P[((e1 - 2) * (e1 - 3) / 2 + e2 - 3) * 3 + 2]){
	         Q[c] = 1;
         }
         else{
            Q[c] = 0;
         }
			c += 1;
	    }
	}
	"""
   weave.inline(src, var)
   del e
   np.save('/p/tmp/boers/Q_trmm7_global_wd_perc%d_tm%d_season%d_bool_nb_sig005_slice%d'%(perc, tm, season, core), Q)
   print "batch %d saved ..."%core
   del Q
   return 0



def master():
   c = 0
   for perc in xrange(80, 100):
      print "perc", perc
      for tm in [10]:
         print "tm", tm
         for j in [2]:
            print "season", j
            dat = np.load('/p/tmp/boers/trmm7_global_wd_bursts_cor_seasonal_rain_perc%d_season%d.npy'%(perc, j))
            P = np.array(np.load('trmm7_P_2000_3ms_mnoe675_thresholds_005_tm%d_2.npy'%tm).flatten(), dtype = 'int')
            noep = np.load('trmm7_global_wd_nob_cor_seasonal_rain_perc%d_season%d.npy'%(perc, j))
            nodes = dat.shape[0]
            noe = dat.shape[1]
            dat = dat.reshape(nodes * noe)
            noc = 15
            print "calculating event synchronization ..."
            for core in range(noc):
               mpi.submit_call("EvSync",(dat, noe, nodes, core, noc, perc, tm, j, P, noep), id = c)
               print "batch %d submitted ..."%core
               c += 1
            del dat
mpi.run()
