# Copyright (C) 2018 by
# Niklas Boers <boers@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

import numpy as np
from scipy import weave
import scipy.stats as st


tlen = 5844
mnoe = 270
nodes = 2
y = 16

index_seasons = np.zeros(4, dtype = 'object')
index_seasons[0] = np.arange(0, 59, 1)
index_seasons[1] = np.arange(59, 151, 1)
index_seasons[2] = np.arange(151, 243, 1)
index_seasons[3] = np.arange(243, 334, 1)
for i in range(0, y - 1):
	index_seasons[0] = np.concatenate((index_seasons[0], np.arange(334 + 365 * i, 424 + 365 * i, 1)))
	index_seasons[1] = np.concatenate((index_seasons[1], np.arange(424 + 365 * i, 516 + 365 * i, 1)))
	index_seasons[2] = np.concatenate((index_seasons[2], np.arange(516 + 365 * i, 608 + 365 * i, 1)))
	index_seasons[3] = np.concatenate((index_seasons[3], np.arange(608 + 365 * i, 699 + 365 * i, 1)))


def EvSync(dat, tlen, nodes, core, noc, tm):
	"""Calculates triangle of Event Synchronization Matrix Q as list"""
	# tlen = length of time series
	# workspace (tmp event series)
	ex = np.empty(tlen, dtype='int32')
	ey = np.empty(tlen, dtype='int32')
	# output
	l = nodes * ((nodes / noc) - 1) / 2 + (core )  * (nodes / noc)
	#l = .5 * nodes * (nodes - 1)
	Q = np.zeros(l, dtype='uint16')
	# delay tau < taumax
	taumax = tm
	#parameters for parallelization
	var = ['dat', 'Q', 'nodes', 'tlen', 'taumax', 'core', 'noc','ex','ey']
	src = r"""
	int tau, tmp, count, dst;
	int i, k, m, n, t, c;
	int lx, ly;
	float q;
	/*ex = malloc(tlen*sizeof(int));
	ey = malloc(tlen*sizeof(int));*/
	/* generate event series ex and ey */
	c = 0;
	for(i = core; i < nodes; i += noc) {
	    lx = 0;
	    for(t = 0; t < tlen; t++)
	        if(dat[i*tlen + t])
	            ex[lx++] = t;
	    for(k = 0; k < i; k++) {
	        ly = 0;
	        for(t = 0; t < tlen; t++)
	            if(dat[k*tlen + t])
	                 ey[ly++] = t;
	        /* count event synchronisations */
	        count = 0;
	        for(m = 1; m < lx - 1; m++) {
	            for(n = 1; n < ly - 1; n++) {
				dst = ex[m] -ey[n];
				if(dst > taumax)
					continue;
	                tmp = ex[m+1] - ex[m];
	                if(tmp > ex[m] - ex[m-1])
	                    tmp = ex[m] - ex[m-1];
	                tau = ey[n+1] - ey[n];
	                if(tau > ey[n] - ey[n-1])
	                    tau = ey[n] - ey[n-1];
	                if(tau > tmp)
	                    tau = tmp;
	                tau /= 2;
	                if(abs(ex[m] - ey[n]) <= taumax && abs(ex[m] - ey[n]) < tau)
	                 count++;
					if(dst < -taumax)
						break;
	            }
	        }
			if(lx<3)
			{
				q = 0;
			}
			else if(ly<3)
			{
				q = 0;
			}
			else
			{
	        	q = count;
			}
	        Q[c] = q;
			c += 1;
	    }
	}
	"""
	weave.inline(src, var)
	del dat, ex, ey
	return Q

# P = np.zeros((7875,3)) --> xrange(3, 128)
# 270 events max:

# P1 = np.zeros((36046,3))
# P2 = np.zeros((36046,3))
# P3 = np.zeros((36046,3))
#
# for tm in xrange(3,4):
#    a = 0
#    for i in xrange(3,271):
#       for j in xrange(3,i + 1):
#          l = index_seasons[0].shape
#          season1 = np.zeros(l, dtype = "bool")
#          season2 = np.zeros(l, dtype = "bool")
#          index_seasons[0] = np.array(index_seasons[0], dtype = 'uint32')
#          season1[:i] = 1
#          season2[:j] = 1
#          dat = np.zeros((nodes, tlen), dtype = "bool")
#          cor = np.zeros(1000)
#          for k in xrange(1000):
#             dat[0, index_seasons[0]] = np.random.permutation(season1)
#             dat[1, index_seasons[0]] = np.random.permutation(season2)
#             core = 0
#             noc = 1
#             cor[k] = EvSync(dat, tlen, nodes, core, noc, tm)
#          th05 = st.scoreatpercentile(cor, 95)
#          th02 = st.scoreatpercentile(cor, 98)
#          th01 = st.scoreatpercentile(cor, 99)
#          P1[a] = [i, j, th05]
#          P2[a] = [i, j, th02]
#          P3[a] = [i, j, th01]
#          a += 1
#    np.save('trmm7_P_1000_3ms_mnoe270_thresholds_05_tm%d'%tm,P1)
#    np.save('trmm7_P_1000_3ms_mnoe270_thresholds_02_tm%d'%tm,P2)
#    np.save('trmm7_P_1000_3ms_mnoe270_thresholds_01_tm%d'%tm,P3)

P1 = np.zeros((226801,3))
P2 = np.zeros((226801,3))
P3 = np.zeros((226801,3))
P4 = np.zeros((226801,3))
P5 = np.zeros((226801,3))


for tm in [30]:
   a = 0
   for i in xrange(3,676):
      for j in xrange(3,i + 1):
         l = index_seasons[2].shape
         season1 = np.zeros(l, dtype = "bool")
         season2 = np.zeros(l, dtype = "bool")
         index_seas = np.array(index_seasons[2], dtype = 'uint32')
         season1[:i] = 1
         season2[:j] = 1
         dat = np.zeros((nodes, tlen), dtype = "bool")
         cor = np.zeros(2000)
         for k in xrange(2000):
            dat[0, index_seas] = np.random.permutation(season1)
            dat[1, index_seas] = np.random.permutation(season2)
            core = 0
            noc = 1
            cor[k] = EvSync(dat, tlen, nodes, core, noc, tm)
         th05 = st.scoreatpercentile(cor, 95)
         th02 = st.scoreatpercentile(cor, 98)
         th01 = st.scoreatpercentile(cor, 99)
         th005 = st.scoreatpercentile(cor, 99.5)
         th001 = st.scoreatpercentile(cor, 99.9)
         P1[a] = [i, j, th05]
         P2[a] = [i, j, th02]
         P3[a] = [i, j, th01]
         P4[a] = [i, j, th005]
         P5[a] = [i, j, th001]
         a += 1
   np.save('trmm7_P_2000_3ms_mnoe675_thresholds_05_tm%d_2'%tm, P1)
   np.save('trmm7_P_2000_3ms_mnoe675_thresholds_02_tm%d_2'%tm, P2)
   np.save('trmm7_P_2000_3ms_mnoe675_thresholds_01_tm%d_2'%tm, P3)
   np.save('trmm7_P_2000_3ms_mnoe675_thresholds_005_tm%d_2'%tm, P4)
   np.save('trmm7_P_2000_3ms_mnoe675_thresholds_001_tm%d_2'%tm, P5)
