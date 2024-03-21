cimport numpy as np
cimport cython
import numpy as np
from libc.stdlib cimport rand, RAND_MAX
# from numpy cimport ndarray

@cython.boundscheck(False)
def ESreg(np.ndarray[np.float_t, ndim = 2] e1, np.ndarray[np.float_t, ndim = 2] e2, np.ndarray[np.float_t, ndim = 1] t, np.ndarray[np.float_t, ndim = 1] t12, np.ndarray[np.float_t, ndim = 1] t21, int taumax, double tlen):
    cdef int i, j, m, n, n1, n2
    cdef int delay12 = 0
    cdef int delay21 = 0
    mnoe = e1.shape[1]
    n1 = e1.shape[0]
    n2 = e2.shape[0]
    for i in xrange(n1):
        for j in xrange(n2):
            for m in xrange(1, mnoe - 1):
                if e1[i, m] > 0:
                    for n in xrange(1, mnoe - 1):
                        if e2[j, n] > 0:
                            dst = e1[i, m] - e2[j, n]
                            if dst > taumax:
                                continue
                            tau = min([e1[i, m] - e1[i, m - 1], e1[i, m + 1] - e1[i, m], e2[j, n] - e2[j, n - 1], e2[j, n + 1] - e2[j, n]]) / 2.
                            if dst < 0 and abs(dst) < tau and abs(dst) < taumax:
                                t12[int(e1[i, m])] += 1
                                delay12 += abs(dst)
                            elif dst == 0:
                                t[int(e1[i, m])] += 1
                            elif dst > 0 and abs(dst) < tau and abs(dst) < taumax:
                                t21[int(e2[j, n])] += 1
                                delay21 += abs(dst)
                            if dst < -taumax:
                                break
    return t, t12, t21, delay12, delay21
