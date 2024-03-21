import numpy as np
import cython_func

def ESreg(e1, e2, taumax, tlen):
    e1 = np.array(e1, dtype = 'int')
    e2 = np.array(e2, dtype = 'int')
    mnoe = e1.shape[1]
    n1 = e1.shape[0]
    n2 = e2.shape[0]
    t12 = np.zeros(tlen + 1)
    t21 = np.zeros(tlen + 1)
    t = np.zeros(tlen + 1)
    delay12 = 0
    delay21 = 0
    for i in xrange(n1):
        for j in xrange(n2):
            for m in xrange(1, mnoe - 1):
                if e1[i, m] > 0:
                    for n in xrange(1, mnoe - 1):
                        if e2[j, n] > 0:
                            dst = e1[i, m] - e2[j, n]
                            if dst > taumax:
                                continue
                            tau = np.min([e1[i, m] - e1[i, m - 1], e1[i, m + 1] - e1[i, m], e2[j, n] - e2[j, n - 1], e2[j, n + 1] - e2[j, n]]) / 2.
                            if dst < 0 and np.abs(dst) < tau and np.abs(dst) < taumax:
                                t12[e1[i, m]] += 1
                                delay12 += np.abs(dst)
                            elif dst == 0:
                                t[e1[i, m]] += 1
                            elif dst > 0 and np.abs(dst) < tau and np.abs(dst) < taumax:
                                t21[e2[j, n]] += 1
                                delay21 += np.abs(dst)
                            if dst < -taumax:
                                break
    return t, t12, t21, delay12 / np.sum(t12), delay21 / np.sum(t21)



n1 = 2
n2 = 2

mnoe = 70
tlen = 6000
e1 = np.zeros((n1, mnoe))
e2 = np.zeros((n2, mnoe))


for i in xrange(n1):
    a = np.sort(np.unique((np.random.randint(0, tlen, mnoe))))
    b = np.sort(np.unique((np.random.randint(0, tlen, mnoe))))
    e1[i, :a.shape[0]] = a
    # e2[i, :b.shape[0]] = b
    e2[i, :a.shape[0]] = a + 5

print e1
print e2
taumax = 10
# t, t12, t21, delay12, delay21 =  ESreg(e1, e2, taumax, tlen)
#
#
# print t12.sum()
# print delay12
#
#
# print t21.sum()
# print delay21


t = np.zeros(tlen)
t12 = np.zeros(tlen)
t21 = np.zeros(tlen)

t, t12, t21, delay12, delay21 = cython_func.ESreg(e1, e2, t, t12, t21, taumax, tlen)

print t12.sum()
print delay12 / t12.sum()


print t21.sum()



print delay21 / t21.sum()
