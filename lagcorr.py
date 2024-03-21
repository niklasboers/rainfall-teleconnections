import numpy as np
import scipy.stats as st


def lagcorr(x,y,lag=None,verbose=False, cor = 'P'):
    '''Compute lead-lag correlations between 2 time series.

    <x>,<y>: 1-D time series.
    <lag>: lag option, could take different forms of <lag>:
          if 0 or None, compute ordinary correlation and p-value;
          if positive integer, compute lagged correlation with lag
          upto <lag>;
          if negative integer, compute lead correlation with lead
          upto <-lag>;
          if pass in an list or tuple or array of integers, compute
          lead/lag correlations at different leads/lags.

    Note: when talking about lead/lag, uses <y> as a reference.
    Therefore positive lag means <x> lags <y> by <lag>, computation is
    done by shifting <x> to the left hand side by <lag> with respect to
    <y>.
    Similarly negative lag means <x> leads <y> by <lag>, computation is
    done by shifting <x> to the right hand side by <lag> with respect to
    <y>.

    Return <result>: a (n*2) array, with 1st column the correlation
    coefficients, 2nd column correpsonding p values.

    Currently only works for 1-D arrays.
    '''

    import numpy
    from scipy.stats import pearsonr, spearmanr

    if len(x)!=len(y):
        raise('Input variables of different lengths.')

    #--------Unify types of <lag>-------------
    if numpy.isscalar(lag):
        if abs(lag)>=len(x):
            raise('Maximum lag equal or larger than array.')
        if lag<0:
            lag=-numpy.arange(abs(lag)+1)
        elif lag==0:
            lag=[0,]
        else:
            lag=numpy.arange(lag+1)
    elif lag is None:
        lag=[0,]
    else:
        lag=numpy.asarray(lag)

    #-------Loop over lags---------------------
    result=[]
    if verbose:
        print '\n#<lagcorr>: Computing lagged-correlations at lags:',lag
    if cor == 'P':
        for ii in lag:
            if ii<0:
                result.append(pearsonr(x[:ii],y[-ii:]))
            elif ii==0:
                result.append(pearsonr(x,y))
            elif ii>0:
                result.append(pearsonr(x[ii:],y[:-ii]))
    elif cor == 'S':
        for ii in lag:
            if ii<0:
                result.append(spearmanr(x[:ii],y[-ii:]))
            elif ii==0:
                result.append(spearmanr(x,y))
            elif ii>0:
                result.append(spearmanr(x[ii:],y[:-ii]))

    result=numpy.asarray(result)

    return result
