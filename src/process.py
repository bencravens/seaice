"""script file that processes the data imported from NETCDF4 files"""

#import libraries
import grab
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy
from scipy import stats
import os
import sys

def t_test_gridpoint(lons,lats,modelvar,controlvar,pval_filter, t_or_p):
    """with data given, performs the student's t-test on each gridpoint of an area for the new model (modelval) and the control model(controlvar) and returns either the t-statistic where the corresponding pvalue is below pval filter or p value at each point depending on t_or_p (if true, returns tstat, if false returns pvalues)"""
    #first checking that modelvar and control var are the same shape
    assert controlvar.shape == modelvar.shape
    [z,x,y] = modelvar.shape
    
    # now making an empty array to chuck results into
    tstats = np.ma.zeros([x, y])
    pvals = np.ma.zeros([x, y])
    for i in xrange(x):
        for j in xrange(y):
            tstats[i, j], pvals[i, j] = scipy.stats.mstats.ttest_ind(
                controlvar[:, i, j], modelvar[:, i, j])

    tstats = np.ma.masked_where(pvals > pval_filter, tstats)
    if t_or_p:
        return tstats
    else:
        return pvals

