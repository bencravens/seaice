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
import csv

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

def read_lims(varname,month,csvdir):
    """reads upper and lower limits for plot from the dir csvdir in the file varname.txt and returns them"""
    #first changing into correct directory
    os.chdir(csvdir)
    print os.listdir("./")    
    
    with open("{}_{}.txt".format(varname,month)) as File:
        reader = csv.DictReader(File)
        for row in reader:
            return row #there should only be one row

def anom_limit_setup(varname,month,models,csvdir):
    """sets upper and lower bound for a variable for anomaly plots based on the max/min value TOTAL across all anomalies"""
    for model in models:
        lons,lats,myvar,totaldiff,units = grab.month_map_anom("/media/windowsshare",model,month,varname,False)
        maxes = []
        mins = []
        maxes.append(np.ma.max(myvar))
        mins.append(np.ma.min(myvar))
    os.chdir(csvdir)
    with open("{}_{}.txt".format(varname,month),"w+") as csv_file:
        fieldnames=['Max','Min']
        writer = csv.DictWriter(csv_file,fieldnames=fieldnames) 
        writer.writeheader()
        absmax = abs(max(maxes))
        absmin = abs(min(mins))
        if absmax > absmin:
            data = {'Max':absmax, 'Min':-1.0*absmax}
        elif absmin > absmax:
            data = {'Max':absmin, 'Min':-1.0*absmin}
        else:
            data = {'Max':absmin, 'Min':-1.0*absmin}
        writer.writerow(data)
    print "finished writing limits. They are {} and {}".format(max(maxes),min(mins))










