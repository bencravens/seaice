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

def save_lims(modelname,monthnum,varname,csvdir,myvar):
    """with a given variable myvar and modelname, saves the upper and lower limits for the plots corresponding to the given model and variable into a csv file modelname.csv in the given directory csvdir. This is so that later on month_map_anom can make plots with consistent scaling."""
    #first changing into correct directory
    os.chdir(csvdir)
    print os.listdir("./")
    
    #deleting old file if it exists
    if os.path.exists("{}_{}_{}.txt".format(modelname,monthnum,varname)):
        os.remove("{}_{}_{}.txt".format(modelname,monthnum,varname))
    
    #now making a blank file
    f  = open("{}_{}_{}.txt".format(modelname,monthnum,varname),"w+")
    f.close()
    
    #now making data dictionary to write to csv file
    data = {'Model': modelname, 'Month': monthnum, 'Variable': varname, 'Max': np.ma.max(myvar)*1.1, 'Min': 0.9*np.ma.min(myvar)}

    #now opening the made file with csv
    with open("{}_{}_{}.txt".format(modelname,monthnum,varname),"w") as csv_file:
        fieldnames = ['Model', 'Month', 'Variable', 'Max', 'Min']
        writer = csv.DictWriter(csv_file,fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(data)

    print "lims written to csv file"
