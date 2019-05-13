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

def read_lims(varname,csvdir):
    """reads upper and lower limits for plot from the dir csvdir in the file varname.txt and returns them"""
    #first changing into correct directory
    os.chdir(csvdir)
    print os.listdir("./")    
    
    with open("{}.txt".format(varname)) as File:
        reader = csv.DictReader(File)
        for row in reader:
            return row #there should only be one row

def anom_limit_setup(varname,months,models,csvdir):
    """sets upper and lower bound for a variable for anomaly plots based on the max/min value TOTAL across all anomalies"""
    for model in models:
            maxes = []
            mins = []
            for month in months:
                lons,lats,myvar,totaldiff,units = grab.month_map_anom("/media/windowsshare",model,month,varname,False)
                maxes.append(np.ma.max(myvar))
                mins.append(np.ma.min(myvar))
    os.chdir(csvdir)
    with open("{}.txt".format(varname),"w+") as csv_file:
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

def total_ice_diff(path,modelname,monthnum):
    """function to calculate the total difference in ice volume between a model and a control for a give month"""  
    lons, lats, aice_diff, total_aice_diff, units = grab.month_map_anom(path,modelname,monthnum,"aice",True)
    lons, lats, sithick_diff, total_sithick_diff, units = grab.month_map_anom(path,modelname,monthnum,"sithick",True)
    print "total aice change is {}".format(np.sum(aice_diff))
    print "shape of aice_diff is {}".format(np.shape(aice_diff))
    print "shape of sithick diff is {}".format(np.shape(sithick_diff))
    print "total ice volume difference is {}".format(np.sum(np.multiply(aice_diff,sithick_diff)))
    

if __name__=="__main__":
    total_ice_diff("/media/windowsshare/","u-au866",2) 
