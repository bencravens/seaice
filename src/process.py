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
from scipy.interpolate import griddata

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

def regrid(arr1,lats1,lons1,arr2,lats2,lons2,modelname,monthstr):
    """takes two netCDF arrays and their respective latitudes and longitudes and regrids 
    both arrays to be on the same grid."""
    #NOTE!!!!!!!! Here arr1 should be the model data and arr2 should be the NSIDC data (this is because the model has coarser resolution
    # in the area that we are plotting over. 
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawmapboundary(fill_color='grey')
    
    #projecting both arrays onto grid
    #make mesh grid of latlons for plotting
    ny=arr1.shape[0]; nx=arr1.shape[1]
    print("shape of model data is {},{}".format(nx,ny))
    glons,glats = m.makegrid(nx,ny)

    #interp data of real world to this reg grid 
    gdata1 = griddata((lons1[~arr1.mask].ravel(),lats1[~arr1.mask].ravel()),
                    arr1[~arr1.mask].ravel(),(glons,glats),method='cubic')
    
    gdata2 = griddata((lons2[~arr2.mask].ravel(),lats2[~arr2.mask].ravel()),arr2[~arr2.mask].ravel(),(glons,glats),method='cubic')
   
    # make mask of missing data (ie land)
    mask = np.zeros((ny,nx))
    mask[arr1.mask] = 1
    gmask = griddata((lons1.ravel(),lats1.ravel()),mask.ravel(),
                    (glons,glats),method='cubic')
    pdat = np.ma.masked_array(gdata1-gdata2,gmask>0.5) # masked, gridded array
    return gdata1,gdata2,pdat,glons,glats

def permutation_test(arraystack_1, arraystack_2):
    #takes two n*m*t arrays, where the n*m represents a grid of variables and the t represents that area of gridcells as it changes in time. 
    [x,y] = np.shape(arraystack_1[0])
    pvals = np.zeros([x,y])
    count = 0
    total = x*y
    for i in range(0,x):
        for j in range(0,y):
            t = time.time()
            NSIDC_sample = [] 
            model_sample = []
            for arr in arraystack_1:
                model_sample.append(arr[i,j])
            for arr in arraystack_2:
                NSIDC_sample.append(arr[i,j])
            NSIDC_sample = np.asarray(NSIDC_sample)
            model_sample = np.asarray(model_sample)
            print "doing permutation test {}/{}".format(count,total)
            pvals[i,j] = permutation_test(NSIDC_sample,model_sample,
                                            method='approximate',
                                            num_rounds=1000,
                                            seed=0)
            elapsed = time.time() - t
            print pvals[i,j]
            print "time required to do one operation is {}".format(elapsed)
            count += 1
    np.save('/home/ben/Desktop/pvals.npy',pvals)
    print "permutation_test done"

if __name__=="__main__":
    total_ice_diff("/media/windowsshare/","u-au866",2) 
