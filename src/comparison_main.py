import grab
import process
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy
from scipy import stats
import os
import sys
import plot
import netCDF4
from scipy.interpolate import griddata
from mlxtend.evaluate import permutation_test
import time 

def permutation_test_NSIDC_plot(modelname,monthnum):
    """This function takes the data from a model run (typically 50 years in model time) and real world
    data from the NSIDC website. It then runs a permutation test on these two datasets for each individual gridcell, and returns a p-value for each gridcell. The anomaly of the mean of these two datasets is then calculated, and it is plotted, being masked if the pvalue is >0.05"""
    
    #JUST CALCULATING ANOMALY
    lons,lats,aice_stack,aice = grab.NSIDC_data("/media/windowsshare/NSIDC_ben/ice",str(monthnum))
    print np.shape(lons), np.shape(lats), np.shape(aice)
    
    #grabbing model data
    lons_model,lats_model,aice_model = grab.ice_area_map_mean("/media/windowsshare",modelname,monthnum)  
    #need to reproject longitude
    inds = (lons_model>180.)
    lons_new = lons_model
    lons_new[inds] = lons_model[inds] - 360.0 

    #now we want to do the permutation test on each gridcell... 
    #need to grab stack of ice_area data from model
    lons_stack,lats_stack, aice_model_stack = grab.month_map_data("/media/windowsshare/{}/ice".format(modelname), modelname, monthnum, "aice") 
    
    #we might already have the arrays for this particular instance.. if so it makes sense to load them
    try:
        aice_model_stack_regrid = np.load('/media/windowsshare/NSIDC_ben/aice_model_{}_{}.npy'.format(modelname,monthnum))
        aice_stack_regrid = np.load('/media/windowsshare/NSIDC_ben/aice_NSIDC_{}.npy'.format(monthnum))
        pvals = np.load('/home/ben/Desktop/pvals/pvals_{}_{}.npy'.format(modelname,monthnum))
        print "arrays loaded successfully"
    #if not, we need to generate them 
    except:
        print np.shape(aice_model_stack), np.shape(aice_stack)
        print "regridded arrays are not already stored.. need to create them for the first time"
        aice_model_stack_regrid = []
        aice_stack_regrid = []
        #now need to regrid all 
        for array in aice_model_stack:
            aice_model_temp, NSIDC_placeholder, anom, glons, glats = process.regrid(array,lats_model,lons_model,aice_stack[0],lats,lons,modelname,monthnum)
            aice_stack_regrid.append(aice_model_temp)
        for array in aice_stack:
            aice_model_placeholder, NSIDC_temp, anom, glons, glats = process.regrid(aice_model_stack[0],lats_model,lons_model,array,lats,lons,modelname,monthnum)
            aice_model_stack_regrid.append(NSIDC_temp)
        np.save('/media/windowsshare/NSIDC_ben/aice_model_{}_{}.npy'.format(modelname,monthnum),aice_model_stack_regrid)
        np.save('/media/windowsshare/NSIDC_ben/aice_NSIDC_{}.npy'.format(monthnum),aice_stack_regrid)
        # we have calculated the regridded arrays, so we will now store them for convinience. 
        # now we need to actually do the permutation test
        pvals = gridcell_history(aice_stack_regrid,aice_model_stack_regrid)
        np.save('/home/ben/Desktop/pvals/pvals_{}_{}.npy'.format(modelname,monthnum),pvals)

    #loaded arrays. Lets now plot them, masking the areas where there is no significance
    print "testing"
    print "making grid to plot on"
    fig,ax = plt.subplots(figsize=(8,8))  
    aice_model, aice_NSIDC, anom, glons, glats = process.regrid(aice_model,lats_model,lons_model,aice,lats,lons,modelname, monthnum)
    #mask the values where the pvalue is >0.05 (insignificant) or == 0 (masked entry)
    m = Basemap(resolution='h', projection='spstere',lat_0=-90, lon_0=-180, boundinglat=-55)
    
    #plot this
    print "anom, glons, glats"
    print np.shape(anom), np.shape(glons), np.shape(glats)
    anom = np.ma.masked_where(pvals>0.05,anom)
    cs=m.pcolormesh(glons,glats,anom,latlon=True,cmap='seismic')
    plt.clim(-1.0,1.0)
    m.drawcoastlines()
    cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
    plt.title("Mean ice concentration anomaly\n Between model scheme {} and NSIDC data in the month of {}".format(modelname,monthnum))
    fig.savefig("/home/ben/Desktop/permutation_test_{}_{}".format(modelname,monthnum))
    plt.show()
    plt.close()
    
def permutation_test_model_plot(modelname1,modelname2,monthnum):
    """This function takes the data from two model runs (typically 50 years in model time). It then runs a permutation test on these two datasets for each individual gridcell, and returns a p-value for each gridcell. The anomaly of the mean of these two datasets is then calculated, and it is plotted, being masked if the pvalue is >0.05"""
    
    #grabbing model data
    lons_model,lats_model,aice_model_1 = grab.ice_area_map_mean("/media/windowsshare",modelname1,monthnum)  
    lons_model,lats_model,aice_model_2 = grab.ice_area_map_mean("/media/windowsshare",modelname2,monthnum)  
    
    #need to reproject longitude
    #inds = (lons_model>180.)
    #lons_new = lons_model
    #lons_new[inds] = lons_model[inds] - 360.0 

    #now we want to do the permutation test on each gridcell... 
    #need to grab stack of ice_area data from model
    lons_stack,lats_stack, aice_model_1_stack = grab.month_map_data("/media/windowsshare/{}/ice".format(modelname1), modelname1, monthnum, "aice") 
    lons_stack,lats_stack, aice_model_2_stack = grab.month_map_data("/media/windowsshare/{}/ice".format(modelname2), modelname2, monthnum, "aice") 
    
    #we might already have the arrays for this particular instance.. if so it makes sense to load them
    try:
        #file could be saved as format modelname1_modelname2 or modelname2_modelname1
        try:
            pvals = np.load('/home/ben/Desktop/pvals/model_pvals_{}_{}_{}.npy'.format(modelname1,modelname2,monthnum))
            print "arrays loaded successfully"
        except:
            pvals = np.load('/home/ben/Desktop/pvals/model_pvals_{}_{}_{}.npy'.format(modelname2,modelname1,monthnum))
            print "arrays loaded successfully"
    #if not, we need to generate them 
    except:
        # now we need to actually do the permutation test
        pvals = gridcell_history(aice_model_1_stack,aice_model_2_stack)
        np.save('/home/ben/Desktop/pvals/model_pvals_{}_{}_{}.npy'.format(modelname1,modelname2,monthnum),pvals)

    #loaded arrays. Lets now plot them, masking the areas where there is no significance
    print "testing"
    print "making grid to plot on"
    fig,ax = plt.subplots(figsize=(8,8))  
    ax.set_facecolor('xkcd:grey')
    #mask the values where the pvalue is >0.05 (insignificant) or == 0 (masked entry)
    m = Basemap(resolution='h', projection='spstere',lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines()
    m.drawlsmask(land_color='grey',ocean_color='grey',lakes=True)
    m.drawmapboundary(linewidth=1)
    #plot this
    anom = np.ma.asarray(aice_model_1 - aice_model_2)
    anom = np.ma.masked_where(pvals>0.05,anom)
    cs=m.pcolormesh(lons_model,lats_model,anom,latlon=True,cmap='seismic')
    plt.clim(-1.0,1.0)
    cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
    plt.title("Mean ice concentration anomaly\n Between model schemes {} and {} in the month of {}".format(modelname1,modelname2,monthnum))
    plt.show()
    plt.close()

def gridcell_history(array_stack_1,array_stack_2):
    """this is a support function that gathers the instances of each gridcell from the n*m grid for every year and collects
    them in a vector, for both of the n*m*t array stacks given in the input. the vector is t long, where t is the third dimension of the array stacks
    (note that the n*m dimensions of the array stack must be the same but the t axes of the array stacks may not be the same.) It then does the 
    permutation test on each of these gridcell histories and returns the corresponding pvalue"""
    #building "history" of each grid cell
    [x,y] = np.shape(array_stack_1[0])
    pvals = np.zeros([x,y])
    count = 0
    total = x*y
    for i in range(0,x,):
        for j in range(0,y):
            t = time.time()
            array_stack_1_sample = [] 
            array_stack_2_sample = []
            for arr in array_stack_1:
                array_stack_1_sample.append(arr[i,j])
            for arr in array_stack_2:
                array_stack_2_sample.append(arr[i,j])
            array_stack_1_sample = np.asarray(array_stack_1_sample)
            array_stack_2_sample = np.asarray(array_stack_2_sample)
            print "doing permutation test {}/{}".format(count,total)
            pvals[i,j] = permutation_test(array_stack_1_sample,array_stack_2_sample,
                                        method='approximate',
                                        num_rounds=1000,
                                        seed=0)
            elapsed = time.time() - t
            print pvals[i,j]
            print "time required to do one operation is {}".format(elapsed)
            count += 1
    return pvals

if __name__=="__main__":
    #permutation_test_NSIDC_plot("u-at053",2)
    permutation_test_model_plot("u-au866","u-at053",2)
        
        

        

        
                
