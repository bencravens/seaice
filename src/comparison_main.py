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

def permutation_test(modelname,monthnum):
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
        pvals = np.load('/home/ben/Desktop/pvals_{}_{}.npy'.format(modelname,monthnum))
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

        #building "history" of each grid cell
        [x,y] = np.shape(aice_model_stack_regrid[0])
        pvals = np.zeros([x,y])
        count = 0
        total = x*y
        for i in range(0,x):
            for j in range(0,y):
                t = time.time()
                NSIDC_sample = [] 
                model_sample = []
                for arr in aice_model_stack_regrid:
                    model_sample.append(arr[i,j])
                for arr in aice_stack_regrid:
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
        np.save('/home/ben/Desktop/pvals_{}_{}.npy'.format(modelname,monthnum),pvals)
        
    #loaded arrays. Lets now plot them, masking the areas where there is no significance
    print "testing"
    print "making grid to plot on"
    aice_model, aice_NSIDC, anom, glons, glats = process.regrid(aice_model,lats_model,lons_model,aice,lats,lons,modelname, monthnum)
    #mask the values where the pvalue is >0.05 (insignificant) or == 0 (masked entry)
    anom = np.ma.masked_where(pvals>0.05,anom)
    m = Basemap(resolution='h', projection='spstere',lat_0=-90, lon_0=-180, boundinglat=-55)
    
    #plot this
    print "anom, glons, glats"
    print np.shape(anom), np.shape(glons), np.shape(glats)
    cs=m.pcolormesh(glons,glats,anom,latlon=True,cmap='seismic')
    plt.clim(-1.0,1.0)
    m.drawcoastlines()
    cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
    plt.title("Mean ice concentration anomaly\n Between model scheme {} and NSIDC data in the month of {}".format(modelname,monthnum))
    plt.show()
    plt.close()

if __name__=="__main__":
    permutation_test("u-at053",2)
        
        

        

        
                
