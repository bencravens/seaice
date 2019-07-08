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

if __name__=="__main__":
    #grabbing NSIDC data
    #models = ["u-at053","u-au866","u-au872","u-au874","u-av231"]
    models = ["u-at053"]
    for model in models:
        #JUST CALCULATING ANOMALY
        lons,lats,aice_stack,aice = grab.NSIDC_data("/media/windowsshare/NSIDC_ben/ice","2")
        print np.shape(lons), np.shape(lats), np.shape(aice)
        #grabbing model data
        lons_model,lats_model,aice_model = grab.ice_area_map_mean("/media/windowsshare",model,2)  
        inds = (lons_model>180.)
        lons_new = lons_model
        lons_new[inds] = lons_model[inds] - 360.0
        
        #PRINTING TOTAL NUMBER OF POINTS THAT ARE NOT MASKED IN EITHER ARRAY
        print np.sum(aice>0)
        
        #chucking into plot function 
        #plot.regrid(aice_model,lats_model,lons_new,aice,lats,lons,model,"February")    
        
        #now we want to do the permutation test on each gridcell... 
        #need to grab stack of ice_area data from model
        lons_stack,lats_stack, aice_model_stack = grab.month_map_data("/media/windowsshare/u-at053/ice", model, 2, "aice")
        print np.shape(aice_model_stack), np.shape(aice_stack)
        
        try:
            aice_model_stack_regrid = np.load('/media/windowsshare/NSIDC_ben/aice_model_stack_regrid.npy')
            aice_stack_regrid = np.load('/media/windowsshare/NSIDC_ben/aice_stack_regrid.npy')
            print "arrays loaded successfully"
        except:
            print "regridded arrays are not already stored.. need to create them for the first time"
            aice_model_stack_regrid = []
            aice_stack_regrid = []
            #now need to regrid all 
        
            for array in aice_model_stack:
                aice_model_temp, NSIDC_placeholder, anom, glons, glats = process.regrid(array,lats_model,lons_model,aice_stack[0],lats,lons,"u-at053","2")
                aice_stack_regrid.append(aice_model_temp)
            for array in aice_stack:
                aice_model_placeholder, NSIDC_temp, anom, glons, glats = process.regrid(aice_model_stack[0],lats_model,lons_model,array,lats,lons,"u-at053","2")
                aice_model_stack_regrid.append(NSIDC_temp)
            np.save('/media/windowsshare/NSIDC_ben/aice_model_stack_regrid.npy',aice_model_stack_regrid)
            np.save('/media/windowsshare/NSIDC_ben/aice_stack_regrid.npy',aice_stack_regrid)
        
        #first testing that the arrays I've grabbed are actually fine
        fig,ax=plt.subplots(figsize=(8,8))
        m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
        m.drawmapboundary(fill_color='grey')
        cs=m.pcolormesh(lons,lats,aice_stack[i],latlon=True)
        plt.clim(0,1.0)
        m.drawcoastlines()
        cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
        plt.title("TESTING AICE STACK {}/{}".format(i,len(aice_stack)))
        plt.show()
        plt.close("all")
    
        #first testing that the arrays I've grabbed are actually fine
        fig,ax=plt.subplots(figsize=(8,8))
        m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
        m.drawmapboundary(fill_color='grey')
        cs=m.pcolormesh(lons_stack,lats_stack,aice_model_stack[0],latlon=True)
        plt.clim(0,1.0)
        m.drawcoastlines()
        cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
        plt.title("TESTING AICE MODEL STACK")
        plt.show()
        plt.close("all")
        
        #building "history" of each grid cell
        [x,y] = np.shape(aice_model_stack_regrid[0])
        pvals = []
