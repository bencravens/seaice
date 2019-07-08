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

if __name__=="__main__":
        #grabbing NSIDC data
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
        try:
            aice_model_stack_regrid = np.load('/media/windowsshare/NSIDC_ben/aice_model_stack_regrid.npy')
            aice_stack_regrid = np.load('/media/windowsshare/NSIDC_ben/aice_stack_regrid.npy')
            pvals = np.load('/home/ben/Desktop/pvals.npy')
            print "arrays loaded successfully"
        except:
            print np.shape(aice_model_stack), np.shape(aice_stack)
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
            np.save('/home/ben/Desktop/pvals.npy',pvals)
        
        #loaded arrays. Lets now plot them, masking the areas where there is no significance
        print "testing"
        print "making grid to plot on"
        m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
        [nx,ny] = np.shape(aice_model)
        glons,glats = m.makegrid(nx,ny)
        
        #now we just want to take the mean anomaly
        aice_model_regrid = np.mean(aice_model_stack_regrid,axis=0)
        aice_regrid = np.mean(aice_stack_regrid,axis=0)
        anom = aice_regrid - aice_model_regrid
        
        #mask the values where the pvalue is >0.05 (insignificant) or == 0 (masked entry)
        anom = np.ma.masked_where(pvals>0.05,anom)
        anom = np.ma.masked_where(pvals==0,anom)
        
        #plot this
        cs=m.pcolormesh(glons,glats,anom,latlon=True,cmap='seismic')
        plt.clim(-1.0,1.0)
        m.drawcoastlines()
        cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
        plt.title("Mean ice concentration anomaly\n Between model scheme {} and NSIDC data in the month of {}".format("u-at053","2"))
        plt.show()
        fig.savefig("/home/ben/Desktop/anom_perm")
        plt.close()
        
        

        

        
                
