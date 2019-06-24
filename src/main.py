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

if __name__=="__main__":
    #grabbing NSIDC data
    lons,lats,aice = grab.NSIDC_data("/media/windowsshare/NSIDC_ben/ice","2")
    print np.shape(lons), np.shape(lats), np.shape(aice)
    #grabbing model data
    lons_model,lats_model,aice_model = grab.ice_area_map_mean("/media/windowsshare","u-at053",2)  
    #chucking into plot function 
    plot.regrid(aice_model,lats_model,lons_model,aice,lats,lons)    
