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
    #models = ["u-at053","u-au866","u-au872","u-au874","u-av231"]
    models = ["u-at053"]
    for model in models:
        lons,lats,aice = grab.NSIDC_data("/media/windowsshare/NSIDC_ben/ice","2")
        print np.shape(lons), np.shape(lats), np.shape(aice)
        #grabbing model data
        lons_model,lats_model,aice_model = grab.ice_area_map_mean("/media/windowsshare",model,2)  
        inds = (lons_model>180.)
        lons_new = lons_model
        lons_new[inds] = lons_model[inds] - 360.0

        #chucking into plot function 
        plot.regrid(aice_model,lats_model,lons_new,aice,lats,lons,model,"February")    

    #plot grabbed things to test
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, aice, latlon=True,cmap="jet")
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    plt.show()
    """
