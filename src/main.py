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

if __name__=="__main__":
    lons,lats,aice = grab.NSIDC_data("/media/windowsshare/NSIDC_ben/ice","2")
    print np.shape(lons), np.shape(lats), np.shape(aice)
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, aice, latlon=True,cmap="jet")
    plt.show()
    
