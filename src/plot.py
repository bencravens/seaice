"""script file to plot grabbed variables over time"""
import grab
import grabtest
import mystats
from netCDF4 import *
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import statistics as stats
from mpl_toolkits.basemap import Basemap
import cmocean.cm as cmo

def ice_area_seasonal_main():
	"""grabbing dataset for all models, calculating seasonal ice area, plotting."""
	#vectorizing plotting routine
	models = ["at053", "au866", "av231", "au872", "au874"]
	models2 = ["au866","au874","av231","at053","au872"]
	testmodels = ["au874"]
	for i, model in enumerate(models):
		tempstd,tempmean,tempmax,tempmin = grab.ice_area_seasonal("/media/windowsshare", "u-{}".format(model))
		print "the minimum areas are {}".format(tempmin)
		plot(range(1,13),tempmean,'Months of the year','Sea ice area (m^2)', '{} - Seasonal'.format(model),[3,2,i+1],tempstd,tempmax,tempmin)
	elapsed = time.time() - t
	print "Elapsed time is {}".format(elapsed)
	plt.show()

def ice_area_tseries_main():
	"""makes a nice time series plot for ice area of a given model..."""
	ice_area = grab.ice_area_tseries("/media/windowsshare","u-at053")
	plt.plot(ice_area)
	plt.title("Sea ice area from 1990-2009 for control model")
	plt.xlabel("Time (months)")
	plt.ylabel("Total antarctic sea ice area (m^2)")
	plt.show()

def ice_area_month_main():
	"""makes a time series plot of ice area for a given month"""
	ice_area = grab.ice_area_month("/media/windowsshare","u-at053",2) #grabbing data for feb, (i.e) month=2	
	plt.plot(ice_area)
	plt.title("Sea ice area from 1990-2009 for the month of February")
	plt.xlabel("Years since 1990")
	plt.ylabel("Total antarctic sea ice area (m^2)")
	std_dev = stats.stdev(ice_area)*np.ones(len(ice_area))
	xvals = np.linspace(0,19,20,dtype='int')
	plt.fill_between(xvals,ice_area+std_dev,ice_area-std_dev,facecolor='green',alpha=0.2,linestyle="--")
	plt.show()

def ice_area_month_map_main(modelname,monthnum):
	"""function called when map plots are wanted..."""
	lons,lats,aice = grab.ice_area_month_map("/media/windowsshare",modelname,monthnum) #grabbing data for feb
	fig,ax=plt.subplots(figsize=(8,8))
	m = Basemap(resolution='h',projection='spstere',lat_0=-90,lon_0=-180,boundinglat=-55)
	m.drawcoastlines(linewidth=1)
	m.fillcontinents(color='grey')
	m.drawmapboundary(linewidth=1)
	cm = m.pcolormesh(lons,lats,aice,latlon=True)
	monthdict={1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr', 5:'May', 6:'Jun', 7:'Jul', 8:'Aug', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}	
	plt.title("map plot for {} mean ice area in the month of {}".format(modelname,monthdict[monthnum])) 
	cbar = m.colorbar(cm,location='bottom',pad="5%")
	cbar.set_label('% ice area in grid cell (0.0-1.0)')
	fig.savefig('/home/ben/Desktop/mapplots/{}-{}'.format(modelname,monthnum))	

def ice_area_month_map_anom_main(modelname,monthnum):
	"""function called when anomaly map plots are wanted..."""
	lons,lats,aice = grab.ice_area_month_map_anom("/media/windowsshare",modelname,monthnum) #grabbing data for feb
	fig,ax=plt.subplots(figsize=(8,8))
	m = Basemap(resolution='h',projection='spstere',lat_0=-90,lon_0=-180,boundinglat=-55)
	m.drawcoastlines(linewidth=1)
	m.fillcontinents(color='grey')
	m.drawmapboundary(linewidth=1)
	cm = m.pcolormesh(lons,lats,aice,latlon=True,cmap='seismic')
	monthdict={1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr', 5:'May', 6:'Jun', 7:'Jul', 8:'Aug', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}	
	plt.title("map plot for mean difference in ice area in the month of {}\n between control model u-at053 and {}".format(monthdict[monthnum],modelname)) 
	cbar = m.colorbar(cm,location='bottom',pad="5%",cmap=cmo.ice)
	cbar.set_label('% ice area in grid cell (0.0-1.0)')
	plt.clim(-1.0,1.0) #we want an absolute scale for the colorbar
	fig.savefig('/home/ben/Desktop/anomplots/{}-{}'.format(modelname,monthnum))	


def plot(input_x,input_y,xlab,ylab,title,plotarr,std_devs,maxes,mins):
	"""plots a given dataset given inputs, titles and graph labels etc"""
	plt.subplot(plotarr[0],plotarr[1], plotarr[2])
	plt.plot(input_x,input_y)
	plt.plot(input_x,maxes,color='r',linestyle=":")
	plt.plot(input_x,mins,color='r',linestyle=":")
	plt.ylabel(ylab)
	plt.xlabel(xlab)
	plt.title(title)
	xmin = min(input_x)
	xmax = max(input_x)
	ymin = min(input_y)
	ymax = max(input_y)
	plt.xlim([xmin,xmax])
	plt.ylim([ymin*0.9,ymax*1.1])
	#now plotting the fill between the mean value +- the standard deviation
	plt.fill_between(input_x,input_y-std_devs,input_y+std_devs,facecolor='green',alpha=0.2,linestyle="--")
	plt.tight_layout() #making sure the plots do not overlap...

if __name__=='__main__':
	models = ["au866", "av231", "au872", "au874"]
	testmodels = ["at053"]
	for model in models:
		ice_area_month_map_anom_main("u-{}".format(model),2)
		ice_area_month_map_anom_main("u-{}".format(model),9)
