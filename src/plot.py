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

def ice_area_seasonal_main():
	"""grabbing dataset for all models, calculating seasonal ice area, plotting."""
	#timing to see how much numpy speeds up
	t = time.time()

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
	ice_area = grab.ice_area_tseries("/media/windowsshare","u-at053")
	ice_area_std = stats.stdev(ice_area)
	plt.plot(ice_area)
	plt.title("Sea ice area from 1990-2009 for control model")
	plt.xlabel("Time (months)")
	plt.ylabel("Total antarctic sea ice area (m^2)")
	plt.show()	

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
#	ice_area_tseries_main()
	ice_area_seasonal_main()
