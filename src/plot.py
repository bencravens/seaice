"""script file to plot grabbed variables over time"""
import grab
import grabtest
import mystats
from netCDF4 import *
import os
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
	"""grabbing dataset for all models, calculating seasonal ice area, plotting."""
	#timing to see how much numpy speeds up
	t = time.time()

	#vectorizing plotting routine
	models = ["at053", "au866", "av231", "au872", "au874"]
	testmodels = ["at053"]
	for i, model in enumerate(testmodels):
		tempstd,tempmean = grab.ice_area_seasonal("/media/windowsshare", "u-{}".format(model))
		print "the size of tempmean is {}".format(tempmean)
		print type(tempmean)
		plot(range(1,13),tempmean,'Months of the year','Sea ice area (m^2)', '{} - Seasonal'.format(model),[3,2,i+1],tempstd)
		del tempstd
		del tempmean
	elapsed = time.time() - t
	print "Elapsed time is {}".format(elapsed)
	plt.show()

def plot(input_x,input_y,xlab,ylab,title,plotarr,std_devs):
	"""plots a given dataset given inputs, titles and graph labels etc"""
	plt.subplot(plotarr[0],plotarr[1], plotarr[2])
	plt.plot(input_x,input_y)
	plt.ylabel(ylab)
	plt.xlabel(xlab)
	plt.title(title)
	xmin = min(input_x)
	xmax = max(input_x)
	ymin = min(input_y)
	ymax = max(input_y)
	plt.xlim([xmin,xmax])
	plt.ylim([ymin*0.9,ymax*1.1])
	
	#now comes more stuff... adding max/min lines
	plt.axhline(y=ymax,color='r',linestyle="-")
	plt.axhline(y=ymin,color='r',linestyle="-")

	#now plotting the fill between the mean value +- the standard deviation
	plt.fill_between(input_x,input_y-std_devs,input_y+std_devs,facecolor='green',alpha=0.2)

if __name__=='__main__':
	main()
