"""script file to plot grabbed variables over time"""
import grab
import grabtest
import mystats
from netCDF4 import *
import os
import numpy as np
import matplotlib.pyplot as plt

def main():
	"""grabbing dataset for all models, calculating seasonal ice area, plotting."""

	#vectorizing plotting routine
	models = ["at053", "au866", "av231", "au872", "au874"]
	ice_areas = []
	for i, model in enumerate(models):
		tempmodel = grab.ice_area_seasonal("/media/windowsshare", "u-{}".format(model))
		ice_areas.append(tempmodel)
		plot(range(1,13),ice_areas[i],'Months of the year','Sea ice area (m^2)', '{} - Seasonal'.format(model),[3,2,i+1])
	plt.show()

def plot(input_x,input_y,xlab,ylab,title,plotarr):
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
	plt.ylim([ymin,ymax])
	
if __name__=='__main__':
	main()
