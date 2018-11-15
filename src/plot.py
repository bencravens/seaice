"""script file to plot grabbed variables over time"""
import grab
import mystats
from netCDF4 import *
import os
import numpy as np
import matplotlib.pyplot as plt

def main():
	"""grabbing dataset for all models, calculating seasonal ice area, plotting."""
	au866 = grab.get("/media/windowsshare","u-au866")	
	ice_area_au866 = mystats.seasonal_ice_area(au866)
	plot(range(1,13),ice_area_au866,'Months of the year', 'Ice area (m^2)', 'Seasonal')

def plot(input_x,input_y,xlab,ylab,title):
	"""plots a given dataset given inputs, titles and graph labels etc"""
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
	plt.show()
	
if __name__=='__main__':
	main()
