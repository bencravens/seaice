"""script file to plot grabbed variables over time"""
import grab
import mystats
from netCDF4 import *
import os
import numpy as np
import matplotlib.pyplot as plt

def main():
	"""grabbing dataset for model au866"""
	au866 = grab.get("/media/windowsshare","u-au866")	
	seasonal_ice_area(au866)

def seasonal_ice_area(model_data):
	"""calculates the mean seasonal ice area for a model. Returns a 1x12 vector
	the total ice area is calculate by aice*tarea, where aice is the proportion that ice fills up the gridbox (goes from 0-1)
	and t-area is the area of each gridbox."""
	
	month_mean = np.zeros(12) # a vector that will be populated with the mean ice area per month...	
	
	#Going through each entry in each month. Finding the mean total area for each month.
	for i, month in enumerate(model_data):	
		run_num = len(month) #number of runs per month
		total_area = 0
		for run in month:
			tarea = np.ma.array(run.variables['tarea'][:,:],dtype='float64')
			aice = np.ma.squeeze(np.ma.array(run.variables['aice'][:,:],dtype='float64'))
			total_area += np.ma.sum(aice*tarea)
		month_mean[i] = (total_area/run_num)
	
	#now plotting 
	plt.plot(month_mean)
	plt.show()
	plt.ylabel('Ice area (m^2)')
	plt.xlabel('Months of the year')
	plt.title('Seasonal')

if __name__=='__main__':
	main()
