"""script to perform statistical operations on netCDF file variable data"""
import numpy as np

def monthlymean(var, startyear, endyear):
	"""calculates the monthly mean for a specified netCDF variable"""
	pass

def seasonal_ice_area(model_data):
	"""calculates the mean seasonal ice area for a model. Returns a 1x12 vector
	the total ice area is calculate by aice*tarea, where aice is the proportion that ice fills up the gridbox (goes from 0-1)
	and t-area is the area of each gridbox."""
	
	month_mean = np.zeros(12) # a vector that will be populated with the mean ice area per month...	
	
	#Going through each entry in each month. Finding the mean total area for each month.
	for i, month in enumerate(model_data):	
		run_num = len(month) #number of runs per month
		total_area = 0
		if (i==0):
			#simply loading the tarea. It is the same for all models so we simply load it from the first file.
			tempdata = model_data[0][0]
			tarea = np.ma.array(tempdata.variables['tarea'][:,:],dtype='float64')
		for run in month:
			aice = np.ma.squeeze(np.ma.array(run.variables['aice'][:,:],dtype='float64'))
			total_area += np.ma.sum(aice*tarea)
			run.close()
		month_mean[i] = (total_area/run_num)
	
	#returning calculated means
	return month_mean


