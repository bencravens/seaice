"""script to import variables from NetCDF files on high capacity storage."""
import os
from netCDF4 import *
import numpy as np

def ice_area_seasonal(path,modelname):
	"""imports variables from NetCDF files with specified path and variable name"""
	os.chdir("../../../../")
	os.chdir("{}/{}/{}".format(path,modelname,"ice"))
	filecount=0
	#simply counting number of files in directory.
	for filename in os.listdir('./'):
		filecount +=1
	#now we will sort the files based on month.
	monthareas= np.zeros(12)
	for filenum,filename in enumerate(os.listdir('./')):
		print filename	
		testdata = Dataset(filename)
		#grabbing month value from filename... format does not vary
		monthstr = filename[19:21]
		month = int(monthstr)-1
		print "Grabbing {}, file {} of {}".format(filename,filenum,filecount)
		aice = np.ma.squeeze(np.ma.array(testdata.variables['aice'][:,:],dtype='float64'))
		tarea = np.ma.array(testdata.variables['tarea'][:,:],dtype='float64')			
		#adding total ice area into its respective month bin
		monthareas[month] += np.ma.sum(aice*tarea)
		print "area added is {}".format(np.ma.sum(aice*tarea))
		testdata.close()
	#now calculating the mean for each month and returning that value
	for month in monthareas:
		month /= filecount/12 #we know there is an even amount for each month as we have model runs in whole years
	return monthareas
