"""script to import variables from NetCDF files on high capacity storage."""
import os
from netCDF4 import *
import numpy as np
import statistics as stats
import pickle

def ice_area_seasonal(path,modelname):
	"""imports variables from NetCDF files with specified path and variable name"""
	os.chdir("../../../../")
	os.chdir("{}/{}/{}".format(path,modelname,"ice"))
	filecount=len(os.listdir('./'))
	
	#now we will sort the files based on month.
	monthareas=np.ma.zeros([12,filecount/12],dtype='float64') # assumes same number of files for each month
	monthcount = np.zeros(12)	
	stdevs = np.ma.zeros(12,dtype='float64') 
	means = np.ma.zeros(12,dtype='float64')
	for filenum,filename in enumerate(os.listdir('./')):
		print filename	
		testdata = Dataset(filename)
		#now grabbing latitude just to check if we are in the southern hemisphere (some data is full world)
		lats = np.ma.array(testdata.variables['TLAT'][:,:],dtype='float64')
		cond = lats <= 0 #if we are below the equator, grab
		if filenum==0:
			tarea = np.ma.array(testdata.variables['tarea'][:,:],dtype='float64')[cond]
		#grabbing month value from filename... format does not vary
		monthstr = filename[19:21]
		monthnum= int(monthstr)-1
	
		print "the month we are grabbing is {}".format(monthstr)
		print "Grabbing {}, file {} of {}".format(filename,filenum,filecount)
		aice = np.ma.squeeze(np.ma.array(testdata.variables['aice'][:,:],dtype='float64'))[cond] 
		#adding total ice area into its respective month bin
		monthareas[monthnum,monthcount[monthnum]] = np.ma.sum(aice*tarea)
		monthcount[monthnum] =+1	
		print "area added is {}".format(np.ma.sum(aice*tarea))
		testdata.close()
		#now converting to numpy array
	
	#now calculating the mean for each month and returning that value, as well as the standard deviation for each month.
	for i in xrange(12):
		stdevs[i] = np.ma.std(monthareas[i,:])
		means[i] = np.ma.mean(monthareas[i,:])


	return stdevs, means

def ice_area_tseries(path,modelname):
	"""imports variables from NetCDF files with specified path and variable name"""
	os.chdir("../../../../")
	os.chdir("{}/{}/{}".format(path,modelname,"ice"))
	filecount=len(os.listdir('./'))
	
	#declaring array to hold ice areas for each timestep
	ice_area = []
	
	#sorting the files and appending the mean of each to an array	
	for i,filename in enumerate(sorted(os.listdir('./'))):
		print filename #just making sure we are grabbing files in order...
		testdata = Dataset(filename)
		if i==0:
			tarea = np.ma.array(testdata.variables['tarea'][:,:],dtype='float64')
		aice = np.ma.squeeze(np.ma.array(testdata.variables['aice'][:,:],dtype='float64')) 
		ice_area.append(np.ma.sum(aice*tarea))
		testdata.close()
	
	#Now that we have all of the data we will return it
	return ice_area
	
