"""script to import variables from netCDF files on hcs, and then calculate the total ice area"""
import os
from netCDF4 import *
import numpy as np

def ice_area_seasonal(path,modelname):
	os.chdir('../../../../')
	os.chdir("{}/{}/{}".format(path,modelname,"ice"))
	print os.listdir('./')
	
	#now trying to open each file
	for filename in os.listdir('./'):
		testdata = Dataset(filename)
		print "successful"
