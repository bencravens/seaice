"""script to import variables from NetCDF files on high capacity storage."""
import os
from netCDF4 import *

def get(path,modelname):
	"""imports variables from NetCDF files with specified path and variable name"""
	os.chdir("../../../../")
	os.chdir("{}/{}/{}".format(path,modelname,"ice"))
	alldata = []
	for filename in os.listdir('./'):
		testdata = Dataset(filename)
		alldata.append(testdata)
		print "storing {}".format(filename)
	return alldata
