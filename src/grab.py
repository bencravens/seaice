"""script to import variables from NetCDF files on high capacity storage."""
import os
from netCDF4 import *

def get(path,modelname):
	"""imports variables from NetCDF files with specified path and variable name"""
	os.chdir("../../../../")
	os.chdir("{}/{}/{}".format(path,modelname,"ice"))
	filecount=0
	#simply counting number of files in directory.
	for filename in os.listdir('./'):
		filecount +=1
	#now we will sort the files based on month.
	alldata=[[],[],[],[],[],[],[],[],[],[],[],[]]
	for filename in os.listdir('./'):
		testdata = Dataset(filename) #grabbing dataset from file
		monthstr = filename[19:21] #grabbing month value from filename... format does not vary
		month = int(monthstr)-1 #chucking dataset into alldata's "month" boxes.. alldata[0] is jan, alldata [1] is feb ect.
		print monthstr
		alldata[month].append(testdata)
	return alldata
