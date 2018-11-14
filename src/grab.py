"""script to import variables from NetCDF files on high capacity storage."""
import os

def get(path, var, modelname,date):
	"""imports variables from NetCDF files with specified path and variable name"""
	os.chdir("../../../../")
	os.chdir("{}/{}".format(path,modelname))
	print os.listdir("./")	

