"""script file to plot grabbed variables over time"""
import grab
import mystats
from netCDF4 import *
import os

def main():
	"""grabbing dataset for model au866"""
	au866 = grab.get("/media/windowsshare","u-au866")	
	seasonal_ice_area(au866)

def seasonal_ice_area(model_data):
	"""calculates the mean seasonal ice area for a model. Returns a 1x12 vector
	the total ice area is calculate by aice*tarea, where aice is the proportion that ice fills up the gridbox (goes from 0-1)
	and t-area is the area of each gridbox."""
	for run in model_data:
		aice = run.variables['aice']
		tarea = run.variables['tarea']
		print aice, tarea
	
if __name__=='__main__':
	main()
