"""script file containing tools for plotting data grabbed from NETCDF4 files"""

# IMPORT LIBRARIES
import grab
import process
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy
from scipy import stats
import os
import sys

######## SPECIFIC DATA GRABBING FUNCTIONS ###############################


def ice_area_seasonal_main():
    """grabbing dataset for all models, calculating seasonal ice area, plotting."""
    models = ["at053", "au866", "av231", "au872", "au874"]
    fig, ax = plt.subplots(figsize=(8, 8))
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_area_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum areas are {}".format(tempmin)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice area (m^2)')
        #plt.fill_between(range(1, 13), tempmean-tempstd, tempmean +
        #                 tempstd, facecolor='green', alpha=0.2, linestyle="--")
        plt.plot(range(1, 13), tempmean,label=model)
    plt.legend()
    plt.title('Seasonal map plot of mean total ice area at each month')
    plt.show()    
    fig.savefig('/home/ben/Desktop/seasonal/allmodel-seasonal-nostd')

def ice_volume_seasonal_main():
    """grabbing dataset for all models, calculating seasonal ice volume, plotting."""
    models = ["at053", "au866", "av231", "au872", "au874"]
    fig, ax = plt.subplots(figsize=(8, 8))
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_volume_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum volumes  are {}".format(tempmin)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice volume (m^2)')
        plt.plot(range(1, 13), tempmean,label=model)
        plt.fill_between(range(1, 13), tempmean-tempstd, tempmean +
                         tempstd, facecolor='green', alpha=0.2, linestyle="--")
        plt.title('Seasonal map plot of mean total ice volume for model {}'.format(model))
        plt.show()    
        fig.savefig('/home/ben/Desktop/seasonal/seasonal-volume-{}'.format(model))

def ice_area_tseries_main():
    """makes a nice time series plot for ice area of a given model..."""
    ice_area = grab.ice_area_tseries("/media/windowsshare", "u-at053")
    plt.plot(ice_area)
    plt.title("Sea ice area from 1990-2009 for control model")
    plt.xlabel("Time (months)")
    plt.ylabel("Total antarctic sea ice area (m^2)")
    plt.show()


def ice_area_month_main(modelname,monthnum):
    """makes a time series plot of ice area for a given month"""
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
             7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    ice_area = grab.ice_area_month(
        "/media/windowsshare",modelname,monthnum)
    print ice_area
    plt.title("Sea ice area runs of control model {}\nIn the month of {}, year 2000 forcing".format(
        modelname,monthdict[monthnum]))
    plt.xlabel("Run number")
    plt.ylabel("Total antarctic sea ice area (m^2)")
    stddev = np.std(ice_area)
    plt.errorbar(range(1,21),ice_area,yerr=stddev,linestyle="None", marker="o",color='b')
    plt.show()

def ice_volume_month_main(modelname,monthnum):
    """makes a time series plot of ice volume for a given month"""
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
             7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    ice_volume = grab.ice_volume_month(
        "/media/windowsshare",modelname,monthnum)
    print ice_volume
    plt.title("Sea ice volume runs of control model {}\nIn the month of {}, year 2000 forcing".format(
        modelname,monthdict[monthnum]))
    plt.xlabel("Run number")
    plt.ylabel("Total antarctic sea ice volume (m^2)")
    stddev = np.std(ice_volume)
    plt.errorbar(range(1,21),ice_volume,yerr=stddev,linestyle="None", marker="o",color='b')
    plt.show()


############### GENERAL FUNCTIONS FOR MAP PLOTTING ##################################


def month_map_mean_main(modelname,monthnum,varname,csvdir,isice):
    """function called when map plots of mean of variable varname are wanted..."""
    lons, lats, myvar,units = grab.month_map_mean(
        "/media/windowsshare",modelname,monthnum,varname,isice)  # grabbing data
    #now saving limits of plot to csv file so that month_map_anom_main can use them
    process.save_lims(modelname,monthnum,varname,csvdir,myvar)
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    #m.drawlsmask(land_color='grey', ocean_color='black', lakes=True)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True)
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    if units=="1":
        plt.title("map plot for {} mean {} in the month of {}.\n With units of fractional area (0.0-1.0)".format(
            modelname, varname, monthdict[monthnum]))
    else:
        plt.title("map plot for {} mean {} in the month of {}.\n With units {}".format(
        modelname, varname, monthdict[monthnum],units))
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    cbar.set_label('{} in grid cell'.format(varname))
    lims = process.read_lims(modelname,monthnum,varname,csvdir)
    #plt.clim(float(lims["Min"]),float(lims["Max"]))
    plt.show()
    fig.savefig(
        '/home/ben/Desktop/mapplots/{}-{}-{}'.format(modelname, varname, monthnum))
    plt.close()

def month_map_anom_main(modelname, monthnum, varname,csvdir,isice):
    """function called when anomaly map plots of variable varname are wanted..."""
    lons, lats, myvar, total_diff, units = grab.month_map_anom(
        "/media/windowsshare", modelname, monthnum, varname,isice)  # grabbing data
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawmapboundary(linewidth=1)
    #m.drawlsmask(land_color='grey', ocean_color='black', lakes=True)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True, cmap='seismic')
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    if units=="1":
        plt.title("{} - Control in the month of {} \n variable plotted is {} \n total difference is {}, units are fractional area (0.0-1.0)".format(
            modelname, monthdict[monthnum], varname, total_diff))
    else:
        plt.title("{} - Control in the month of {} \n variable plotted is {} \n total difference is {}, units are {}".format(
            modelname, monthdict[monthnum], varname, total_diff,units))
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    cbar.set_label('change in variable {}'.format(varname))
    #reading in plot limit dictionary from csv file to be consistent with month_map_mean (assumes month_map_mean has been run...)
    lims = process.read_lims(modelname,monthnum,varname,csvdir)
    print "LIMS INCOMING"
    print lims
    plt.clim(float(lims["Min"]),float(lims["Max"]))
    plt.show()
    fig.savefig(
        '/home/ben/Desktop/anomplots/{}-{}-{}-anom'.format(modelname, varname, monthnum))
    plt.close()

def month_map_variance_main(modelname, monthnum, varname):
    """plots the variance of a seasonal variable myvar for a given model and month"""
    lons, lats, myvar = grab.month_map_stddev(
        "/media/windowsshare", modelname, monthnum, varname)
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.fillcontinents(color='grey')
    m.drawmapboundary(linewidth=1)
    # doing np.multiply to get the variance, which is the square of the std-dev
    cm = m.pcolormesh(lons, lats, np.multiply(
        myvar, myvar), latlon=True, cmap='Blues')
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    plt.title("variance of {} for model {} in the month of {}".format(
        varname, modelname, monthdict[monthnum]))
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    cbar.set_label('variance in variable {}'.format(varname))
    fig.savefig(
        '/home/ben/Desktop/varplots/{}-{}-{}'.format(modelname, varname, monthnum))


def t_test_main(modelname, monthnum, varname, pval_filter, t_or_p):
    """performs the student's t-test on each gridpoint. if t_or_p is true, plots the t-statistic scalar field, else plots the p value."""
    lons, lats, modelvar = grab.month_map_data(
        "/media/windowsshare", modelname, monthnum, varname)
    lons, lats, controlvar = grab.month_map_data(
        "/media/windowsshare", "u-at053", monthnum, varname)
    #formatting info for plot title based on whether we want to plot tstat or pvalue
    if t_or_p:
        varstring = "tstatistic"
        plotinfo = "plotted where pval < {}".format(pval_filter)
    else:
        varstring = "pvalue"
        plotinfo = ""
    result = process.t_test_gridpoint(lons,lats,modelvar,controlvar,pval_filter,t_or_p)
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawlsmask(land_color='grey', ocean_color='aqua', lakes=True)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, result, latlon=True, cmap='seismic')
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    plt.title("{} of {} for model {} in the month of {}\n{}".format(
        varstring,varname, modelname, monthdict[monthnum], plotinfo))
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    cbar.set_label('tstatistic of {}'.format(varname))
    fig.savefig(
        '/home/ben/Desktop/tstatplots/{}-{}-{}_{}-TEST'.format(modelname, varname, monthnum,varstring))


def t_test_area_main(modelname, monthnum, varname, latrange, lonrange, outputdir):
    """performs the student's t-test on a given rectangular feature area.
    compares with control model u-at053 by stripping points in selected area of
    spatial property and treats them as a sequence of data. plots the selected area as a 
    method of double checking, and then writes the result of each t-test to a file ttests.txt
    in directory given by outputdir"""
    lons, lats, modelvar = grab.month_map_mean(
        "/media/windowsshare", modelname, monthnum, varname)
    lons1, lats1, controlvar = grab.month_map_mean(
        "/media/windowsshare", "u-at053", monthnum, varname)
    latmin = latrange[0]
    latmax = latrange[1]
    lonmin = lonrange[0]
    lonmax = lonrange[1]

    # making a mask so that we only plot the values which fall inside the bounded range for the model
    latcond = np.logical_or(lats < latmin, lats > latmax)
    loncond = np.logical_or(lons < lonmin, lons > lonmax)
    cond = np.logical_or(latcond, loncond)
    # now doing it for the control
    latcond1 = np.logical_or(lats1 < latmin, lats1 > latmax)
    loncond1 = np.logical_or(lons1 < lonmin, lons1 > lonmax)
    cond1 = np.logical_or(latcond1, loncond1)
    # now masking the out of bounds values
    modelvar_masked = np.ma.masked_where(cond, modelvar)
    controlvar_masked = np.ma.masked_where(cond1, controlvar)

    # now we will plot these to see if they make sense
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawlsmask(land_color='grey', ocean_color='aqua', lakes=True)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, modelvar_masked, latlon=True, cmap='Reds')
    plt.title("{} in selected area for {}".format(varname, modelname))
    plt.show()
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawlsmask(land_color='grey', ocean_color='aqua', lakes=True)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons1, lats1, controlvar_masked,
                      latlon=True, cmap='Reds')
    plt.title("{} in selected area for {}".format(varname, "u-at053"))
    plt.show()

    # now that we've done that, it's time to do the area-wise t-test.
    # first stripping modelvar and controlvar of sparial data.. converting them into a sequence

    # now we want the entries which do NOT match the prior condition of being outside of the selected area.
    modelvar_seq = modelvar[np.logical_not(cond)]
    controlvar_seq = controlvar[np.logical_not(cond1)]

    # now finally performing the t-test and printing the results.
    [tstat, pval] = scipy.stats.ttest_ind(modelvar_seq, controlvar_seq)
    print "Results of area-wide t-test for {}: tstat {}, pval {}".format(
        modelname, tstat, pval)

    # now writing output to file
    file = open("{}/ttest_{}_{}.txt".format(outputdir, modelname, varname), "w")
    file.write(
        "Results of area-wide t-test for {}: tstat {}, pval {}".format(modelname, tstat, pval))
    file.close()

def scatterplot_area_main(modelname, monthnum, varname, latrange, lonrange, outputdir):
    """ visually compares modelname with control model u-at053 by stripping points in selected area of
    spatial property and treats them as a sequence of data. Then creates a scatterplot of the two arrays
    and saves figure in a given output directory outputdir"""
    lons, lats, modelvar = grab.month_map_mean(
        "/media/windowsshare", modelname, monthnum, varname)
    lons1, lats1, controlvar = grab.month_map_mean(
        "/media/windowsshare", "u-at053", monthnum, varname)
    latmin = latrange[0]
    latmax = latrange[1]
    lonmin = lonrange[0]
    lonmax = lonrange[1]

    # making a mask so that we only plot the values which fall inside the bounded range for the model
    latcond = np.logical_or(lats < latmin, lats > latmax)
    loncond = np.logical_or(lons < lonmin, lons > lonmax)
    cond = np.logical_or(latcond, loncond)
    # now doing it for the control
    latcond1 = np.logical_or(lats1 < latmin, lats1 > latmax)
    loncond1 = np.logical_or(lons1 < lonmin, lons1 > lonmax)
    cond1 = np.logical_or(latcond1, loncond1)
    # now masking the out of bounds values
    modelvar_masked = np.ma.masked_where(cond, modelvar)
    controlvar_masked = np.ma.masked_where(cond1, controlvar)

    # first stripping modelvar and controlvar of spatial data.. converting them into a sequence
    # now we want the entries which do NOT match the prior condition of being outside of the selected area.
    modelvar_seq = modelvar[np.logical_not(cond)]
    controlvar_seq = controlvar[np.logical_not(cond1)]
    
    #now making scatter plot
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.scatter(modelvar_seq,controlvar_seq)
    plt.title("scatterplot comparing variable {} for {} with control".format(varname,modelname))
    plt.xlabel(modelname)
    plt.ylabel("u-at053")
    plt.show()
    fig.savefig("{}/{}-{}-{}-scatter.png".format(outputdir,modelname,varname,monthnum))

def plot(input_x, input_y, xlab, ylab, title, plotarr, std_devs, maxes, mins):
    """plots a given dataset given inputs, titles and graph labels etc"""
    plt.subplot(plotarr[0], plotarr[1], plotarr[2])
    plt.plot(input_x, input_y)
    plt.plot(input_x, maxes, color='r', linestyle=":")
    plt.plot(input_x, mins, color='r', linestyle=":")
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(title)
    xmin = min(input_x)
    xmax = max(input_x)
    ymin = min(input_y)
    ymax = max(input_y)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin*0.9, ymax*1.1])
    # now plotting the fill between the mean value +- the standard deviation
    plt.fill_between(input_x, input_y-std_devs, input_y+std_devs,
                     facecolor='green', alpha=0.2, linestyle="--")
    plt.tight_layout()  # making sure the plots do not overlap...

if __name__=="__main__":
    month_map_mean_main("u-au866",2,"fhocn_ai","/home/ben/Documents/summer2019/plotlims",False)
    month_map_anom_main("u-au866",2,"fhocn_ai","/home/ben/Documents/summer2019/plotlims",False)
