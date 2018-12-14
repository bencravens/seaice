"""script file containing tools for plotting data grabbed from NETCDF4 files"""

# IMPORT LIBRARIES
import grab
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy
from scipy import stats
import os

######## ICE AREA SPECIFIC DATA GRABBING FUNCTIONS ###############################


def ice_area_seasonal_main():
    """grabbing dataset for all models, calculating seasonal ice area, plotting."""
    # vectorizing plotting routine
    models = ["at053", "au866", "av231", "au872", "au874"]
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_area_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum areas are {}".format(tempmin)
        fig, ax = plt.subplots(figsize=(8, 8))
        plt.plot(range(1, 13), tempmean)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice area (m^2)')
        plt.title('{} - Seasonal'.format(model))
        plt.fill_between(range(1, 13), tempmean-tempstd, tempmean +
                         tempstd, facecolor='green', alpha=0.2, linestyle="--")
        plt.tight_layout()  # making sure the plots do not overlap...
        plt.show()
        fig.savefig('/home/ben/Desktop/seasonal/{}-seasonal'.format(model))


def ice_area_tseries_main():
    """makes a nice time series plot for ice area of a given model..."""
    ice_area = grab.ice_area_tseries("/media/windowsshare", "u-at053")
    plt.plot(ice_area)
    plt.title("Sea ice area from 1990-2009 for control model")
    plt.xlabel("Time (months)")
    plt.ylabel("Total antarctic sea ice area (m^2)")
    plt.show()


def ice_area_month_main():
    """makes a time series plot of ice area for a given month"""
    ice_area = grab.ice_area_month(
        "/media/windowsshare", "u-at053", 2)  # grabbing data for feb, (i.e) month=2
    plt.plot(ice_area)
    plt.title("Sea ice area from 1990-2009 for the month of February")
    plt.xlabel("Years since 1990")
    plt.ylabel("Total antarctic sea ice area (m^2)")
    std_dev = stats.stdev(ice_area)*np.ones(len(ice_area))
    xvals = np.linspace(0, 19, 20, dtype='int')
    plt.fill_between(xvals, ice_area+std_dev, ice_area-std_dev,
                     facecolor='green', alpha=0.2, linestyle="--")
    plt.show()

############### GENERAL FUNCTIONS FOR MAP PLOTTING ##################################


def month_map_mean_main(modelname, monthnum, varname):
    """function called when map plots of mean of variable varname are wanted..."""
    lons, lats, myvar = grab.month_map_mean(
        "/media/windowsshare", modelname, monthnum, varname)  # grabbing data
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.fillcontinents(color='grey')
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True)
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    plt.title("map plot for {} mean {} in the month of {}".format(
        modelname, varname, monthdict[monthnum]))
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    cbar.set_label('{} in grid cell'.format(varname))
    fig.savefig(
        '/home/ben/Desktop/mapplots/{}-{}-{}'.format(modelname, varname, monthnum))


def month_map_anom_main(modelname, monthnum, varname):
    """function called when anomaly map plots of variable varname are wanted..."""
    lons, lats, myvar, total_diff = grab.month_map_anom(
        "/media/windowsshare", modelname, monthnum, varname)  # grabbing data
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.fillcontinents(color='grey')
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True, cmap='seismic')
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    plt.title("Control - {} in the month of {} \n variable plotted is {} \n total difference is {}".format(
        modelname, monthdict[monthnum], varname, total_diff))
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    cbar.set_label('change in variable {}'.format(varname))
    plt.clim(-1.0, 1.0)  # one may have to adjust clim depending on variable...
    fig.savefig(
        '/home/ben/Desktop/anomplots/{}-{}-{}'.format(modelname, varname, monthnum))


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

    # now performing the t-test
    # modelvar and controlvar should be the same shape
    assert controlvar.shape == modelvar.shape
    [z, x, y] = modelvar.shape

    # now making an empty array to chuck results into
    tstats = np.ma.zeros([x, y])
    pvals = np.ma.zeros([x, y])
    for i in xrange(x):
        for j in xrange(y):
            tstats[i, j], pvals[i, j] = scipy.stats.mstats.ttest_ind(
                controlvar[:, i, j], modelvar[:, i, j])

    tstats = np.ma.masked_where(pvals > pval_filter, tstats)
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawlsmask(land_color='grey', ocean_color='aqua', lakes=True)
    m.drawmapboundary(linewidth=1)
    if t_or_p:
        cm = m.pcolormesh(lons, lats, tstats, latlon=True, cmap='seismic')
        monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                     7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
        plt.title("tstatistic of {} for model {} in the month of {}\n plotted where pval <{}".format(
            varname, modelname, monthdict[monthnum], pval_filter))
        cbar = m.colorbar(cm, location='bottom', pad="5%")
        cbar.set_label('tstatistic of {}'.format(varname))
        fig.savefig(
            '/home/ben/Desktop/tstatplots/{}-{}-{}_tstat'.format(modelname, varname, monthnum))
    else:
        cm = m.pcolormesh(lons, lats, pvals, latlon=True, cmap='seismic')
        monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                     7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
        plt.title("pvalues of {} for model {} in the month of {}".format(
            varname, modelname, monthdict[monthnum]))
        cbar = m.colorbar(cm, location='bottom', pad="5%")
        cbar.set_label('pvalues of {}'.format(varname))
        fig.savefig(
            '/home/ben/Desktop/tstatplots/{}-{}-{}_pval'.format(modelname, varname, monthnum))


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


if __name__ == '__main__':
    models = ["u-au866", "u-au872", "u-au874", "u-av231"]
    months = [2]
    for model in models:
        for month in months:
            t_test_area_main(
                model, month, "aice", [-86.53, -63.73], [160, 230], "/home/ben/Documents/summer2019")
