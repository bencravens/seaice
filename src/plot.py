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
from scipy.interpolate import griddata

######## SPECIFIC DATA GRABBING FUNCTIONS ###############################

def ice_area_seasonal_main_all():
    """grabbing dataset for all models, calculating seasonal ice area, plots all on same graph."""
    # vectorizing plotting routine
    models = ["at053", "au866", "av231", "au872", "au874"]
    testmodels = ["au866"]
    fig, ax = plt.subplots(figsize=(8, 8))
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_area_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum areas are {}".format(tempmin)
        plt.plot(range(1, 13), tempmean,label=model)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice area (m^2)')
        plt.title('Seasonal mean total sea ice area'.format(model))
        plt.tight_layout()  # making sure the plots do not overlap...
        plt.ylim(0,1.8e13)
    plt.legend(bbox_to_anchor=(0.9, 0.3),bbox_transform=plt.gcf().transFigure)
    fig.savefig("/home/ben/Desktop/allmodel.png")
 
def ice_volume_seasonal_main_all():
    """grabbing dataset for all models, calculating seasonal ice volume, plots all on same graph."""
    # vectorizing plotting routine
    models = ["at053", "au866", "av231", "au872", "au874"]
    testmodels = ["au866"]
    fig, ax = plt.subplots(figsize=(8, 8))
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_volume_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum volumes are {}".format(tempmin)
        plt.plot(range(1, 13), tempmean,label=model)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice volume (m^3)')
        plt.title('Seasonal mean total sea ice volume'.format(model))
        plt.tight_layout()  # making sure the plots do not overlap...
        plt.ylim(0,2.5e13)
    plt.legend(bbox_to_anchor=(0.9, 0.3),bbox_transform=plt.gcf().transFigure)
    fig.savefig("/home/ben/Desktop/allmodel-volume.png")

def ice_volume_seasonal_main():
    """grabbing dataset for all models, calculating seasonal ice volume, plotting."""
    # vectorizing plotting routine
    models = ["at053", "au866", "av231", "au872", "au874"]
    testmodels = ["au866"]
    fig, ax = plt.subplots(figsize=(8, 8))
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_volume_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum volumes are {}".format(tempmin)
        plt.plot(range(1, 13), tempmean,label=model)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice volume (m^3)')
        plt.title('Seasonal mean total sea ice volume of {}'.format(model))
        plt.fill_between(range(1, 13), tempmean-tempstd, tempmean +
                         tempstd, facecolor='green', alpha=0.2, linestyle="--")
        plt.tight_layout()  # making sure the plots do not overlap...
        plt.ylim(0,2.5e13)
        plt.show()

def ice_area_seasonal_main():
    """grabbing dataset for all models, calculating seasonal ice area, plotting."""
    # vectorizing plotting routine
    models = ["at053", "au866", "av231", "au872", "au874"]
    testmodels = ["au866"]
    fig, ax = plt.subplots(figsize=(8, 8))
    for i, model in enumerate(models):
        tempstd, tempmean, tempmax, tempmin = grab.ice_area_seasonal(
            "/media/windowsshare", "u-{}".format(model))
        print "the minimum areas are {}".format(tempmin)
        plt.plot(range(1, 13), tempmean,label=model)
        plt.xlabel('Months of the year')
        plt.ylabel('Sea ice area (m^2)')
        plt.title('Seasonal mean total sea ice area of {}'.format(model))
        plt.fill_between(range(1, 13), tempmean-tempstd, tempmean +
                         tempstd, facecolor='green', alpha=0.2, linestyle="--")
        plt.tight_layout()  # making sure the plots do not overlap...
        plt.ylim(0,1.8e13)
        plt.show()

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
    plt.xlabel("Years in model time")
    plt.ylabel("Total antarctic sea ice area (m^2)")
    stddev = np.std(ice_area)
    plt.plot(range(1,21),ice_area,linestyle="None", marker="o",color='b')
    plt.show()

def ice_area_month_main_all(monthnum):
    """makes a time series plot of ice area for a given month for all models"""
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
             7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    symboldict = {"u-at053": "o", "u-au866": "s", "u-au872": "p", "u-au874":"v", "u-av231":"D"}
    fig, ax = plt.subplots(figsize=(8, 8))
    for model in ["u-at053", "u-au866","u-au872","u-au874","u-av231"]:
        ice_area = grab.ice_area_month(
            "/media/windowsshare",model,monthnum)
        print ice_area
        plt.title("Sea ice area runs of all models\nIn the month of {}, year 2000 forcing".format(
            monthdict[monthnum]))
        plt.xlabel("Years in model time")
        plt.ylabel("Total antarctic sea ice area (m^2)")
        plt.plot(range(1,21),ice_area,label=model,linestyle="None", marker=symboldict[model],markersize=5.0)
    plt.legend(loc='upper right', numpoints=1,prop={'size': 8})
    fig.savefig('/home/ben/Documents/summer2019docs/metoffice/{}_month_area'.format(monthdict[monthnum]))
    print "saving as /home/ben/Documents/summer2019docs/metoffice/{}_month_area".format(monthdict[monthnum])

def ice_volume_month_main(modelname,monthnum):
    """makes a time series plot of ice volume for a given month"""
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
             7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    ice_volume = grab.ice_volume_month(
        "/media/windowsshare",modelname,monthnum)
    print ice_volume
    plt.title("Sea ice volume runs of control model {}\nIn the month of {}, year 2000 forcing".format(
        modelname,monthdict[monthnum]))
    plt.xlabel("Years in model time")
    plt.ylabel("Total antarctic sea ice volume (m^3)")
    stddev = np.std(ice_volume)
    plt.plot(range(1,21),ice_volume,linestyle="None", marker="o",color='b')
    plt.show()

def ice_volume_month_main_all(monthnum):
    """makes a time series plot of ice volume for a given month for all models"""
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
             7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    symboldict = {"u-at053": "o", "u-au866": "s", "u-au872": "p", "u-au874":"v", "u-av231":"D"}
    fig,ax = plt.subplots(figsize=(8,8))
    for model in ["u-at053", "u-au866","u-au872","u-au874","u-av231"]:
        ice_volume = grab.ice_volume_month(
            "/media/windowsshare",model,monthnum)
        print ice_volume
        plt.title("Sea ice volume runs of all models\nIn the month of {}, year 2000 forcing".format(
            monthdict[monthnum]))
        plt.xlabel("Years in model time")
        plt.ylabel("Total antarctic sea ice volume (m^3)")
        plt.plot(range(1,21),ice_volume,label=model,linestyle="None", marker=symboldict[model],markersize=5.0)
    plt.legend(loc='upper right', numpoints=1,prop={'size':8})
    fig.savefig('/home/ben/Documents/summer2019docs/metoffice/{}_month_volume'.format(monthdict[monthnum]))
    print "saving as /home/ben/Documents/summer2019docs/metoffice/{}_month_volume".format(monthdict[monthnum])

############### GENERAL FUNCTIONS FOR MAP PLOTTING ##################################


def month_map_mean_main(modelname,monthnum,varname,csvdir,isice):
    """function called when map plots of mean of variable varname are wanted..."""
    #dictionary to set limits on some plots manually as outliers obscure detail of data. quick fix, will change later
    limitdict = {"ardg":[0.0,0.72],"fhocn_ai":[-80.0,0],"fsurf_ai":[-60.0,0],"siflcondtop":[-80.0,0.0],"siflsensupbot":[-2400.0,0],"sihc":[-1.6e9,0.0], "sithick":[0,6], "dardg1dt":[0,5], "opening":[0,12.5]}
    lons, lats, myvar,units = grab.month_map_mean(
        "/media/windowsshare",modelname,monthnum,varname,isice)  # grabbing data
    #now saving limits of plot to csv file so that month_map_anom_main can use them
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawlsmask(land_color='grey',ocean_color='grey',lakes=True)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True,cmap="jet")
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    if units=="1":
        plt.title("{} mean {} in the month of {}".format(
            modelname, varname, monthdict[monthnum])) 
        cbar.set_label('{}[fractional area]'.format(varname))
    else:
        plt.title("{} mean {} in the month of {}".format(
        modelname, varname, monthdict[monthnum]))
        cbar.set_label('{}[{}]'.format(varname,units))
    try:
        plt.clim(limitdict[varname][0],limitdict[varname][1])
    except:
        print "this variable does not have preset plot limits. Allowing matplotlib to set them."
    fig.savefig(
        '/home/ben/Desktop/mapplots/{}_{}_{}'.format(modelname, varname, monthnum))
    plt.close()

def month_map_anom_main(modelname, monthnum, varname,csvdir,isice):
    """function called when anomaly map plots of variable varname are wanted..."""
    lons, lats, myvar, total_diff, units = grab.month_map_anom_test(
        "/media/windowsshare", modelname, monthnum, varname,isice)  # grabbing data
    print("total difference in variable {} is {}".format(varname,total_diff))
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawmapboundary(linewidth=1)
    m.drawlsmask(land_color='grey', ocean_color='grey', lakes=True)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True, cmap='seismic')
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    plt.title("{} {} anomaly in the month of {}".format(modelname,varname,monthdict[monthnum]))
    if units=="1":
        cbar.set_label('$\Delta$ {}[fractional area]'.format(varname))
    else:
        cbar.set_label('$\Delta$ {}[{}]'.format(varname,units))
    #generating limits for plots
    try:
        #try to read them. May have already been generate
        lims = process.read_lims(varname,csvdir)
    except:
        #If they haven't been generated, generate then read them
        print "error with reading limits.. attemping to create limits"
        process.anom_limit_setup(varname,[2,9],["u-au866","u-au872","u-au874","u-av231"],csvdir)
        lims = process.read_lims(varname,csvdir)
    print lims
    plt.clim(float(lims["Min"]),float(lims["Max"]))
    plt.show()
    fig.savefig(
        '/home/ben/Desktop/{}-{}-{}'.format(modelname, varname, monthnum))
    plt.close()

def month_map_anom_volume(modelname, monthnum):
    """function called when anomaly map plots of variable varname are wanted..."""
    lons, lats, aice, total_diff_aice, units = grab.month_map_anom_test(
        "/media/windowsshare", modelname, monthnum, "aice", True)  # grabbing data
    print("total difference in variable {} is {}".format("aice",total_diff_aice))
    lons, lats, sithick, total_diff_sithick, units = grab.month_map_anom_test(
        "/media/windowsshare", modelname, monthnum, "sithick", False)
    lonstest, latstest, myvartest, tarea = grab.month_map_test("/media/windowsshare",modelname,"aice")
    concentration = np.multiply(aice,tarea)
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawmapboundary(linewidth=1)
    m.drawlsmask(land_color='grey', ocean_color='grey', lakes=True)
    cm = m.pcolormesh(lons, lats, concentration, latlon=True, cmap='jet')
    monthdict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    plt.title("{} total sea ice concentration anomaly in the month of {}".format(modelname,monthdict[monthnum]))
    plt.clim([-0.2,0.2])
    fig.savefig(
        '/home/ben/Desktop/totalanom_{}_{}.png'.format(modelname,monthnum))
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
    # first stripping modelvar and controlvar of spatial data.. converting them into a sequence

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
    lons, lats, modelvar = grab.ice_area_map_mean(
        "/media/windowsshare",modelname,monthnum)
    lons1, lats1, controlvar = grab.ice_area_map_mean(
        "/media/windowsshare","u-at053",monthnum)
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
    plt.title("scatterplot comparing variable aice*tarea for {} with control".format(varname,modelname))
    plt.xlabel(modelname)
    plt.ylabel("u-at053")
    #plt.xlim([0.0,2.5e9])
    #plt.ylim([0.0,2.5e9])
    plt.plot(modelvar_seq,modelvar_seq,"r")
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

def testplot(modelname,varname):
    """function called to test whether or not grabbing functions are grabbing correct data, comparing to Panoply plots..."""
    lons, lats, myvar = grab.month_map_test(
        "/media/windowsshare",modelname,varname)  # grabbing data
    #now saving limits of plot to csv file so that month_map_anom_main can use them
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawcoastlines(linewidth=1)
    m.drawlsmask(land_color='grey',ocean_color='grey',lakes=True)
    m.drawmapboundary(linewidth=1)
    cm = m.pcolormesh(lons, lats, myvar, latlon=True)
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    plt.show()
    plt.close()

def plot_area_main(modelname, monthnum, outputdir,xlim,ylim):
    """plots the selected snow + ice thickness for the given area"""
    #grabbing sea ice thickness
    lons, lats, sithick, units = grab.month_map_mean(
        "/media/windowsshare", modelname, monthnum, "hi",True)
    #grabbing snow thickness 
    lons, lats, sisnthick, units = grab.month_map_mean(
        "/media/windowsshare", modelname, monthnum, "hs",True)
    totalthick= sithick + sisnthick
    print np.ma.mean(sithick)
    print np.ma.mean(sisnthick)
    
    #now we will make the plot
    fig, ax = plt.subplots(figsize=(8, 8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawmapboundary(linewidth=0.3)
    cm = m.pcolormesh(lons, lats, totalthick, latlon=True, cmap='gist_rainbow')
    cbar = m.colorbar(cm, location='bottom', pad="5%")
    m.shadedrelief()
    plt.clim([0.0,4.0])
    plt.xlim([xlim[0],xlim[1]])
    plt.ylim([ylim[0],ylim[1]])
    plt.title("hi + hs in selected area for {}".format(modelname))
    plt.show()
    #fig.savefig("/home/ben/Desktop/{}_pat".format(modelname))

def plot_all(modelname,monthnum,var1,varcond,xlim,ylim,var2="aice"):
    """plots all files of a model for a given variable and saves to an output directory. Optionally plots varname+varname2 if varcond==True"""
    if varcond:
        lons,lats,vars1,varnames1 = grab.month_map_data("/media/windowsshare",modelname,monthnum,var1)
        lons,lats,vars2,varnames2 = grab.month_map_data("/media/windowsshare",modelname,monthnum,var2)
        for i in xrange(20):
            var1_slab = vars1[i,:,:]
            var2_slab = vars2[i,:,:]
            total_slab = var1_slab + var2_slab
            #now we will make the plot
            fig, ax = plt.subplots(figsize=(8, 8))
            m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
            m.drawmapboundary(linewidth=0.3)
            cm = m.pcolormesh(lons, lats, total_slab, latlon=True, cmap='gist_rainbow')
            cbar = m.colorbar(cm, location='bottom', pad="5%")
            m.shadedrelief()
            plt.clim([0.0,4.0])
            plt.xlim([xlim[0],xlim[1]])
            plt.ylim([ylim[0],ylim[1]])
            plt.title("plot of total of variables {}+{}\n for file {} \n month# {}".format(var1,var2,varnames1[i],monthnum))
            print "saving picture number {}/20".format(i+1)
            fig.savefig("/home/ben/Desktop/{}_{}_{}".format(modelname,monthnum,i))
            plt.close()

def regrid(arr1,lats1,lons1,arr2,lats2,lons2,modelname,monthstr):
    """takes two netCDF arrays and their respective latitudes and longitudes and regrids 
    both arrays to be on the same grid. It then plots them and saves the plot, masking the areas where the grid does not overlap."""
    #NOTE!!!!!!!! Here arr1 should be the model data and arr2 should be the NSIDC data (this is because the model has coarser resolution
    # in the area that we are plotting over. 
    
    #NEW VERSION
    fig,ax=plt.subplots(figsize=(8,8))
    m = Basemap(resolution='h', projection='spstere',
                lat_0=-90, lon_0=-180, boundinglat=-55)
    m.drawmapboundary(fill_color='grey')
    
    #projecting both arrays onto grid
    #make mesh grid of latlons for plotting
    ny=arr1.shape[0]; nx=arr1.shape[1]
    print("shape of model data is {},{}".format(nx,ny))
    glons,glats = m.makegrid(nx,ny)

    #interp data of real world to this reg grid 
    gdata1 = griddata((lons1[~arr1.mask].ravel(),lats1[~arr1.mask].ravel()),
                    arr1[~arr1.mask].ravel(),(glons,glats),method='cubic')
    
    gdata2 = griddata((lons2[~arr2.mask].ravel(),lats2[~arr2.mask].ravel()),arr2[~arr2.mask].ravel(),(glons,glats),method='cubic')
    
    # make mask of missing data (ie land)
    mask = np.zeros((ny,nx))
    mask[arr1.mask] = 1
    gmask = griddata((lons1.ravel(),lats1.ravel()),mask.ravel(),
                    (glons,glats),method='cubic')

    pdat = np.ma.masked_array(gdata1-gdata2,gmask>0.5) # masked, gridded array

    cs=m.pcolormesh(glons,glats,pdat,latlon=True,cmap='seismic')
    plt.clim(-1.0,1.0)
    m.drawcoastlines()
    cbar = m.colorbar(cs,location='bottom',pad="5%",extend='both')
    plt.title("Mean ice concentration anomaly\n Between model scheme {} and NSIDC data in the month of {}".format(modelname,monthstr))
    plt.show()
    fig.savefig("/home/ben/Desktop/{}".format(modelname))
    plt.close()
    

if __name__=="__main__":
    models = ["u-au866","u-av231"]
#month_map_mean_main(models[0],2,"aice","/home/ben/Desktop/summer2019/plotlims",True)
