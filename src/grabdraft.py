rly(path,modelname,month):
    """imports variables from NetCDF files with specified path and variable name"""
    os.chdir("../../../../")
    os.chdir("{}/{}/{}".format(path,modelname,"ice"))
    filecount=0
    #simply counting number of files in directory.
    for filename in os.listdir('./'):
        filecount +=1
    #now we will sort the files based on month.
    monthdata = []
    stddevs = []
    for filenum,filename in enumerate(os.listdir('./')):
        #grabbing month value from filename... format does not vary
        monthstr = filename[19:21]
        thismonth = int(monthstr)
        if month==thismonth:
            print "Grabbing {}, file {} of {}".format(filename,filenum,filecount)
            aice = np.ma.squeeze(np.ma.array(testdata.variables['aice'][:,:],dtype='float64'))
        tarea = np.ma.array(testdata.variables['tarea'][:,:],dtype='float64')
        #adding total ice area into its respective month bin
        monthareas[month].append(np.ma.sum(aice*tarea))
        print "area added is {}".format(np.ma.sum(aice*tarea))
        testdata.close()
    #now calculating the mean for each month and returning that value, as well as the standard deviation for each month.
    for i,month in enumerate(monthareas):
        stddevs[i] = stats.stdev(month)
        means[i] = stats.mean(month)
    return stddevs, means
