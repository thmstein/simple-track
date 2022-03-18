#!/usr/local/sci/bin/python2.7

from netCDF4 import Dataset as ncfile
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import h5py as h5py
import cartopy.feature as cfeat

###################################################
# loadfile IS A USER SPECIFIED FUNCTION TO LOAD THE DATA AND TIME STAMP INFORMATION
# OUTPUT
# datad = data (2D array)
# fidd = file time identifier yyyymmdd
# hh = file hour
# mm = file minute stamp
###################################################

def loadfile(filename,firsttime):
    nc = ncfile(filename)
    datad = nc.variables['T_b'][0,:,:]
    fidd = filename[-16:-3]
    rr = int(fidd[-13:-9])
    oo = int(fidd[-9:-7])
    dd = int(fidd[-7:-5])
    hh = int(fidd[-4:-2])
    mm = int(fidd[-2::])
    xx=[]
    yy=[]
    if firsttime==True:
        xx = nc.variables['longitude'][:]
        yy = nc.variables['latitude'][:]
    
    return datad, fidd, rr, oo, dd, hh, mm, xx, yy

###################################################
# loadgpm IS A USER SPECIFIED FUNCTION TO LOAD THE GPM DATA AND TIME STAMP INFORMATION
# OUTPUT
# datad = rainfall (2D array) in mm/hr
# fidd = file time identifier yyyymmdd
# hh = file hour
# mm = file minute stamp
###################################################

def loadgpm(filename,firsttime):
    dataset = h5py.File(filename, 'r')
    datad = dataset['Grid/precipitationCal'][:]
    datad = np.transpose(datad)
    datad = datad[:,:,0]
    rr = int(filename[-39:-35])
    oo = int(filename[-35:-33])
    dd = int(filename[-33:-31])
    hh = int(filename[-29:-27])
    mm = int(filename[-27:-25])
    xx=[]
    yy=[]
    if firsttime==True:
        xx = dataset['Grid/lon'][:]
        yy = dataset['Grid/lat'][:]
    
    return datad, rr, oo, dd, hh, mm, xx, yy

###################################################
# timediff IS A USER SPECIFIED FUNCION TO CALCULATE TIME SEPARATION BETWEEN CONSECUTIVE IMAGES
# OUTPUT
# tdif = time difference in units relevant to user specification (to be divided by "dt" in wrapper.py)
###################################################

def timediff(oldh,oldm,newh,newm):
    hdif=newh-oldh
    mdif=newm-oldm
    tdif=60.*hdif+mdif
    
    return tdif

###################################################
# plot_example IS ONLY USED AS AN ILLUSTRATION
# OF THE EXAMPLE DATA
###################################################

def plot_example(write_file_ID, nt, rain, xmat, ymat, newumat, newvmat, num_dt, wasarray, lifearray, threshold, IMAGES_DIR, do_vectors):
    '''
    PLOT FIGURES WITH RAINFALL RATE AND STORM LABELS
    FOR ILLUSTRATIVE AND TESTING PURPOSES
    '''
        
    lrain=rain+0.0
    maskwas = ma.masked_array(wasarray, mask=np.where(wasarray==0,1,0))
    masklife = ma.masked_array(lifearray, mask=np.where(lifearray==0,1,0))
    #lrain[np.where(lrain<=0.)]=0.01
        
    figb=plt.figure(figsize=(6, 5))
    #ax = figa.add_subplot(111)
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    con = plt.pcolor(xmat,ymat,maskwas,vmin=0,vmax=2000)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    colorbar_axes = figb.add_axes([left, bottom, width, 0.05])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Storm ID')
    plt.savefig(IMAGES_DIR + 'Stormid_' + write_file_ID + '.png')
    plt.close()

    figa=plt.figure(figsize=(12, 5))
    #ax = figa.add_axes()
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    axa=plt.axes(projection=ccrs.PlateCarree(),position=[0.125,0.11,0.352,0.77])
    axa.coastlines()
    con = plt.pcolor(xmat,ymat,lrain,vmin=180,vmax=310)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    colorbar_axes = figa.add_axes([left, bottom, width, 0.05])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('10.8 Brightness Temperature [K]')
    #plt.savefig(IMAGES_DIR + 'Temperature_' + write_file_ID + '.png')
    #plt.close()
        
    #figc=plt.figure(figsize=(6, 5))
    #ax2 = figa.add_axes()
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    axb=plt.axes(projection=ccrs.PlateCarree(),position=[0.548,0.11,0.352,0.77])
    axb.coastlines()
    con = plt.pcolor(xmat,ymat,20.*masklife/60.,vmin=0,vmax=12)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    colorbar_axes = figa.add_axes([left, bottom, width, 0.05])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Life time [hr]')
    plt.savefig(IMAGES_DIR + 'Lifetime_' + write_file_ID + '.png')
    plt.close()
    
    if do_vectors==True:
        figd=plt.figure(figsize=(6, 5))
        #ax = figa.add_subplot(111)
        #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
	ax = plt.axes(projection=ccrs.PlateCarree())
	ax.coastlines()
	con = plt.contour(xmat,ymat,lrain,levels=[threshold])
        plt.quiver(xmat[::10],ymat[::10],newumat[::10,::10]/num_dt,newvmat[::10,::10]/num_dt,pivot='mid',units='width')
        plt_ax = plt.gca()
        left, bottom, width, height = plt_ax.get_position().bounds
        posnew=[left,bottom+height/5,width,height*4/5]
        plt_ax.set_position(posnew)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.savefig(IMAGES_DIR + 'Vectors_' + write_file_ID + '.png')
        plt.close()

###################################################
# plot_gpm IS A USER SPECIFIED FUNCTION TO PLOT THE GPM ACCUMULATION OVER THE TRACKING PERIOD
###################################################

def plot_gpm(lons,lats,accumulation,latmin,latmax,lonmin,lonmax,IMAGES_DIR):

     # Plot the figure, define the geographic bounds

     ind = np.where((lats > latmin) & (lats < latmax) & (lons > lonmin) & (lons < lonmax))

     accumulation[np.where(accumulation<=0.1)]=0.1
     masked_array = np.ma.masked_where(accumulation <= 0.1,accumulation)
     maxval = np.nanmax(masked_array[ind])
     maxint = int(np.log2(maxval)*2)
     
     bluemask = np.copy(accumulation)
     bluemask[np.where(accumulation<maxval/4.)]=0.1
          
     clevs = np.arange(0,maxint/2+.5,.25)
     base2 = 2**clevs[::4]
     
     figx = plt.figure(figsize=(6, 5))
     ax = plt.axes(projection=ccrs.PlateCarree())
     ax.coastlines()
     #cmap = plt.get_cmap()
     #cmap.set_bad('gray')
     
     # Plot every masked value as white
     #cs = plt.pcolormesh(lons,lats,np.log2(masked_array),vmin=0,vmax=maxint/2.+.5)
     extent = np.min(lons), np.max(lons), np.min(lats), np.max(lats)
     cs = plt.imshow(np.log2(np.flipud(accumulation)),vmin=0,vmax=maxint/2.+.5,cmap=plt.cm.Greys,extent=extent)
     bs = plt.contourf(lons,lats,np.log2(bluemask),levels=clevs,cmap=plt.cm.viridis)
     
     plt.xlabel('Longitude')
     plt.ylabel('Latitude')
     plt.xlim(lonmin,lonmax)
     plt.ylim(latmin,latmax)

     plt_ax = plt.gca()
     left, bottom, width, height = plt_ax.get_position().bounds
     posnew=[left,bottom+height/5,width,height*4/5]
     plt_ax.set_position(posnew)
     colorbar_axes = figx.add_axes([left, bottom, width, 0.03])
     colorbar_axes2 = figx.add_axes([left, bottom+0.035, width, 0.03])
     # add a colourbar with a label
     abar = plt.colorbar(bs,colorbar_axes2, orientation='horizontal',ticks=[])
     cbar = plt.colorbar(cs, colorbar_axes, orientation='horizontal', ticks=clevs[::4])
     cbar.ax.set_xticklabels(base2)
     cbar.set_label('Total rainfall [mm]')
     
     plt.savefig(IMAGES_DIR + 'IMERG.png')
     plt.close()
     # Set the title and fonts

###################################################
# plot_example IS ONLY USED AS AN ILLUSTRATION
# OF THE EXAMPLE DATA
###################################################

def plot_spisea(write_file_ID, nt, rain, xmat, ymat, newumat, newvmat, num_dt, wasarray, lifearray, x0array, y0array, x1array, y1array, threshold, lons, lats, accumulation, latmin,latmax,lonmin,lonmax, IMAGES_DIR, do_vectors):
    '''
    PLOT FIGURES WITH RAINFALL RATE AND STORM LABELS
    FOR ILLUSTRATIVE AND TESTING PURPOSES
    '''
        
    lrain=rain+0.0
    maskwas = ma.masked_array(wasarray, mask=np.where(wasarray==0,1,0))
    masklife = ma.masked_array(lifearray, mask=np.where(lifearray==0,1,0))
    #lrain[np.where(lrain<=0.)]=0.01
        
    figb=plt.figure(figsize=(6, 5))
    #ax = figa.add_subplot(111)
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    con = plt.pcolor(xmat,ymat,maskwas,vmin=0,vmax=2000)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks(np.arange(80,140,4))
    plt.yticks(np.arange(-12,10,2))
    plt.xlim(xmat.min(),xmat.max())
    plt.ylim(ymat.min(),ymat.max())
    colorbar_axes = figb.add_axes([left, bottom, width, 0.05])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Storm ID')
    plt.savefig(IMAGES_DIR + 'Stormid_' + write_file_ID + '.png')
    plt.close()

    figa=plt.figure(figsize=(12, 10))
    #ax = figa.add_axes()
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    axa=plt.axes(projection=ccrs.PlateCarree(),position=[0.125,0.555,0.352,0.385])
    axa.coastlines()
    con = plt.pcolor(xmat,ymat,lrain,vmin=180,vmax=310)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks(np.arange(80,140,2))
    plt.yticks(np.arange(-12,10,2))
    plt.xlim(xmat.min(),xmat.max())
    plt.ylim(ymat.min(),ymat.max())
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    colorbar_axes = figa.add_axes([left, bottom, width, 0.03])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('10.8 Brightness Temperature [K]')
    #plt.savefig(IMAGES_DIR + 'Temperature_' + write_file_ID + '.png')
    #plt.close()
        
    #figc=plt.figure(figsize=(6, 5))
    #ax2 = figa.add_axes()
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    axb=plt.axes(projection=ccrs.PlateCarree(),position=[0.548,0.555,0.352,0.385])
    axb.coastlines()
    con = plt.pcolor(xmat,ymat,20.*masklife/60.,vmin=0,vmax=12)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks(np.arange(80,140,2))
    plt.yticks(np.arange(-12,10,2))
    plt.xlim(xmat.min(),xmat.max())
    plt.ylim(ymat.min(),ymat.max())
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    colorbar_axes = figa.add_axes([left, bottom, width, 0.03])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Life time [hr]')

    #figc=plt.figure(figsize=(6, 5))
    #ax2 = figa.add_axes()
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    axc=plt.axes(projection=ccrs.PlateCarree(),position=[0.125,0.055,0.352,0.385])
    axc.coastlines()
    templrain=lrain+0.0
    templrain[np.where(lrain>233.)]=300. 
    con = plt.pcolor(xmat,ymat,templrain,vmin=180,vmax=310)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks(np.arange(80,140,2))
    plt.yticks(np.arange(-12,10,2))
    plt.xlim(xmat.min(),xmat.max())
    plt.ylim(ymat.min(),ymat.max())
    
    # Add lines for each storm
    for sij in range(len(x0array)):
		distance = np.sqrt((x0array[sij]-x1array[sij])**2 + (y0array[sij]-y1array[sij])**2)
		if distance > 1.:
			plt.plot([x0array[sij],x1array[sij]],[y0array[sij],y1array[sij]],'-k.')
    
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    colorbar_axes = figa.add_axes([left, bottom, width, 0.03])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Motion')

    #figc=plt.figure(figsize=(6, 5))
    #ax2 = figa.add_axes()
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')

    ind = np.where((lats > latmin) & (lats < latmax) & (lons > lonmin) & (lons < lonmax))

    accumulation[np.where(accumulation<=1.0)]=1.0
    #masked_array = np.ma.masked_where(accumulation <= 0.1,accumulation)
    maxval = np.nanmax(accumulation[ind])
    maxint = int(np.log2(maxval)*2)
     
    #bluemask = np.copy(accumulation)
    #bluemask[np.where(accumulation<maxval/4.)]=0.1
          
    clevs = np.arange(0,maxint/2+1.,.25)
    base2 = 2**clevs[::4]

    axd=plt.axes(projection=ccrs.PlateCarree(),position=[0.548,0.055,0.352,0.385])
    axd.coastlines()
    con = plt.contourf(lons,lats,np.log2(accumulation),levels=clevs,cmap=plt.cm.viridis)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks(np.arange(80,140,2))
    plt.yticks(np.arange(-12,10,2))
    plt.xlim(xmat.min(),xmat.max())
    plt.ylim(ymat.min(),ymat.max())
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/5,width,height*4/5]
    plt_ax.set_position(posnew)
    colorbar_axes = figa.add_axes([left, bottom, width, 0.03])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal', ticks=clevs[::4])
    cbar.ax.set_xticklabels(base2)
    cbar.set_label('Total rainfall [mm]')

    plt.savefig(IMAGES_DIR + 'lifetime/Lifetime_' + write_file_ID + '.png')
    plt.close()
    
    if do_vectors==True:
        figd=plt.figure(figsize=(6, 5))
        #ax = figa.add_subplot(111)
        #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
	ax = plt.axes(projection=ccrs.PlateCarree())
	ax.coastlines()
	con = plt.contour(xmat,ymat,lrain,levels=[threshold])
        plt.quiver(xmat[::10],ymat[::10],newumat[::10,::10]/num_dt,newvmat[::10,::10]/num_dt,pivot='mid',units='width')
        plt_ax = plt.gca()
        left, bottom, width, height = plt_ax.get_position().bounds
        posnew=[left,bottom+height/5,width,height*4/5]
        plt_ax.set_position(posnew)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.savefig(IMAGES_DIR + 'Vectors_' + write_file_ID + '.png')
        plt.close()


