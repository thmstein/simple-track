#!/usr/local/sci/bin/python2.7

from netCDF4 import Dataset as ncfile
import numpy as np
import matplotlib.pyplot as plt

###################################################
# loadfile IS A USER SPECIFIED FUNCTION TO LOAD THE DATA AND TIME STAMP INFORMATION
# OUTPUT
# datad = data (2D array)
# fidd = file time identifier yyyymmdd
# hh = file hour
# mm = file minute stamp
###################################################

def loadfile(filename):
    nc = ncfile(filename)
    datad = nc.variables['var'][200:600,250:550]/32
    datad = np.flipud(np.transpose(datad))
    fidd = filename[-9:-5]
    hh = float(fidd[0:2])
    mm = float(fidd[2:4])
    
    return datad, fidd, hh, mm

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
    lrain[np.where(lrain<=0.)]=0.01
        
    figa=plt.figure(figsize=(6, 7))
    #ax = figa.add_subplot(111)
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    con = plt.pcolor(xmat,ymat,np.log2(lrain),vmin=-1,vmax=5)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/7,width,width*6/7]
    plt_ax.set_position(posnew)
    plt.xlabel('Distance from Chilbolton [km]')
    plt.ylabel('Distance from Chilbolton [km]')
    colorbar_axes = figa.add_axes([left, bottom, width, 0.01])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Rainfall rate [log2 mm hr^{-1}]')
    plt.savefig(IMAGES_DIR + 'Rainrate_' + write_file_ID + '.png')
    plt.close()
        
    figb=plt.figure(figsize=(6, 7))
    #ax = figa.add_subplot(111)
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    con = plt.pcolor(xmat,ymat,wasarray,vmin=-10,vmax=200)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/7,width,width*6/7]
    plt_ax.set_position(posnew)
    plt.xlabel('Distance from Chilbolton [km]')
    plt.ylabel('Distance from Chilbolton [km]')
    colorbar_axes = figb.add_axes([left, bottom, width, 0.01])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Storm ID')
    plt.savefig(IMAGES_DIR + 'Stormid_' + write_file_ID + '.png')
    plt.close()
        
    lifearray[np.where(lifearray==0)]=-6
    figc=plt.figure(figsize=(6, 7))
    #ax = figa.add_subplot(111)
    #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
    con = plt.pcolor(xmat,ymat,5*lifearray,vmin=-30,vmax=60)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    posnew=[left,bottom+height/7,width,width*6/7]
    plt_ax.set_position(posnew)
    plt.xlabel('Distance from Chilbolton [km]')
    plt.ylabel('Distance from Chilbolton [km]')
    colorbar_axes = figc.add_axes([left, bottom, width, 0.01])
    # add a colourbar with a label
    cbar = plt.colorbar(con, colorbar_axes, orientation='horizontal')
    cbar.set_label('Life time [mins]')
    plt.savefig(IMAGES_DIR + 'Lifetime_' + write_file_ID + '.png')
    plt.close()
    if do_vectors==True:
        figd=plt.figure(figsize=(6, 7))
        #ax = figa.add_subplot(111)
        #con = ax.imshow(f, cmap=cm.jet, interpolation='nearest')
        con = plt.contour(xmat,ymat,lrain,levels=[threshold])
        plt.quiver(xmat[::10,::10],ymat[::10,::10],newumat[::10,::10]/num_dt,newvmat[::10,::10]/num_dt,pivot='mid',units='width')
        plt_ax = plt.gca()
        left, bottom, width, height = plt_ax.get_position().bounds
        posnew=[left,bottom+height/7,width,width*6/7]
        plt_ax.set_position(posnew)
        plt.xlabel('Distance from Chilbolton [km]')
        plt.ylabel('Distance from Chilbolton [km]')
        plt.savefig(IMAGES_DIR + 'Vectors_' + write_file_ID + '.png')
        plt.close()

def plot_poi(HTML_DIR, DATA_DIR, filelist, poi_array, poi_boxdown, poi_boxleft, poi_boxup, poi_boxright, maxgifs, padt, padxy, poi_increasing, poi_string, xmathtml, ymathtml):
	if poi_increasing:
		poi_max_list = np.nanmax(poi_array,0) # Find maximum value for each storm (will need to allow for minima for instance in case of brightness temperatures
		poi_max_list[np.where(np.isnan(poi_max_list)==1)] = 0.
		poi_indices = np.argsort(poi_max_list) # Sorts maxima from small to large and returns indices
	else:
		poi_min_list = np.nanmin(poi_array,0)
		poi_min_list[np.where(np.isnan(poi_min_list)==1)] = 999.
		poi_indices = np.argsort(poi_min_list) 
	for gifnum in range(0,maxgifs):
		if gifnum>9:
			gifstring = str(gifnum)
		else:
			gifstring = '0' + str(gifnum)
		if poi_increasing==1:
			poi_index = poi_indices[-1-gifnum]
		else:
			poi_index = poi_indices[gifnum]
		# Establish times of interest for this storm, adding 1 hour each side of the storm period.
		timeind=np.where(np.isnan(poi_array[:,poi_index])==0)
		mintime = np.max([np.min(timeind)-padt,0])
		maxtime = np.min([np.max(timeind)+padt,len(filelist)])
		# Establish domain for each storm, padding it with some sensible amount on all sides 
		minx = np.max([np.nanmin(poi_boxleft[:,poi_index])-padxy,np.min(xmathtml)])
		maxx = np.min([np.nanmax(poi_boxright[:,poi_index])+padxy,np.max(xmathtml)])
		miny = np.max([np.nanmin(poi_boxdown[:,poi_index])-padxy,np.min(ymathtml)])
		maxy = np.min([np.nanmax(poi_boxup[:,poi_index])+padxy,np.max(ymathtml)])
		# Set domain x and y limits to ensure a square domain
		if maxy-miny > maxx-minx:
			maxx = maxx + (maxy-miny-(maxx-minx))/2
			minx = minx - (maxy-miny-(maxx-minx))/2
		else:
			maxy = maxy + (maxx-minx-(maxy-miny))/2
			miny = miny - (maxx-minx-(maxy-miny))/2
		# Load relevant arrays and plot storms with contours
		for tij in range(mintime,maxtime):
			var,file_ID,hourval,minval = loadfile(DATA_DIR + filelist[tij])
			if hourval<10:
				hourstr = '0' + str(hourval)
			else:
				hourstr = str(hourval)
			if minval<10:
				minstr = '0' + str(minval)
			else:
				minstr = str(minval)
			wasfile = HTML_DIR + 'wasarray' + file_ID + '.npy'
			wasarray = np.load(wasfile)
			wasmask = np.where(wasarray!=poi_index,0,1)
			fig = plt.figure(figsize=(6,6))
			lrain=var+0.0
			lrain[np.where(lrain<=0.)]=0.01
			plt.pcolormesh(xmathtml,ymathtml,np.log2(lrain),vmin=-1,vmax=5) 
			plt.contour(xmathtml,ymathtml,wasmask,levels=[0.5,2.],colors='red',linewidths=2)
			plt.xlim(minx,maxx)
			plt.ylim(miny,maxy)
			plt.xlabel('Distance from Chilbolton [km]')
			plt.ylabel('Distance from Chilbolton [km]')
			plt.title('Object ' + gifstring + ' at ' + hourstr + minstr)
			plt.savefig(HTML_DIR + poi_string + '/' + 'poi_' + poi_string + '_' + gifstring + '_' + file_ID)
			plt.close()

