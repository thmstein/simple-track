import object_tracking
import numpy as np
import datetime
import os
import user_functions

from scipy import interpolate 

##################################################################
# WHAT PERIOD DO YOU WANT TO RUN TRACKING FOR?
##################################################################

startyear = 2021
startmonth = 1
startday = 1
starthour = 0
#startmin = 0

totalhours = 24

start_time = datetime.datetime(startyear,startmonth,startday,starthour,0,0,0)

##################################################################
# BRIEFLY LOAD FIRST HIMAWARI FILE FOR LATS AND LONS
# NEEDED TO INTERPOLATE GPM
##################################################################

# NB squarelength is important for correlations later
squarelength = 50. ## Size in pixels of individual squares to run fft for (dx,dy) displacement. Must divide (x,y) lengths of the array!

IMAGES_DIR = './output/'
DATA_DIR = './himawari/java/13/2021/01/01/'
filelist = os.listdir(DATA_DIR)
filelist = np.sort(filelist)
firsttime=True
var,file_ID,yearval,monthval,dayval,hourval,minval,xx,yy = user_functions.loadfile(DATA_DIR + filelist[0],firsttime)
ylen=np.size(yy)/np.int(squarelength)
xlen=np.size(xx)/np.int(squarelength)
xmat,ymat = np.meshgrid(range(0,xlen*np.int(squarelength)),range(0,ylen*np.int(squarelength)))
lat = yy[0:ylen*np.int(squarelength)]
lon = xx[0:xlen*np.int(squarelength)]

##################################################################
# THE FOLLOWING GENERATES GPM RAINFALL ACCUMULATION FOR THE PERIOD OF INTEREST.
# THIS CAN BE USED TO IDENTIFY REGIONS OF INTEREST, E.G. RAINFALL MAXIMUM, TO FOCUS ANIMATIONS.
##################################################################

GPM_DIR = './gpm/'
IMAGES_DIR = './output/'
filelist = os.listdir(GPM_DIR)
filelist = np.sort(filelist)
firsttime = True

for nt in range(len(filelist)):
	# Load new image
	var,yearval,monthval,dayval,hourval,minval,xx,yy = user_functions.loadgpm(GPM_DIR + filelist[nt],firsttime)
	now_time = datetime.datetime(yearval,monthval,dayval,hourval,minval,0,0)
	delta_t = now_time - start_time
	elap = 24*delta_t.days + delta_t.seconds/3600.
	if elap > totalhours:
	    break
	if firsttime==True:
	    totalrain = var/2.
	    theLats = yy
	    theLons = xx
	    firsttime=False
	else:
	    totalrain = totalrain + var/2.	

latmin=-15
latmax=0
lonmin=100
lonmax=120

inda = np.where((theLons>lonmin) & (theLons<lonmax))
indb = np.where((theLats>latmin) & (theLats<latmax))
newrain = totalrain[indb[0][0]:indb[0][-1]+1,inda[0][0]:inda[0][-1]+1]
x, y = np.float32(np.meshgrid(theLons[inda], theLats[indb]))
xhim, yhim = np.float32(np.meshgrid(lon, lat))

himrain=interpolate.griddata((x.flatten(),y.flatten()),newrain.flatten() , (xhim,yhim),method='linear')

#x, y = np.float32(np.meshgrid(theLons,theLats))
user_functions.plot_gpm(xhim,yhim,himrain,latmin,latmax,lonmin,lonmax,IMAGES_DIR)

##################################################################
# THE FOLLOWING PARAMETERS SHOULD BE CHANGED BASED ON THE DATA (RESOLUTION ETC.)
##################################################################

# Integer number (dimensions user-defined) to identify MINIMUM time difference between consecutive data files
# Example 1: Radar data 5-minutes apart with time stamp in filename, dt = 5
# Example 2: Satellite brightness temperatures hourly with time stamp in filename, dt = 1
# NB. When writing storms, (dx,dy) will have units PIXELS per TIME STEP (specified by dt), so already scaled by number of missing files
dt = 20. 
dt_tolerance = 60. # Maximum separation in time allowed between consecutive images

under_t=True ## True = labelling areas *under* the threshold (e.g. brightness temperature), False = labelling areas *above* threshold (e.g. rainfall)
threshold = 233. ## Threshold used to identify objects (with value of variable greater than this threshold)
minpixel = 4. ## Minimum object size in pixels
rafraction = 0.01 ## Minimum fractional cover of objects required for fft to obtain (dx,dy) displacement
dd_tolerance = 3. # Maximum difference in displacement values between adjacent squares (to remove spurious values) - scaled by num_dt if necessary
halopixel = 5. ## Radius of halo in pixels for orphan storms - big halo assumes storms may spawn "children" at a distance multiple pixels away

flagwrite = True ## For writing storm history data in a text file 
doradar = False ## doradar=True is calculate range and azimuth for real-time tracking with radar (e.g. Chilbolton). doradar=False any other use, radar coordinates not relevant
misval = -999 ## Missing value
struct2d = np.ones((3,3)) ## np.ones((3,3)) is 8-point connectivity for labelling storms. Can be changed to user preference.

flagplot = True ## For plotting images (vectors and IDs). Also set plot_type...
flagplottest = False ## For plotting fft correlations (testing only, very slow, lots of plots)

if flagplot or flagplottest:
    plot_type = '.png'
    if plot_type == '.eps':
        my_dpi_global = 150 # for rasterized plots, this controls the size and quality of the final plot
    elif plot_type == '.png':
        my_dpi_global = 300

##################################################################
# THE FOLLOWING PARAMETERS CAN BE CHANGED, BUT SHOULD NOT BE
##################################################################

lapthresh = 0.6 ## Minimum fraction of overlap (0.6 in TITAN)

##################################################################
# AUTOMATIC SET UP OF DISCRETE VALUES BASED ON USER INPUT PARAMETERS
# ESSENTIAL - NO NEED TO CHANGE THESE
##################################################################
# squarehalf: To determine grid spacing for coarse grid of (dx,dy) estimates
# areastr: For filename identifier of area threshold used
# thr_str: For filename identifier of variable threshold used
# fftpixels: Minimum number of thresholded pixels needed to calculate (dx,dy) 
# halosq: To identify if new cell is nearby existing cell
##################################################################

squarehalf = int(squarelength/2)   
areastr = str(int(minpixel))
thr_str = str(int(threshold))
sql_str = str(int(squarelength))
fftpixels = squarelength**2/int(1./rafraction)
halosq = halopixel**2

##################################################################
# AUTOMATIC SET UP OF TEXT STRING FOR INFORMATION ON LABELLING
# NOT ESSENTIAL
##################################################################

label_method='Brightness Temperature < ' + thr_str + 'K'    

#################################################################
# For test data, try the following:
# 
# All data (5-minute intervals)
#DATA_DIR = './data/'
#IMAGES_DIR = './output/'
# Sparse data (10-minute intervals), to test similarity in vector fields (scaling by num_dt working correctly)
#DATA_DIR = './data/'
#IMAGES_DIR = './output/'
# Missing data (10-minute intervals, 1 file missing), to test dt_tolerance
#DATA_DIR = './data/'
#IMAGES_DIR = './output/'
#################################################################
filelist = os.listdir(DATA_DIR)
filelist = np.sort(filelist)
if doradar:
	rarray=np.sqrt(xmat**2+ymat**2);
	azarray=np.arctan(xmat/ymat);
	azarray[np.where((xmat<0) & (ymat>=0))]=azarray[np.where((xmat<0) & (ymat>=0))]+2*np.pi;
	azarray[np.where(ymat<0)]=azarray[np.where(ymat<0)]+np.pi; 
	azarray=180*azarray/np.pi;    
	azarray[np.where(np.isnan(azarray)==1)]=0         

#   Initialise variables
OldData, OldLabels, oldvar, newvar, prev_time = [], [], [], [], []
newwas = 1
plot_vectors = False

oldhourval = []
oldminval = []
oldmask = []
newmask = []
num_dt = []
firsttime = True

for nt in range(len(filelist)):
	# Load new image
	var,file_ID,yearval,monthval,dayval,hourval,minval,xx,yy = user_functions.loadfile(DATA_DIR + filelist[nt],firsttime)
	now_time = datetime.datetime(yearval,monthval,dayval,hourval,minval,0,0)
	delta_t = now_time - start_time
	elap = 24*delta_t.days + delta_t.seconds/3600.
	if elap > totalhours:
	    break
	if firsttime==True:
	    ylen=np.size(yy)/np.int(squarelength)
	    xlen=np.size(xx)/np.int(squarelength)
	    xmat,ymat = np.meshgrid(range(0,xlen*np.int(squarelength)),range(0,ylen*np.int(squarelength)))
	    lat = yy[0:ylen*np.int(squarelength)]
	    lon = xx[0:xlen*np.int(squarelength)]
	    #var=var[0:xlen*np.int(squarelength),0:ylen*np.int(squarelength)]
	    #print(np.shape(xmat))
	    #print(np.shape(ymat))
	    #print(np.shape(var))
	var=var[0:ylen*np.int(squarelength),0:xlen*np.int(squarelength)]
	hourval=float(hourval)
	minval=float(minval)
	firsttime=False
	print(file_ID)
	write_file_ID = 'S' + sql_str + '_T'+ thr_str +'_A'+ areastr +'_'+ file_ID
	NewLabels=object_tracking.label_storms(var,minpixel,threshold,struct2d,under_t)
	# oldmask, newmask, USED FOR DERIVING (dx,dy)
	# THESE CAN BE CHANGED USING EXPERT KNOWLEDGE (e.g. use raw data rather than binary masks, if displacement information is contained in structures within objects)
	# !!! NB If raw data are used (i.e. not zeros and ones) then fftpixels needs to be changed to remain sensible !!!
	if len(OldLabels) > 1:
	    # CHECK TIME DIFFERENCE BETWEEN CONSECUTIVE IMAGES 
		dtnow = user_functions.timediff(oldhourval,oldminval,hourval,minval)
		num_dt = dtnow/dt
		if dtnow > dt_tolerance:
			print('Data are too far apart in time --- Re-initialise objects')
			OldData, OldLabels, oldvar, newvar, prev_time = [], [], [], [], []
			newwas = 1
			plot_vectors = False
			continue
		oldmask = np.where(OldLabels>=1,1,0)
		newmask = np.where(NewLabels>=1,1,0)
	# Call object tracking routine
	# NewData = list of objects and properties
	# newwas = final label number
	# NewLabels = array with object IDs from [1, nummax] as found by label_storms
	# newumat, newvmat = arrays with (dx,dy) displacement between two images (NB not displacement per dt!!!) 
	# wasarray = array with object IDs consistent across images (i.e. tracked IDs)
	# lifearray = array with object lifetime consistent across images
	NewData, newwas, NewLabels, newumat, newvmat, wasarray, lifearray = object_tracking.track_storms(OldData, var, newwas, NewLabels, OldLabels, xmat, ymat, fftpixels, dd_tolerance, halosq, squarehalf, oldmask, newmask, num_dt, lapthresh, misval, doradar, under_t, IMAGES_DIR, write_file_ID, flagplottest)
	# Write tracked storm information (see object_tracking.write_storms)
	if flagwrite:
		object_tracking.write_storms(write_file_ID, start_time, now_time, label_method, squarelength, rafraction, newwas, NewData, doradar, misval, IMAGES_DIR)
	# Plot tracked storm information (see user_functions.plot_example)
	if flagplot:
		#user_functions.plot_example(write_file_ID, nt, var, lon, lat, newumat, newvmat, num_dt, wasarray, lifearray, x0array, y0array, threshold, IMAGES_DIR, plot_vectors)
		numstorms = NewLabels.max()
		x0array=[]
		y0array=[]
		x1array=[]
		y1array=[]
		for sij in range(numstorms):
			x0array.append(lon[np.int(NewData[sij].centroidx0)])
			x1array.append(lon[np.int(NewData[sij].centroidx)])
			y0array.append(lat[np.int(NewData[sij].centroidy0)])
			y1array.append(lat[np.int(NewData[sij].centroidy)])
		user_functions.plot_spisea(write_file_ID, nt, var, lon, lat, newumat, newvmat, num_dt, wasarray, lifearray, x0array, y0array, x1array, y1array, threshold, xhim, yhim, himrain, latmin,latmax,lonmin,lonmax, IMAGES_DIR, plot_vectors)
	# Save tracking information in preparation for next image
	OldData = NewData
	OldLabels = NewLabels
	oldvar = var
	oldhourval = hourval
	oldminval = minval
	plot_vectors = True

##################################################################
# TRACKING IS NOW FINISHED.
##################################################################

