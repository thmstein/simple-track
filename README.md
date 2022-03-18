# simple-track
Threshold-based object tracking algorithm for 2D data.
Casestudy branch is designed to for case study analysis.
The prepared parameters will download Himawari data for 1-3 January 2020 (himawari.py).
Data for a similar period will need to be downloaded from GPM (or elements related to GPM should be commented out).
The code will track storm objects in Himawari data for given area and brightness temperature thresholds.
A figure is produced after each time step, showing (1) brightness temperatures (2) storm object life times (3) storm object displacements (4) GPM total rainfall accumulation.

# file description

"wrapper.py" is the script that runs the tracking algorithm with the desired input and output directories and parameter choices. This should be the only Python file you need to modify. Relevant parameters are described below.

"object_tracking.py" is the Python program that actually does the tracking. This script should not need modification, unless the user wishes to store additional information about tracked objects. That additional information will need to be added to the "class StormS()" section and to the "write_storms" function in this script.

"user_functions.py" is the user-specified script to load files, calculate time differences, and plot output. Other user-specified functions should be added here.

"himawari.py" automatically downloads Himawari data for a specific date. Data are sliced to the domain specified and interpolated to a 4km horizontal resolution. 

# parameters

"wrapper.py" contains a set of parameters that all need changing in relation to the user preferences and data sets. 

Data-relevant parameters are:
* dt:    		The minimum separation in time possible between consecutive images
* dt_tolerance:	The maximum separation in time allowed between consecutive images (to allow for missing data)
* under_t:	Is the variable of interest smaller than a threshold (e.g. Brightness temperatures, under_t=True) or greater (e.g. Rainfall, under_t=False)

Tracking-relevant parameters are:
* struct2d:	Defines the neighbour-searching function, np.ones((3,3)) is 8-point connectivity.
* minpixel:	The minimum number of pixels for an object to be tracked
* squarelength:	The size in pixels of individual square regions for which displacement vectors will be calculated (should be large enough to cover several mid-sized objects)
* rafraction:	The minimum fractional cover of objects required in order for displacement vectors to be calculated
* dd_tolerance:	The maximum difference (in number of pixels) allowed between adjacent displacement vectors
* halopixel:	Radius of halo in pixels to look for orphaned objects (objects that do not overlap but that are within this radius of a "parent" will still be classed as "child" and to have spawned off the original object)
* lapthresh:	Minimum overlap fraction required for objects to be considered potentially the same between consecutive images [Default is 0.6]

Output-relevant parameters are:
* flagwrite:	If False, then no text files with object information is included in the output. [Default should be True]
* misval:		Preferred value to used for missing values.
* flagplot:		If True, a few images are included in the output (plotting function defined in "user_functions.py" [Trials should set this to True, long runs could set it to False to save time]
* flagplottest:	If True, numerous test images are included to check the displacement vector calculations [Default should be False]
* flaghtml: 	If True, arrays with object IDs are stored at each time step, arrays with object properties are kept in memory, and at the end of the run, images are generated for objects having the largest or smallest values for the property of interest. A separate bash script need to be run to generate a HTML page with animated gifs using these images.

# input

A function to read the data and a function to calculate time difference between consecutive files should be provided in "user_functions.py"

The directory of the input data should be specified in "wrapper.py" under DATA_DIR

Note that the corresponding (x,y) arrays should be defined in "wrapper.py" as xmat and ymat. 

Further changes may be necessary to the listing and date/time specifications in "wrapper.py"

# output

The tracking algorithm outputs a text file with a name related to the matching input filename and
containing details of the objects identified in each timestep, eg.

The number following storm/cell is the storm/cell id
* area:     the number of grid cells that met the threshold in this storm
* centroid: the centre of the rectangle around the storm or the centre of a cell [latix, lonix]
* box:      defines the rectangle around the storm [minlatix, minlonix, nlats, nlons]
* life:     the number of timesteps this storm has been seen minus 1, i.e. if this is zero the storm
          has only been seen in one image so far; this should increment in each successive timestep until the storm disappears.
		  If the storm was created by splitting from a parent, life will be the life of the parent when this storm was created
		  and will increment thereafer.
* dx,dy:       the velocity of this storm at this timestep as determined from pattern correlation
* meanv:    the average of the variable of interest for the object
* extreme:  the minimum value of the variable of interest for the object
* accreted: the storm ids that merged with this one at this timestep
* parent:   if this storm split from another storm at this timestep this is the id of the parent storm 
* child:    the ids of storms that split off this storm at this timestep.

Additional properties can be added by experienced users by editing "object_tracking.py" (see above).

Plots can be generated based on the output (e.g. in "user_functions.py" see plot_example function) but this will slow down the code significantly. 

# generate animations

The main difference in this branch is that it generates an HTML page with animations of the objects of interest.

The following parameters should be considered when generating the page:
* obj_num_est:		the number of unique objects expected in the total images to be tracked (you may need to do a trial run to find this number)
* maxgifs: 			the total number of animations to be created for the HTML page
* padxy: 			each animation will be focused on the object of interest. Add some padding to the edges (in units of the (x,y) grid, not in number of pixels)
* padt: 			number of images to be shown before object first appears and after object disappears
* poi_string		poi = "property of interest"; poi_string = String of property to be used for selection of objects of interest CURRENT OPTIONS are 'area', 'extreme', and 'meanvar' NEW OPTIONS require these to be added to the function write_html_files in object_tracking.
* poi_increasing:	if True, looking for maxima
* HTML_DIR:			directory to store the HTML page, animations, and images

