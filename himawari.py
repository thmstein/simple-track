#!/usr/bin/env python3

"""Download Himawari data from ICARE at the University of Lille for
every 20 minutes and convert from HDF to NetCDF, transforming to a
Cartesian grid and extracting a chosen domain.

Himawari documentation
======================
https://www.data.jma.go.jp/mscweb/en/himawari89/space_segment/spsg_ahi.html
http://www.icare.univ-lille1.fr/extract/subset/order

The following acknowledgement should be made when using Hiumawari data:
"Data provided by the Japan Meteorological Agency. We thank the ICARE
Data and Services Center for providing access to the data used in this
study."

Usage
=====
To see this message and all command line arguments:
$ himawari.py -h

Example: To download channel 13 (10.4 microns, recommended for tracking
MCSs) from 1 January to 3 January 2021 over Java:

$ himawari.py -e -12 102 -4 118 -d java 202101010000 202101032340 13

(Note, however, that 13 is the default channel so this may be omitted.)

By default, the full disc of data will be deleted to save storage space.
If you want to extract multiple domains for the same times, use the -k
(or --keep-full-disc) option. Next time the command is run, the domain
will be extracted from the full-disc data already saved, unless the -O
(or --overwrite-full-disc) option is used.

Python environment
==================
The following is a recommended way of setting up a conda environment to
run this code:

$ conda create --name himawari python=3.9
$ conda activate himawari
$ conda install -c conda-forge python-magic
$ conda install -c conda-forge gdal
$ conda install -c conda-forge netcdf4

If you see an error relating to the libpoppler.so shared library, you
may need to downgrade the poppler package (which is installed as a
dependency of gdal) manually, e.g.:
$ conda install poppler=0.81.0


Author: Simon Peatman (www.simonpeatman.me.uk)
"""

# Imports
import argparse
from pathlib import Path
import warnings
import hashlib
import datetime
import signal
import subprocess
import ftplib
import magic
import numpy as np
from osgeo import gdal
import netCDF4

# By default, GDAL just prints any error messages to sys.stdout; force
# it to raise errors as Python Exceptions instead
gdal.UseExceptions()

# Default location to save downloaded data
DDIR = Path.home() / 'himawari'

# File containing ICARE login information (should contain user name and
# password only, on separate lines)
#
# WARNING: Use a unique password for ICARE and restrict read permissions
# on your .icare file as appropriate!
ICARE_FILE = Path.home() / '.icare'

# Start of time series on ICARE
START_TSERIES = datetime.datetime(2015, 6, 1, 0, 0)

# Timeout limits in seconds
FTP_TIMEOUT = 60    # If FTP request fails, try HTTP instead
HTTP_TIMEOUT = 60   # If HTTP request fails, give up

# FTP host and path to data
FTP_HOST = 'ftp.icare.univ-lille1.fr'
FTP_DIR = 'SPACEBORNE/GEO/HIMAWARI+1407/GEO_L1B.v1.06'

# URL for HTTP download
#
# We really want relx=rely=0 and relw=relh=1 but the values we use here
# (0.000000001 and 0.999999999) are close enough to make no difference.
# We have to use these values to trick ICARE into thinking we are
# subsetting the data, because downloading via HTTP without subsetting
# is forbidden - users are told to use FTP instead.  However, the only
# reason we want to use HTTP in the first place is because some files
# are missing from the FTP server, including everything before 5th
# December 2017!
HTTP_URL = ('http://www.icare.univ-lille1.fr/extract/subset/get_data.php?relx='
            '0.000000001&rely=0.000000001&relw=0.999999999&relh=0.999999999&'
            'file=GEO/HIMAWARI+1407/GEO_L1B.v1.06/%Y/%Y_%m_%d/GEO_L1B-HIMA08_'
            '%Y-%m-%dT%H-%M-%S_G_{}_V1-06.hdf&Xdim=NbColumns&Ydim=NbLines')

# Wavelengths in microns for each channel
WAVELENGTHS = {1: 0.47, 2: 0.51, 3: 0.64, 4: 0.86, 5: 1.6, 6: 2.3, 7: 3.9,
               8: 6.2, 9: 6.9, 10: 7.3, 11: 8.6, 12: 9.6, 13: 10.4, 14: 11.2,
               15: 12.4, 16: 13.3}

# Channel codes in URL for each channel
CHANNEL_CODES = {4: 'VIS08', 5: 'IR016', 6: 'IR022', 7: 'IR038', 8: 'WV062',
                 9: 'WV069', 10: 'WV073', 11: 'IR085', 13: 'IR104',
                 14: 'IR112', 15: 'IR123', 16: 'IR132'}

# Time interval between files, and times (hours, minutes) to skip
TIME_INT = datetime.timedelta(minutes=20)
SKIP_TIMES = [(2, 40), (14, 40)]  # Housekeeping at 02:40 and 14:40 UTC

# File type we want to download
CORRECT_FILE_TYPE = 'Hierarchical Data Format (version 4) data'

# Constants
# Options for RESAMPLING: Average, Bilinear, Cubic, CubicSpline,
#                         Lanczos, Max, Med, Min, Mode,
#                         NearestNeighbour, Q1, Q3
# Options for ORDER: 1, 2, 3
EARTH_RADIUS_AT_EQUATOR = 6378137.0  # metres
EARTH_RADIUS_AT_POLES = 6356752.3    # metres
DEFAULT_RESOLUTION = 2000            # metres
FILL_VALUE = 65535
RESAMPLING = 'Bilinear'
ORDER = 3
LON0_DISC = 59.5147689316962
LON1_DISC = 221.885224964788
LAT0_DISC = -80.8302170404957
LAT1_DISC = 80.8302170404957

# Metadata
LONG_NAME = 'Himawari brightness temperature ({wavelength} microns)'
STANDARD_NAME = 'brightness_temperature'


class PermissionWarning(Warning):
    """Inappropriate permissions"""


class FileFormatError(Exception):
    """Wrong file format"""


def _timeout_handler(signum, frame):
    """Raise Error on timing out"""
    raise TimeoutError


def _parse_cl_args():
    """Parse command line arguments"""

    # Function to parse date-time string
    parse_datetime = lambda s: datetime.datetime.strptime(s, '%Y%m%d%H%M')

    # Declare arguments
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-e', '--extent', nargs=4,
                        help='lat0 lon0 lat1 lon1 (in degrees)')
    parser.add_argument('-o', '--outdir', type=Path, help='output directory')
    parser.add_argument('-d', '--domain-name', help='output subdirectory')
    parser.add_argument('-O', '--overwrite-full-disc', action='store_true')
    parser.add_argument('-k', '--keep-full-disc', action='store_true',
                        help=('keep full disc of Himawari data as well as '
                              'chosen domain'))
    parser.add_argument('-x', '--xres', type=float,
                        help='x resolution for Cartesian grid (metres)')
    parser.add_argument('-y', '--yres', type=float,
                        help='y resolution for Cartesian grid (metres)')
    parser.add_argument('start_time', type=parse_datetime, help='YYYYMMDDHHmm')
    parser.add_argument('end_time', type=parse_datetime, help='YYYYMMDDHHmm')
    parser.add_argument('channel', type=int, nargs='?', default=13,
                        help='channel number')

    # Parse arguments
    return parser.parse_args()


def valid_args(time1, time2, channel, extent=None, outdir=None):
    """Return a valid set of arguments, or raise error if not possible
    to do so.
    """

    msg = None

    # Check times
    if not isinstance(time1, datetime.datetime):
        msg = 'time1 must be a datetime.datetime object'
    if not isinstance(time2, datetime.datetime):
        msg = 'time2 must be a datetime.datetime object'
    if time1 < START_TSERIES:
        msg = START_TSERIES.strftime('time1 cannot be before %Y/%m/%d %H:%M')
    if time2 < time1:
        msg = 'time2 cannot be before time1'

    # Check channel
    if channel not in CHANNEL_CODES:
        msg = (f'channel must be one of {list(CHANNEL_CODES.keys())}')
    if msg is not None:
        raise ValueError(msg)

    # Check extent
    if extent is not None:
        if not isinstance(extent, (tuple, list)):
            raise TypeError('extent must be a tuple or list')
        if len(extent) != 4:
            raise ValueError('extent must have four values')

    # Check outdir
    if outdir is None:
        outdir = DDIR
    else:
        outdir = Path(outdir)

    # Return values
    return time1, time2, channel, extent, outdir


def _permissions_num2str(owner, group, universal):
    """Compute permissions in numerical form (e.g., 664) to string form
    (e.g., rw-rw-r--).
    """

    def num2str(num):
        s = ''
        for permission, value in zip('rwx', (4, 2, 1)):
            if num >= value:
                num -= value
                s += permission
            else:
                s += '-'
        return s
    return ''.join([num2str(int(x)) for x in [owner, group, universal]])


def check_permissions(path, *max_permissions):
    """Issue a warning if the given file has less restrictive
    permissions than those given (as integers 0-7, using the octal
    representation).
    """

    # Get permissions
    permissions = oct(Path(path).stat().st_mode)[-3:]
    owner, group, universal = permissions
    max_owner, max_group, max_universal = max_permissions

    # Check permissions
    warn = False
    if int(owner) > int(max_owner):
        warn = True
    if int(group) > int(max_group):
        warn = True
    if int(universal) > int(max_universal):
        warn = True
    if warn:
        permissions_str = _permissions_num2str(*permissions)
        max_permissions_str = _permissions_num2str(*max_permissions)
        warnings.warn(
            f'Permissions for {str(path)} are {permissions_str}; should be no '
            f'higher than {max_permissions_str}', PermissionWarning)


def get_icare_login_information(icare_file):
    """Get login information (user name, password and md5 hash of
    password) from the given file.

    This file should contain the user name and password only, in that
    order, on separate lines.

    The ICARE file should have limited read/write permissions; will
    raise warning if necessary.
    """

    # Check file permissions
    check_permissions(icare_file, 6, 4, 0)

    # Get user name and password from file
    with open(icare_file, 'r') as f:
        username, passwd = [l.strip() for l in f]

    # Run md5sum on password
    md5passwd = hashlib.md5(passwd.encode()).hexdigest()

    # Return
    return username, passwd, md5passwd


def iterate_datetime(start, end):
    """Iterate for every datetime from *start* to *end* INCLUSIVE, with
    an interval of *TIME_INT*.

    Skip any times in *SKIP_TIMES* (list of 2-tuples giving hours and
    minutes).
    """

    while True:
        if start > end:
            break
        if (start.hour, start.minute) not in SKIP_TIMES:
            yield start
        start += TIME_INT


def delete_file(path, error_if_not_found=False):
    """Delete the given path if it exists"""

    try:
        Path(path).unlink()
    except FileNotFoundError as err:
        if error_if_not_found:
            raise err


def check_file_type(path, file_type):
    """Raise exception if the given file is not of the given type"""

    path_file_type = magic.from_file(str(path))
    if path_file_type != file_type:
        msg = f'{str(path)} is a {path_file_type}, not {file_type}'
        raise FileFormatError(msg)


def hdf_save_path(outdir, channel):
    """Returns the path to which HDF files will be saved when
    downloaded. Path will contain strftime format codes.
    """

    odir = Path(outdir) / f'full_disc/{channel:02d}/%Y/%m/%d'
    ofile = f'himawari_{WAVELENGTHS[channel]}_%Y%m%d_%H%M.hdf'
    return odir / ofile


def nc_full_disc_path(hdf_path):
    """Returns the path to which the full disc of data will be saved as
    a NetCDF file after converting the downloaded HDF file.
    """

    return Path(str(Path(hdf_path).parent / Path(hdf_path).stem) + '.nc')


def nc_domain_path(outdir, channel, domain_name=None, extent=None):
    """Returns the path to which the extracted domain will be saved as a
    NetCDF file.

    Must provide either *domain_name* or *extent*.  If *domain_name* is
    provided, *extent* will be ignored.

    Path will contain strftime format codes.
    """

    hdf_stem = hdf_save_path(outdir, channel).stem
    if domain_name is not None:
        subdir = domain_name
    else:
        lat = lambda y: str(abs(y)) + ['S', 'N'][y >= 0]
        lon = lambda x: str(abs(x)) + ['W', 'E'][x >= 0]
        subdir = '_'.join([f(e) for f, e in zip([lat, lon, lat, lon], extent)])
    return outdir / subdir / f'{channel:02d}/%Y/%m/%d' / (hdf_stem + '.nc')


def download_hdf(time1, time2, channel, outdir=None, overwrite=False):
    """Download HDF Himawari data from ICARE, via FTP where possible
    (faster) or via HTTP as a fallback (not all data are on the FTP
    server).

    *channel* is an integer from 1 to 16 inclusive and corresponds to
    the following wavelengths in microns:

    channel |  wavelength
    --------|------------
        1   |     0.47   <-- unavailable
        2   |     0.51   <-- unavailable
        3   |     0.64   <-- unavailable
        4   |     0.86
        5   |     1.6
        6   |     2.3
        7   |     3.9
        8   |     6.2
        9   |     6.9
       10   |     7.3
       11   |     8.6
       12   |     9.6    <-- unavailable
       13   |    10.4
       14   |    11.2
       15   |    12.4
       16   |    13.3

    *time1* and *time2* are `datetime.datetime` objects.  Data are
    downloaded from *time1* to *time2* inclusive.

    If *overwrite* is False, silently skip over any files which already
    exist.
    """

    # Check arguments
    time1, time2, channel, _, outdir = valid_args(time1, time2, channel,
                                                  outdir=outdir)

    # Get ICARE login information
    username, passwd, md5passwd = get_icare_login_information(ICARE_FILE)

    # Register timeout handler
    signal.signal(signal.SIGALRM, _timeout_handler)

    # Connect to FTP host and log in
    opath = hdf_save_path(outdir, channel)
    ncpath = nc_full_disc_path(opath)
    print('\nDOWNLOADING:')
    with ftplib.FTP(FTP_HOST) as ftp:
        ftp.login(user=username, passwd=passwd)
        ftp.cwd(FTP_DIR)

        # Iterate for each time
        failed_files = []
        successful_files = 0
        for tt in iterate_datetime(time1, time2):

            # Get output location
            opath_t = tt.strftime(str(opath))
            if (Path(opath_t).is_file() or \
                Path(tt.strftime(str(ncpath))).is_file()) and not overwrite:
                continue
            Path(opath_t).parent.mkdir(parents=True, exist_ok=True)
            dfile = tt.strftime('GEO_L1B-HIMA08_%Y-%m-%dT%H-%M-%S_G_'
                                f'{CHANNEL_CODES[channel]}_V1-06.hdf')

            # Try downloading via FTP
            changed_dir = False
            try:
                ftp.cwd(tt.strftime(f'%Y/%Y_%m_%d'))
                changed_dir = True
                signal.alarm(FTP_TIMEOUT)
                with open(opath_t, 'wb') as fout:
                    ftp.retrbinary(f'RETR {dfile}', fout.write)
                check_file_type(opath_t, CORRECT_FILE_TYPE)

            # If FTP fails
            except Exception:

                # Delete local file (it will be corrupt/incomplete)
                delete_file(opath_t)

                # Try downloading via HTTP instead
                cmd = tt.strftime('curl --retry 5 --silent --show-error --data'
                                  f' "user={username}&passwd={md5passwd}" '
                                  f'"{HTTP_URL}" > {opath_t}')
                cmd = cmd.format(CHANNEL_CODES[channel])
                try:
                    signal.alarm(HTTP_TIMEOUT)
                    subprocess.check_output(cmd, shell=True)
                    check_file_type(opath_t, CORRECT_FILE_TYPE)

                # If HTTP fails
                except Exception:
                    print(f'MISSING FROM SERVER: {opath_t}')

                    # Delete local file (it will be corrupt/incomplete)
                    delete_file(opath_t)

                    # Append to list of failed files
                    failed_files.append(opath_t)

                # Successful download via HTTP
                else:
                    print('(HTTP)', opath_t)
                    successful_files += 1

            # Sucessful download via FTP
            else:
                print('(FTP) ', opath_t)
                successful_files += 1

            # If the FTP commands got as far as changing directory, go
            # back up two levels ready for the next file
            finally:
                if changed_dir:
                    ftp.cwd('../..')

    # Cancel timeout alarm
    signal.alarm(0)


def himawari_hdf_to_cartesian_nc(hdf_path_t, xres_m=None, yres_m=None,
                                 overwrite=False):
    """Convert a Himawari HDF file (full disc) to NetCDF, reprojecting
    from the geostationary nadir view to a Cartesian grid.

    *xres_m* and *yres_m* are the grid box lengths in metres.  If
    unspecified, both will default to *DEFAULT_RESOLUTION*.

    Returns the Path to the NetCDF file created.
    """

    # Check whether we need to do this conversion
    nc_path_t = nc_full_disc_path(hdf_path_t)
    if nc_path_t.is_file() and not overwrite:
        return nc_path_t
    temp_path = str(nc_path_t)[:-2] + 'temp{n}.nc'

    # Open file, and read in size and metadata
    ds = gdal.Open(str(hdf_path_t))
    nx, ny = ds.RasterXSize, ds.RasterYSize
    md = ds.GetMetadata()

    # Get details of orbit
    satellite_lon = float(md['Sub_Satellite_Longitude'])
    altitude_above_earth_centre_km = int(md['Altitude'])
    altitude_above_ground_m = altitude_above_earth_centre_km * 1000 - \
                              EARTH_RADIUS_AT_EQUATOR

    # Compute coordinates of corner of image, relative to centre
    if xres_m is None:
        xres_m = DEFAULT_RESOLUTION
    if yres_m is None:
        yres_m = DEFAULT_RESOLUTION
    xcorner = xres_m * nx // 2
    ycorner = yres_m * ny // 2

    # Convert to NetCDF and add metadata
    srs1 = (f'+proj=geos +h={altitude_above_ground_m} '
            f'+a={EARTH_RADIUS_AT_EQUATOR} +b={EARTH_RADIUS_AT_POLES} '
            f'+lon_0={satellite_lon} +x_0=0.0 +y_0=0.0 +ellps=WGS84 '
            '+datum=WGS84 +units=m +no_defs')
    ullr = [-xcorner, ycorner, xcorner, -ycorner]
    opts1 = gdal.TranslateOptions(outputSRS=srs1, outputBounds=ullr,
                                  format='netCDF', unscale=True,
                                  outputType=gdal.GDT_Float32)
    temp_path1 = temp_path.format(n=1)
    gdal.Translate(temp_path1, ds, options=opts1)
    ds = None  # close HDF file

    # Get missing values
    nodata = []
    for att in ['_FillValue', 'Missing_Output']:
        try:
            nodata.append(int(md[att]))
        except KeyError:
            pass

    # Resolution
    lat_range = LAT1_DISC - LAT0_DISC
    lon_range = LON1_DISC - LON0_DISC
    resx = lon_range / nx
    resy = lat_range / ny

    # Warp (i.e., reproject) to regular lat/lon grid
    srs2 = '+proj=latlong +datum=WGS84'
    extent = [LON0_DISC, LAT0_DISC, LON1_DISC, LAT1_DISC]
    opts2 = gdal.WarpOptions(
        srcNodata=nodata, dstNodata=FILL_VALUE, srcSRS=srs1, dstSRS=srs2,
        xRes=resx, yRes=resy, outputBounds=extent, outputType=gdal.GDT_Float32,
        resampleAlg=getattr(gdal, f'GRA_{RESAMPLING}'), polynomialOrder=ORDER, 
        format='netCDF')
    temp_path2 = temp_path.format(n=2)
    gdal.Warp(temp_path2, temp_path1, options=opts2)

    # Warp destroys metadata so add them back in
    for att in ['_FillValue', 'Missing_Output']:
        md[att] = FILL_VALUE
    md['Westernmost_Longitude'] = LON0_DISC
    md['Easternmost_Longitude'] = LON1_DISC
    md['Southernmost_Latitude'] = LAT0_DISC
    md['Northernmost_Latitude'] = LAT1_DISC
    metadata_pairs = [f'{k}={v}' for k, v in md.items()]
    opts3 = gdal.TranslateOptions(metadataOptions=metadata_pairs,
                                  format='netCDF')
    gdal.Translate(str(nc_path_t), temp_path2, options=opts3)

    # Delete HDF file and temporary NetCDF files
    del_paths = [hdf_path_t, str(hdf_path_t) + '.aux.xml']
    for n in range(1, 3):
        del_paths.append(temp_path.format(n=n))
    for del_path in del_paths:
        delete_file(del_path)
    return Path(nc_path_t)


def extract_nc_domain(source, target, extent):
    """Extract the given domain from the full disc of data.

    *source* and *target* = paths to input and output NetCDF files
    *extent* = lat0, lon0, lat1, lon1 (SW and NE corners)
    """

    # Read in NetCDF
    d = netCDF4.Dataset(source)
    Tb_arr = d.variables['Band1'][:]
    atts = {att.replace('GDAL_', ''): d.getncattr(att) for att in d.ncattrs()}
    d.close()

    # Work out domain we want to keep 
    lat_range = LAT1_DISC - LAT0_DISC
    lon_range = LON1_DISC - LON0_DISC
    ny, nx = Tb_arr.shape
    spacing_x = lon_range / (nx - 1)
    spacing_y = lat_range / (ny - 1)
    lats_all = np.arange(LAT0_DISC, LAT1_DISC + spacing_y / 2, spacing_y)
    lons_all = np.arange(LON0_DISC, LON1_DISC + spacing_x / 2, spacing_x)
    lat0, lon0, lat1, lon1 = [float(e) for e in extent]
    lats_keep = (lat0 <= lats_all) & (lats_all <= lat1)
    lons_keep = (lon0 <= lons_all) & (lons_all <= lon1)

    # Create variables
    Path(target).parent.mkdir(exist_ok=True, parents=True)
    root = netCDF4.Dataset(target, 'w')
    root.createDimension('time', 1)
    root.createDimension('latitude', lats_keep.sum())
    root.createDimension('longitude', lons_keep.sum())
    Tb_var = root.createVariable('T_b', 'f', ('time', 'latitude', 'longitude'),
                                 fill_value=Tb_arr.fill_value)
    time_var = root.createVariable('time', 'f', ('time',))
    lat_var = root.createVariable('latitude', 'f', ('latitude',))
    lon_var = root.createVariable('longitude', 'f', ('longitude',))

    # Set values and units
    Tb_var[:] = Tb_arr[None, lats_keep][..., lons_keep]
    channel_codes_dict_rev = {v: k for k, v in CHANNEL_CODES.items()}
    channel = channel_codes_dict_rev[atts['Band']]
    Tb_var.setncattr('long_name',
                     LONG_NAME.format(wavelength=WAVELENGTHS[channel]))
    Tb_var.setncattr('standard_name', STANDARD_NAME)
    Tb_var.setncattr('units', atts['units'])
    lat_var[:] = lats_all[lats_keep]
    lat_var.setncattr('standard_name', 'latitude')
    lat_var.setncattr('units', 'degrees_north')
    lon_var[:] = lons_all[lons_keep]
    lon_var.setncattr('standard_name', 'longitude')
    lon_var.setncattr('units', 'degrees_east')
    time_units = START_TSERIES.strftime('days since %Y-%m-%d %H:%M')
    nominal_time = datetime.datetime.strptime(atts['Nominal_Time'],
                                              '%Y-%m-%dT%H:%M:%SZ')
    seconds = (nominal_time - START_TSERIES).total_seconds()
    units = time_units.split()[0]
    seconds_per_unit = {'days': 86400, 'hours': 3600, 'minutes': 60,
                        'seconds': 1}[units]
    time_var[:] = [seconds / seconds_per_unit]
    time_var.setncattr('standard_name', 'time')
    time_var.setncattr('units', time_units)

    # Set global metadata
    for att in ['Altitude', 'Band', 'Band_Type', 'Beginning_Acquisition_Date',
                'Contact', 'Credits', 'End_Acquisition_Date', 'Icare_ID',
                'Input_Files', 'LibGeostat_Version', 'Missing_Output',
                'Nadir_Pixel_Size', 'GDAL', 'Nominal_Time',
                'Original_Nominal_Time', 'Production_Center',
                'Production_Date', 'Product_Description', 'Product_Name',
                'Product_Version', 'Projection_Longitude',
                'Radiance_Computation', 'Radiance_Computation_Info',
                'Reference', 'Scan_Origin', 'Sensors', 'Software_Version',
                'Sub_Satellite_Longitude', 'Sweep_Angle_Axis']:
        try:
            root.setncattr(att, atts[att])
        except KeyError:
            pass
    root.setncattr('channel', str(channel))
    root.setncattr('wavelength_in_microns', WAVELENGTHS[channel])

    # Save file
    root.close()
    print(target)
    

def check_file_list(time1, time2, channel, outdir=None, domain_name=None,
                    extent=None):
    """Print a list of missing files in the given time range"""

    # Check arguments
    time1, time2, channel, extent, outdir = valid_args(
        time1, time2, channel, extent=extent, outdir=outdir)

    # Iterate for each time 
    dpath = nc_domain_path(outdir, channel, domain_name=domain_name,
                           extent=extent)
    first = True
    for tt in iterate_datetime(time1, time2):

        # Check for file
        dpath_t = tt.strftime(str(dpath))
        if not Path(dpath_t).is_file():
            if first:
                print('\nMISSING FILES:')
                first = False
            print(dpath_t)


if __name__ == '__main__':

    # Parse and check command line arguments
    args = _parse_cl_args()
    time1, time2, channel, extent, outdir = valid_args(
        args.start_time, args.end_time, args.channel, extent=args.extent,
        outdir=args.outdir)

    # Download HDF files
    download_hdf(time1, time2, channel, outdir=outdir,
                 overwrite=args.overwrite_full_disc)

    # For each time, convert to NetCDF and extract domain
    hdf_path = hdf_save_path(outdir, channel)
    this_nc_domain_path = nc_domain_path(outdir, channel, extent=extent,
                                         domain_name=args.domain_name)
    print('\nCONVERTING TO NETCDF:')
    for tt in iterate_datetime(time1, time2):
        try:
            nc_path_t = himawari_hdf_to_cartesian_nc(
                tt.strftime(str(hdf_path)), xres_m=args.xres, yres_m=args.yres,
                overwrite=args.overwrite_full_disc)
        except RuntimeError:
            pass
        else:
            extract_nc_domain(nc_path_t, tt.strftime(str(this_nc_domain_path)),
                              extent)
            if not args.keep_full_disc:
                delete_file(nc_path_t)

    # Check for missing files
    check_file_list(time1, time2, channel, outdir=outdir,
                    domain_name=args.domain_name, extent=extent)
