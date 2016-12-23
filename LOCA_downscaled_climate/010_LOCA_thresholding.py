#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# Created: 11/14/16
# Last modified: 12/22/16


#---------------#
#--DESCRIPTION--#
#---------------#
# This script calculates annual summaries of various
# heat-related metrics at each grid cell:
#   • days over 90, 95, 100, 105, 110°F
#   • longest stretch of days over those same thresholds each year
#   • absolute max
#   • nights over 80°F
#   • number consecutive nights over 80°F
#   • current vs future frequency of NWS heat advisories (assuming advisory criteria remain the same)
#   • current deviant temperature day frequency vs future (assuming criteria remain the same)


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
import fnmatch
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
from pylab import *


#------------#
#-- SET WD --#
#------------#
path1 = "L:/LOCA_summaries_NM"
os.chdir(path1)	#set directory
os.getcwd()	#check directory


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### create netcdfs from arrays of results
## @param output name (e.g. 'filename.nc')
## @param data_source (where the array data for input into the nc file comes from)
def makencdf(name, data_source, lat, lon):
	try:
		nc_new = Dataset(name, 'w', format ='NETCDF4_CLASSIC')
		nc_new.createDimension('lon',  len(lon))
		nc_new.createDimension('lat',  len(lat))
		nc_new.createDimension('time', len(data_source))
		longitudes = nc_new.createVariable('longitude', np.float32, ('lon',))
		latitudes  = nc_new.createVariable('latitude',  np.float32, ('lat',))
		time       = nc_new.createVariable('time',      np.int, ('time',))
		data       = nc_new.createVariable('data',      np.float32, ('time', 'lat', 'lon'))
		nc_new.variables['latitude'][:]  = lat
		nc_new.variables['longitude'][:] = lon
		nc_new.variables['data'][:]      = data_source
	finally:
		# remove file from memory
		nc_new.close()

#### find maximum number of consecutive days above a given temp threshold
## @param array1:  array of temperatures
## @param threshold:  temperature threhold
def max_consec(array1, threshold):
	array1 = array1 >= threshold
	## bound arrays by zeros
	bounded = np.insert(array1, [0, np.shape(array1)[0]], 0, axis=0)
	diffs = np.diff((bounded).astype(int), axis=0)
	return np.apply_along_axis(return_max_consec, 0, diffs)

#### sub-function of max_consec
## @param a:  input array (i.e. diffs, from max_consec function)
def return_max_consec(a):
	try:
		starts = np.argwhere(a == 1)
		stops  = np.argwhere(a == -1)
		maxrun = np.max(np.column_stack((stops[:,0] - starts[:,0])))
	except:
		maxrun = 0
		return maxrun
	else:
		return maxrun

#### farenheit to celsius
def F2C(TempF):
	return (TempF - 32) * 5 / float(9)


#---------------#
#-- PREP DATA --#
#---------------#
# create list of all tmax or tmin files
import_list = []
for root, dirnames, filenames in os.walk('L:/LOCA_data_NM'):
    for filename in fnmatch.filter(filenames, '*.nc'):
        import_list.append(os.path.join(root, filename))

# separate the file paths into categories so they can be merged into period/RCP specific NCDFs
tmax_list   = [i for i in import_list if "tmax" in i]
hist_list   = [i for i in tmax_list if "historical_1986-2005" in i]
mid45_list  = [i for i in tmax_list if "mid-century_2040-2059" in i][:28]
mid85_list  = [i for i in tmax_list if "mid-century_2040-2059" in i][28:]
late45_list = [i for i in tmax_list if "late-century_2080-2099" in i][:28]
late85_list = [i for i in tmax_list if "late-century_2080-2099" in i][28:]
list_a = [hist_list, mid45_list, mid85_list, late45_list, late85_list]


#------------------#
#-- thresholding --#
#------------------#
# For a list of temperature thresholds, this loop calculates the number of days 
# each year at or above those threholds, as well as the longest stretch of days 
# per year above that threshold for each cell.  It does this for each emissions 
# scenario (historical, RCP4.5, RCP8.5), each time period (historical, 
# mid-century, late-century), and each of the 28 GCMs within each period-scenario.  
# Then it combines that summary data into a single NCDF. For example, the output 
# "_____.nc" contains 28 layers that represent the mean annual number of 100F 
# days for each of the 28 GCMs. Those 28 means were in turn calculated from the 
# 20 years of LOCA data from each GCM. The same is true for the average max, min, 
# 10th percentile, and 90th percentile of annual 100F days.

## specify whether you want annual totals or max consecutive days?
agg_type = 'total'   # enter 'consecutive' or 'total'

## define temperature thresholds
thresholds = [90, 95, 100, 105, 110]

## for each temperature threshold:
for thresh in thresholds:

    ## for each time period/emissions scenario:
	for lst in list_a:

	    ## create empty arrays for various traits
	    means = zeros((len(lst), 91, 96), dtype = int)
	    mins  = zeros((len(lst), 91, 96), dtype = int)
	    maxes = zeros((len(lst), 91, 96), dtype = int)
	    p90th = zeros((len(lst), 91, 96), dtype = int)
	    p10th = zeros((len(lst), 91, 96), dtype = int)

	    ## for each file path in lst
		for i, path in enumerate(lst):

			## open netcdf
			nc1 = Dataset(path, 'r')
			nc1.variables.keys()
			data = nc1.variables['data']

			## create empty array for the 20 annual totals
			ntimes, ny, nx = shape(data)
			result = zeros((20, ny, nx), dtype = int)

			## for each of the 20 years
			for j in range(20):
				ll = 0 + 365 * j    #upper limit day
				ul = 365 + 365 * j  #lower limit day

				## calculate most consecutive days over threshold or annual total over threshold
				if agg_type == 'consecutive':
					result[j, :, :] = max_consec(data[ll:ul, :, :], F2C(thresh))
				else:
					for k in range(ll, ul):
						result[j, :, :] += data[k, :, :] >= F2C(thresh)
					
			## calculate summary traits of the 20 annual values
			means[i, :, :] = np.mean(result, axis = 0)
			mins[i, :, :]  = np.min(result, axis = 0)
			maxes[i, :, :] = np.max(result, axis = 0)
			p90th[i, :, :] = np.percentile(result, q = 90, axis = 0) 
			p10th[i, :, :] = np.percentile(result, q = 10, axis = 0)

		## create netcdfs of results
		all_results = [means, mins, maxes, p90th, p10th]
		all_results_names = ["means", "mins", "maxes", "p90th", "p10th"]
		for l in range(len(all_results)):
			output = all_results[l]
			pieces = lst[0].rsplit('\\', )
			outname = 'd' + str(thresh) + 'F_' + pieces[1] + '_' + pieces[2] + '_' + all_results_names[l] + '.nc'
			makencdf(name = outname, 
				     data_source = output,
				     lat = nc1.variables['latitude'][:],
				     lon = nc1.variables['longitude'][:])


#------------------#
#--  MAKE PLOTS  --#
#------------------#
f = Dataset("L:/LOCA_summaries_NM/d90F_historical_1986-2005_consec_means.nc", 'r')
f.variables.keys()
data = f.variables['data']
data = data[4, :, :]

lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
lon_0 = lons.mean()
lat_0 = lats.mean()
lon, lat = np.meshgrid(lons, lats)

m = Basemap(width = 5000000, height = 3500000,
            resolution = 'l', projection = 'stere',
            lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)
xi, yi = m(lon, lat)
cs = m.pcolor(xi, yi, np.squeeze(data))
m.drawparallels(np.arange(-80., 81., 10.), labels = [1,0,0,0], fontsize = 10)
m.drawmeridians(np.arange(-180., 181., 10.), labels = [0,0,0,1], fontsize = 10)
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# colorbar
cbar = m.colorbar(cs, location = 'bottom', pad = "10%")

# title
plt.title('Days over 100F per year')
plt.show()
