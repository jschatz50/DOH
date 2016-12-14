#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 11/14/16


#---------------#
#--DESCRIPTION--#
#---------------#
# This script calculates annual summaries of various
# heat-related metrics at each grid cell:
#   • days over 90, 95, 100, 105, 110°F
#   • absolute max
#   • longest stretch of days over 100°F
#   • nights over 80°F
#   • number consecutive nights over 80°F
#   • current vs future frequency of NWS heat advisories (assuming advisory criteria remain the same)
#   • current deviant temperature day frequency vs future (assuming criteria remain the same)


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
from pylab import *
import netCDF4
import fnmatch
import os
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np


#------------#
#-- SET WD --#
#------------#
path1 = "L:/LOCA_summaries_NM"
os.chdir(path1)	#set directory
os.getcwd()	#check directory


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
# find elements of a list matching a criteria
def find(lst, a):
    return [i for i, x in enumerate(lst) if x == a]


#-------------------#
#-- PROCESS NCDFs --#
#-------------------#
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


# For a list of temperature thresholds (90, 95, 100, 105, 110F), this loop calculates
# the number of days each year at or above those threholds.  It does this for each 
# emissions scenario within each (historical, RCP4.5, RCP8.5), each time period 
# (historical, mid-century, late-century), and each of the 28 GCMs within each
# period-scenario.  Then it combines that summary data into a single NCDF.  

# For example, the output "_____.nc" contains 28 layers that represent the mean
# annual number of 100F days for each of the 28 GCMs. Those 28 means were in turn 
# calculated from the 20 years of LOCA data from each GCM. The same is true for the
# average max, min, 10th percentile, and 90th percentile of annual 100F days.

# for each temp threshold, calculate annual summaries for each model, RCP, 
# and time period, & merge into NCDF
thresholds = [90, 95, 100, 105, 110]

for g in range(len(thresholds)):
	threshold_F = thresholds[g]  # 32.2=90°F; 35=95°F; 37.777=100°F; 40=105°F; 43=110°F
	threshold_C = (threshold_F - 32) * 5 / float(9)

	for h in range(len(list_a)):
		temp_list = list_a[h]

	    #create empty arrays for various traits
	    means = zeros((len(temp_list), 91, 96), dtype = int)
	    mins  = zeros((len(temp_list), 91, 96), dtype = int)
	    maxes = zeros((len(temp_list), 91, 96), dtype = int)
	    p90th = zeros((len(temp_list), 91, 96), dtype = int)
	    p10th = zeros((len(temp_list), 91, 96), dtype = int)

		for i in range(len(temp_list)):
			# open file; get variables
			nc1 = Dataset(temp_list[i], 'r')
			nc1.variables.keys()
			data = nc1.variables['data']

			# create empty array for all 20 annual totals
			ntimes, ny, nx = shape(data)
			result = zeros((20, ny, nx), dtype = int)
			for j in range(20):
				ll = 0 + 365 * j    #upper limit day
				ul = 365 + 365 * j  #lower limit day
				for k in range(ll, ul):
					result[j, :, :] += data[k, :, :] >= threshold_C
					
			# calculate summary traits of the 20 annual totals
			means[i, :, :] = np.mean(result, axis = 0)
			mins[i, :, :]  = np.min(result, axis = 0)
			maxes[i, :, :] = np.max(result, axis = 0)
			p90th[i, :, :] = np.percentile(result, q = 90, axis = 0) 
			p10th[i, :, :] = np.percentile(result, q = 10, axis = 0)

		# create ncdfs of results
		all_results = [means, mins, maxes, p90th, p10th]
		all_results_names = ["means", "mins", "maxes", "p90th", "p10th"]

		for l in range(len(all_results)):
			output = all_results[l]
			pieces = temp_list[0].rsplit('\\', )
			outname = 'd' + str(threshold_F) + 'F_' + pieces[1] + '_' + pieces[2] + '_' + all_results_names[l] + '.nc'
			nc2 = Dataset(outname, 'w', format = 'NETCDF4_CLASSIC')
			nc2.createDimension('lon',  len(nc1.variables['longitude']))
			nc2.createDimension('lat',  len(nc1.variables['latitude']))
			nc2.createDimension('time', len(output))
			longitudes = nc2.createVariable('longitude', np.float32, ('lon',))
			latitudes  = nc2.createVariable('latitude',  np.float32, ('lat',))
			time       = nc2.createVariable('time',      np.float32, ('time',))
			data       = nc2.createVariable('data',      np.float32, ('time','lat', 'lon'))
			nc2.variables['latitude'][:]  = nc1.variables['latitude'][:]
			nc2.variables['longitude'][:] = nc1.variables['longitude'][:]
			nc2.variables['data'][:] = output
			# write to file
			nc2.close()


#------------------#
#--  MAKE PLOTS  --#
#------------------#
f = Dataset("L:/LOCA_data_NM/results/d100F_historical_1986-2005_historical_means.nc", 'r')
f.variables.keys()
data = f.variables['data']
data = data[4, :, :]

lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
lon_0 = lons.mean()
lat_0 = lats.mean()
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

m = Basemap(width = 5000000, height = 3500000,
            resolution = 'l', projection = 'stere',
            lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)
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
