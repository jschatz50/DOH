#------------#
#-- AUTHOR --#
#------------#
# Jason Schatz
# Created:  11/29/16
# Last modified:  11/29/16


#-----------------#
#-- DESCRIPTION --#
#-----------------#
# Processes NLDAS2 forcing data for New Mexico. Script calculates
# summer percentiles for daily daily min/max temperature, apparent 
# temperature, and heat index for 1981-2010 for each grid cell.
# The script then applies those percentiles to calculate the number
# of days in each year that temperatures exceed those percentiles.
# Also calculates number of days above absolute thresholds.
# Thresholds used by the CDC are:
# -days over 90, 95, 100, 105F
# -days over the 90th, 95th, and 98th percentiles
# https://ephtracking.cdc.gov/showIndicatorPages.action?selectedContentAreaAbbreviation=15&selectedIndicatorId=79

# CDC claims to only calculate these metrics May-Sep.  I'll calculate percentiles that way, 
# but not totals. If there's a 90F day during March, I'd like to know about it.  Maybe they
# do, too, but it's not clear on their documentation.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
import datetime as dt
import glob
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
import re
from StringIO import StringIO


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### Farenheit to Celsius
def F2C(tempF):
	tempC = (tempF - 32) * 5/float(9)
	return tempC


#-------------------#
#-- PROCESS FILES --#
#-------------------#
path1="L:/NLDAS/Annual_summaries"
os.chdir(path1)	#set directory

import_list = glob.glob('L:/NLDAS/Daily/*.nc')

filenum = 2   # only need a few for thresholding
nc_file = import_list[filenum]
fh = Dataset(nc_file, mode = 'r')
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
dates = fh.variables['time'][:]
data1 = fh.variables['data'][:]
years = list(range(1981, 1981 + len(data1) / 365 + 1))
fh.close()

# dates dataframe for indexing
df = pd.DataFrame(dates, columns = ['DATE'])
df['DATE'] = df['DATE'].astype(str)
df['YEAR'] = df.DATE.str[0:4].astype(int)
df['MONTH'] = df.DATE.str[4:6].astype(int)
df['DAY'] = df.DATE.str[6:8].astype(int)

# calculate percentiles (1981-2010; May-Sept)
mask = ((df['YEAR'] >= 1981) & (df['YEAR'] <= 2010) & 
	    (df['MONTH'] >= 5) & (df['MONTH'] <= 9)).as_matrix()
subset = data1[mask, :, :]
p90th = np.percentile(subset, q = 90, axis = 0)
p95th = np.percentile(subset, q = 95, axis = 0)
p98th = np.percentile(subset, q = 98, axis = 0)

# calculate days over thresholds for all years/days
#thresholds = [p90th, p95th, p98th, F2C(90), F2C(95), F2C(100), F2C(105)]   #for AT and Tair
thresholds = [p90th, p95th, p98th, 90, 95, 100, 105]    #for HI
names = ['p90th', 'p95th', 'p98th', 'd90F', 'd95F', 'd100F', 'd105F']

for i in range(len(thresholds)):
	try:
		threshold = thresholds[i]
		ntimes, ny, nx = np.shape(data1)
		result = np.zeros((ntimes / 365 + 1, ny, nx), dtype = int)
		for j in range(ntimes / 365 + 1):
			try:
				ll = 0 + 365 * j    #upper limit day
				ul = 365 + 365 * j  #lower limit day
				for k in range(ll, ul):
					result[j, :, :] += data1[k, :, :] >= threshold
			except:
				pass
		# generate netcdfs
		outname = import_list[filenum][15:] + '_' + names[i] + '.nc'
		nc1 = Dataset(outname, 'w', format='NETCDF4_CLASSIC')
		nc1.createDimension('lon',  len(lons))
		nc1.createDimension('lat',  len(lats))
		nc1.createDimension('time', len(years))
		longitudes = nc1.createVariable('longitude', np.float32, ('lon',))
		latitudes  = nc1.createVariable('latitude',  np.float32, ('lat',))
		time       = nc1.createVariable('time',      np.float32, ('time',))
		data       = nc1.createVariable('data',      np.float32, ('time', 'lat', 'lon'))
		nc1.variables['latitude'][:]  = lats
		nc1.variables['longitude'][:] = lons
		nc1.variables['time'][:] = years
		nc1.variables['data'][:] = result
	finally:
		nc1.close()


#---------------#
#-- test plot --#
#---------------#
nc_file = 'L:/NLDAS/Annual_summaries/AT_daily_maxes.nc_d90F.nc'
fh = Dataset(nc_file, mode = 'r')
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
data = fh.variables['data'][:]
dates = fh.variables['time'][:]
fh.close()
data = data[31, :, :]

lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(width = 650000, height = 750000,
            resolution = 'l', projection = 'stere',
            lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi, yi, np.squeeze(data))
m.readshapefile('H:/Jason/GIS/Political_boundaries/NM_counties/NM_counties', 'NM_counties', color = 'w')
m.readshapefile('H:/Jason/GIS/Political_boundaries/NM_state_boundary/NM_state', 'NM_state', color = 'black', linewidth = 3)
cbar = m.colorbar(cs, location = 'bottom', pad = "10%")
plt.show()
