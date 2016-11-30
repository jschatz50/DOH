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

# CDC claims to only calculate these metrics May-Sep, which seems stupid to me.  I'll calculate percentiles that way, 
# but not totals. If there's a 90F day during March, I'd like to know about it.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
import math
from StringIO import StringIO
import pandas as pd
import re
import glob
import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import Basemap


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
path1="L:/NLDAS"
os.chdir(path1)	#set directory

#import_list = glob.glob('L:/NLDAS/Annual_totals/*.nc')

infile = 'AT_daily_maxes'
nc_file = 'L:/NLDAS/' + infile + '.nc'
fh = Dataset(nc_file, mode = 'r')
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
dates = fh.variables['time'][:]
data = fh.variables['data'][:]
years = list(range(1981, 1981+len(data)/365))
fh.close()

# dates dataframe for indexing
df = pd.DataFrame(dates, columns=['DATE'])
df['DATE'] = df['DATE'].astype(str)
df['YEAR'] = df.DATE.str[0:4].astype(int)
df['MONTH'] = df.DATE.str[4:6].astype(int)
df['DAY'] = df.DATE.str[6:8].astype(int)

# calculate percentiles (1981-2010; May-Sept)
mask = ((df['YEAR'] >= 1981) & (df['YEAR'] <= 2010) & 
	    (df['MONTH'] >= 5) & (df['MONTH'] <= 9)).as_matrix()
subset = data[mask, :, :]
p90th = np.percentile(subset, q = 90, axis = 0)
p95th = np.percentile(subset, q = 95, axis = 0)
p98th = np.percentile(subset, q = 98, axis = 0)

# calculate days over thresholds for all years/days
thresholds = [p90th, p95th, p98th, F2C(90), F2C(95), F2C(100), F2C(105)]
names = ['p90th', 'p95th', 'p98th', 'd90F', 'd95F', 'd100F', 'd105F']

for i in range(len(thresholds)):
	threshold = thresholds[i]
	ntimes, ny, nx = np.shape(data)
	result = np.zeros((ntimes/365, ny, nx), dtype = int)
	for j in range(ntimes/365):
		ll = 0 + 365 * j    #upper limit day
		ul = 365 + 365 * j  #lower limit day
		for k in range(ll, ul):
			result[j, :, :] += data[k, :, :] >= threshold
	
	# generate netcdfs
	outname = infile + '_' + names[i] + '.nc'
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
	# write to file
	nc1.close()














#---------------#
#-- test plot --#
#---------------#
nc_file = 'L:/NLDAS/AT_daily_maxes_p90th.nc'
fh = Dataset(nc_file, mode = 'r')
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
data = fh.variables['data'][:]
dates = fh.variables['time'][:]
fh.close()
data = data[0, :, :]

# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(width = 650000, height = 750000,
            resolution = 'l', projection = 'stere',
            lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi, yi, np.squeeze(data))
m.drawcoastlines()
m.drawstates()
m.drawcountries()
cbar = m.colorbar(cs, location = 'bottom', pad = "10%")
plt.show()
