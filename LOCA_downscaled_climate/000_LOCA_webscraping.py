#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 11/08/16


#---------------#
#--DESCRIPTION--#
#---------------#
# This script downlaods LOCA climate projections for Tmax/Tmin/ppt using the
# OPeNDAP interface. The code subsets by NM for all GCMs for three time periods:
# (1) historical:   1986-2005
# (2) mid-century:  2040-2059
# (3) late-century: 2080-2099


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
import datetime as dt
from ftplib import FTP
import glob
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
import re
from StringIO import StringIO
import urllib2
import urllib
from urllib import urlopen
import wget
from zipfile import ZipFile


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### find closest index to specified value within an array
def near(array, value):
	idx = (abs(array - value)).argmin()
	return idx

#### Create netcdfs from spatial/temporal subset
## @param var:   GCM being downloaded
## @param llon:  lower longitude
## @param ulon:  upper longitude
## @param llat:  lower latitude
## @param ulat:  upper latitude
## @param start: start index
## @param end:   end index
def down_to_ncdf(var1, llon, ulon, llat, ulat, start, end):
	try:
		# subset data by variable, then spatially and temporally
		var = nc.variables[var1]
		nm_only = var[start:end, llat:ulat, llon:ulon]
		# create ncdf for NM
		nc1 = Dataset(var1 + ".nc", 'w', format='NETCDF4_CLASSIC')
		nc1.createDimension('lon',  len(range(llon, ulon)))
		nc1.createDimension('lat',  len(range(llat, ulat)))
		nc1.createDimension('time', len(range(start,end)))
		longitudes = nc1.createVariable('longitude', np.float32, ('lon',))
		latitudes  = nc1.createVariable('latitude',  np.float32, ('lat',))
		time       = nc1.createVariable('time',      np.float32, ('time',))
		data       = nc1.createVariable('data',      np.float32, ('time', 'lat', 'lon'))
		nc1.variables['latitude'][:]  = lat[llat:ulat]
		nc1.variables['longitude'][:] = lon[llon:ulon]
		nc1.variables['data'][:] = nm_only
	finally:
		# clear from memory	
		nc1.close()


#---------------------------#
#-- SET WORKING DIRECTORY --#
#---------------------------#
path1 = "L:/LOCA_data_NM"
os.chdir(path1)	#set directory


#---------------------------------#
#-- DOWNLOAD/READ/SUBSET TO NM ---#
#---------------------------------#
#nc       = netCDF4.Dataset('http://cida.usgs.gov/thredds/dodsC/loca_historical')   #historical results
nc       = netCDF4.Dataset('http://cida.usgs.gov/thredds/dodsC/loca_future')   #future results
vars1    = nc.variables.keys()[6:176]  #future:6:176; historical: 6:78
lat      = nc.variables['lat'][:]
lon      = nc.variables['lon'][:]
time_var = nc.variables['time']
dtime    = netCDF4.num2date(time_var[:], time_var.units)

## define time period (#mid-century: 12418:19722; late-century: 27028:34332; historical: 13149:20453)
start = 12418
end   = 19722

## define bounding box (based on state boundaries)
llat = near(lat, 31.3322)
ulat = near(lat, 37.0003)
llon = near(lon, -109.05 + 360)
ulon = near(lon, -103.002 + 360)

## download subsetted data to netcdf
for var in vars1:
	down_to_ncdf(var1 = var, 
		         llon = llon, 
		         ulon = ulon, 
		         llat = llat, 
		         ulat = ulat, 
		         start = start, 
		         end = end)


#-------------------------------------#
#-- DOWNLOAD DIRECTLY FROM FTP SITE --#
#-------------------------------------#
## generate URLs
base_url = "ftp://gdo-dcp.ucllnl.org/pub/dcp/archive/cmip5/loca/LOCA_2016-04-02/"
ftp      = FTP('gdo-dcp.ucllnl.org')
ftp.login()
ftp.cwd('/pub/dcp/archive/cmip5/loca/LOCA_2016-04-02/')
GCMs     = ftp.nlst()
GCMs.remove("tmp")   #not permitted to access
res      = "16th"
RCPs     = ["rcp45","rcp85"]
r1i1p1   = "r1i1p1"
params   = 'tasmax'   #tasmax, tasmin, pr
dates    = [str(i) + "0101-" + str(i) + "1231" for i in range(2040,2060)]

url_list = []
for i in range(len(GCMs)):
	url1 = base_url + GCMs[i] + "/" + res
	for j in range(len(RCPs)):
		url2 = url1 + "/" + RCPs[j] + "/" + r1i1p1
		for k in range(len(params)):
			url3 = url2 + "/" + params[k] + "/" + params[k] + "_day_" + GCMs[i] + "_" + RCPs[j] + "_r1i1p1_"
			for l in range(len(dates)):
				url_list.append(url3 + dates[l] + ".LOCA_2016-04-02.16th.nc")


## download and subset to NM
for i in range(len(url_list)):
	try:
	    # download file
	    wget.download(url_list[i])
	    # open netcdf
	    nc_file = url_list[i].rsplit('/', 1)[-1]
	    fh   = Dataset(nc_file, mode='r')
	    lons = fh.variables['lon'][:]
	    lats = fh.variables['lat'][:]
	    data = fh.variables[params][:]
	    data_units = fh.variables[params].units
	    fh.close()
	    # subset ncdf by NM bounding box
	    lons = lons[271:367]
	    lats = lats[127:218]
	    data = data[ : , 127:218, 271:367]
	    # create new netcdf from subset
	    nc1 = Dataset(url_list[i].rsplit('/', 1)[-1], 'w', format='NETCDF4_CLASSIC')
	    nc1.createDimension('lon',  len(lons))
	    nc1.createDimension('lat',  len(lats))
	    nc1.createDimension('time', len(data))
	    longitudes = nc1.createVariable('longitude', np.float32, ('lon',))
	    latitudes  = nc1.createVariable('latitude', np.float32, ('lat',))
	    time       = nc1.createVariable('time', np.float32, ('time',))
	    data       = nc1.createVariable('data', np.float32, ('time', 'lat', 'lon'))
	    nc1.variables['latitude'][:] = lats
	    nc1.variables['longitude'][:] = lons
	    nc1.variables['data'][:] = data
	    # write file
	    nc1.close()
	except:
		pass


#------------------#
#--  MAKE PLOTS  --#
#------------------#
nc_file = 'L:/LOCA_data_NM/tasmax_CCSM4_r6i1p1_rcp45.nc'
fh = Dataset(nc_file, mode = 'r')
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
data = fh.variables['data'][:]
fh.close()
data = data[180, :, :]

lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(width = 5000000, height = 3500000,
            resolution = 'l', projection = 'stere',
            lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

# Plot Data
cs = m.pcolor(xi, yi, np.squeeze(data))
m.drawparallels(np.arange(-80., 81., 10.), labels = [1, 0, 0, 0], fontsize = 10)
m.drawmeridians(np.arange(-180., 181., 10.), labels = [0, 0, 0, 1], fontsize = 10)
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location = 'bottom', pad = "10%")

# Add Title
plt.title('PPT')

plt.show()


#---------------------------------------#
#-- TIME SERIES FOR A SINGLE LOCATION --#
#---------------------------------------#
nc = netCDF4.Dataset('http://cida.usgs.gov/thredds/dodsC/loca_future')

# specify lats/lons
lon1 = -105.9378 + 360
lat1 = 35.6870

# Find nearest index to desired location
ix = near(lon, lon1)
iy = near(lat, lat1)

# Specify time period
start = dt.datetime(1986,  1,  1, 12, 0)
stop  = dt.datetime(2005, 12, 31, 12, 0)

istart = netCDF4.date2index(start, time_var, select = 'nearest')
istop  = netCDF4.date2index(stop,  time_var, select = 'nearest')

# get time series for a given variable
var = nc.variables[vars1[0]]
hs  = var[istart:istop, iy, ix]
time = dtime[istart:istop]
x = np.vstack(list_a)
y = pd.DataFrame(x).transpose()
y = y.set_index(time)

# time series plot
ts.plot(title = 'Location: Lon=%.2f, Lat=%.2f' % ( lon[ix], lat[iy]), legend = True)
plt.ylabel(var.units)
plt.show()
