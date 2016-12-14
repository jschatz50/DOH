#------------#
#-- AUTHOR --#
#------------#
# Jason Schatz
# Created:  11/21/16
# Last modified:  12/05/16


#-----------------#
#-- DESCRIPTION --#
#-----------------#
# processes NLDAS2 forcing data for New Mexico. Script calculates
# daily min/max temperature, apparent temperature, and heat index 
# at 1/8th degree resolution. This is the data source used
# by the CDC, so we are using it for consistency.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
from datetime import datetime
from datetime import timedelta
import glob
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
from netCDF4 import Dataset
import numpy as np
from operator import itemgetter
import os
import pandas as pd
import re
from StringIO import StringIO
import time


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#

#### 1. Convert specific humidity to relative humidity
##' from Bolton 1980 The computation of Equivalent Potential Temperature 
##' \url{http://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html}
##' @param qair: specific humidity, dimensionless ratio of water mass / total air mass(e.g. kg/kg)
##' @param tair: degrees C
##' @param press: pressure in mb
##' @return rh: relative humidity
def qair2rh(qair, tair, altitudes):
	press = 10 * 101.3 * np.exp(-3.42e-2 * altitudes / (tair + 273.15))
	es = 6.112 * np.exp((17.67 * tair)/(tair + 243.5))
	e  = qair * press / (0.378 * qair + 0.622)
	rh = e / es
	rh[rh > 1] = 1
	rh[rh < 0] = 0
	return rh


#### 2. Calculate apparent temperature
##' AT = −1.3 + 0.92 * T + 2.2 * e (Steadman 1984), where
##' tair is dry bulb temperature (°C) and e is vapor pressure
##' (kPa), calculated after Buck (1981).
def ATcalc(tair, rh):
	es = 0.61078 * np.exp((17.269 * tair)/(tair + 237.3))
	e  = es * rh
	AT = -1.3 + 0.92 * tair + 2.2 * e
	return AT

tair = 22.52
rh = 0.171

#### 3. Calculate heat index
##' from National Weather Service 
##' \url{http://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml}
##' @param qair: specific humidity, dimensionless ratio of water mass / total air mass(e.g. kg/kg)
##' @param tair: degrees C
##' @param press: pressure in mb
##' @return rh: relative humidity

# 3.1. base equation
def HIbase(tair, rh):
	HI = -42.379 + \
	     2.04901523 * tair + \
	     10.14333127 * rh - \
	     0.22475541 * tair * rh - \
	     0.00683783 * tair * tair - \
	     0.05481717 * rh * rh + \
	     0.00122874 * tair * tair * rh + \
	     0.00085282 * tair * rh * rh - \
	     0.00000199 * tair * tair * rh * rh
	return HI

# 3.2. Set conditions for applying different versions of the base equation at different temps/RHs
def conditions(df):
	rh = df['rh']
	tair = df['tair']
	if rh < 13 and tair >= 80 and tair < 112:
		return HIbase(tair, rh) - ((13-rh)/4) * ((17 - abs(tair - 95))/17) ** 0.5
	elif rh > 85 and tair > 80 and tair < 87:
		return HIbase(tair, rh) + ((rh - 85) / 10) * ((87 - tair) / 5)
	elif tair < 80:
		return 0.5 * (tair + 61.0 + ((tair - 68.0) * 1.2) + (rh * 0.094))
	else:
		return HIbase(tair, rh)

# 3.3. final HI calculations
def HIcalc(tair, rh):
	tair = tair * 9/5 + 32   #convert to farenheit
	rh = rh * 100
	df = pd.DataFrame({'tair':tair.flatten()[:],'rh':rh.flatten()[:]})
	df['HI'] = df.apply(conditions, axis = 1)
	HI = np.reshape(df['HI'], tair.shape)   # reshape to match dimensions of array
	return HI


#### 4. Create netcdfs from daily mins/maxes
## @param output name (e.g. 'filename.nc')
## @param data_source (where the array data for input into the nc file comes from)
def makencdf(name, data_source):
	try:
		nc1 = Dataset(name, 'w', format ='NETCDF4_CLASSIC')
		nc1.createDimension('lon',  len( lon))
		nc1.createDimension('lat',  len(lat))
		nc1.createDimension('time', len(dates_list))
		longitudes = nc1.createVariable('longitude', np.float32, ('lon',))
		latitudes  = nc1.createVariable('latitude',  np.float32, ('lat',))
		time       = nc1.createVariable('time',      np.int, ('time',))
		data       = nc1.createVariable('data',      np.float32, ('time', 'lat', 'lon'))
		nc1.variables['time'][:]      = dates
		nc1.variables['latitude'][:]  = lat
		nc1.variables['longitude'][:] = lon
		nc1.variables['data'][:]      = data_source
	finally:
		# remove file from memory
		nc1.close()


#-------------------#
#-- PROCESS FILES --#
#-------------------#
# for each day:
# 1. create 24-hour-length array for both temp and sphum
# 2. compute rh from those arrays
# 3. compute heat index and apparent temperature
# 4. grab daily min/max for that day (for Tair, HI, and AT)
# 5. put those daily tmins/tmaxes into separate netcdfs, which will contain all measurements
# 6. write netcdfs to file

path1 = "L:/NLDAS/Daily"
os.chdir(path1)	#set directory

import_list = glob.glob('L:/NLDAS/raw_data/*.nc')
dates_list = pd.date_range(pd.to_datetime("19810101", format = '%Y%m%d'), periods = 13087).tolist()
dates = [str(dates_list[i].year) + str("%02d" % (dates_list[i].month)) + str("%02d" % (dates_list[i].day)) for i in range(len(dates_list))]
altitudes = np.loadtxt('L:/NLDAS/nldas_grid_altitudes.txt')

# lists for converting GMT to local time (-7 hrs to MST)
# this just grabs the index for the local date; then the index gets applied to the path list
local_dates_times  = [datetime.strptime(import_list[i][36:44] + ' ' + import_list[i][45:49], '%Y%m%d %H%M') \
                     - timedelta(hours = 7) for i in range(len(import_list))]
local_dates        = [str(local_dates_times[i].year) + str("%02d" % (local_dates_times[i].month)) + str("%02d" % (local_dates_times[i].day)) \
                      for i in range(len(local_dates_times))]

# create empty arrays for daily results
tair_mins  = np.zeros((len(dates_list), 49, 53), dtype = float)
tair_maxes = np.zeros((len(dates_list), 49, 53), dtype = float)
HI_mins    = np.zeros((len(dates_list), 49, 53), dtype = float)
HI_maxes   = np.zeros((len(dates_list), 49, 53), dtype = float)
AT_mins    = np.zeros((len(dates_list), 49, 53), dtype = float)
AT_maxes   = np.zeros((len(dates_list), 49, 53), dtype = float)

start = time.clock()
for i in range(len(dates_list)):
	try:
		year  = dates_list[i].year
		month = dates_list[i].month
		day_of_month = dates_list[i].day
		date  = str(year) + str("%02d" % (month)) + str("%02d" % (day_of_month))
		indices = [y for y, x in enumerate(local_dates) if x == date]
		paths = [itemgetter(y)(import_list) for y in indices]
#		paths = [s for s in import_list if date in s]

		#create empty arrays for hourly results
		tair_hourly = np.zeros((len(paths), 49, 53), dtype = float)
		HI_hourly = np.zeros((len(paths), 49, 53), dtype = float)
		AT_hourly = np.zeros((len(paths), 49, 53), dtype = float)

		for j in range(len(paths)):
			nc    = netCDF4.Dataset(paths[j])
			vars1 = nc.variables.keys()
			lat   = nc.variables[vars1[1]][:]
			lon   = nc.variables[vars1[0]][:]
			sphum = nc.variables[vars1[2]][:]
			tair  = nc.variables[vars1[3]][:] - 273.15   #K to C
			rh = qair2rh(sphum, tair, altitudes)
			AT = ATcalc (tair, rh)
			HI = HIcalc(tair, rh)
			# append arrays by hour
			tair_hourly[j, :, :] = tair
			HI_hourly[j, :, :]   = HI
			AT_hourly[j, :, :]   = AT

    	# summary traits by day
		tair_mins[i, :, :]  = np.min(tair_hourly, axis = 0)
		tair_maxes[i, :, :] = np.max(tair_hourly, axis = 0)
		HI_mins[i, :, :]    = np.min(HI_hourly, axis = 0)
		HI_maxes[i, :, :]   = np.max(HI_hourly, axis = 0)
		AT_mins[i, :, :]    = np.min(AT_hourly, axis = 0)
		AT_maxes[i, :, :]   = np.max(AT_hourly, axis = 0)

	except:
		pass
end = time.clock()
print(end - start)


#------------------------------------------#
#--      create NETCDFs from output      --#
#-- (daily timesteps for each parameter) --#
#------------------------------------------#
makencdf('Tair_daily_mins.nc', tair_mins)
makencdf('Tair_daily_maxes.nc', tair_maxes)
makencdf('HI_daily_mins.nc', HI_mins)
makencdf('HI_daily_maxes.nc', HI_maxes)
makencdf('AT_daily_mins.nc', AT_maxes)
makencdf('AT_daily_maxes.nc', AT_maxes)


#---------------#
#-- test plot --#
#---------------#
nc_file = 'L:/NLDAS/Tair_daily_maxes.nc'
fh = Dataset(nc_file, mode = 'r')
lons  = fh.variables['longitude'][:]
lats  = fh.variables['latitude'][:]
data  = fh.variables['data'][:]
dates = fh.variables['time'][:]
fh.close()
data = data[0, :, :]

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
