#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# Created:       12/13/16
# Last modified: 12/13/16 


#---------------#
#--DESCRIPTION--#
#---------------#
# This script samples LOCA downscaled climate projections in order to provide
# climatic covariates for annual fire acreage model.  Specifically, it 
# calculates the spring/summer mean tmax or total ppt for each year, averages the 
# results of all models for each year (i.e. ensemble mean), then samples those
# averages using the coordinates of the stations used to build the fire model.  
# Final file contains: year and median spring/summer tmax and ppt (median across the
# station locations used to build the fire model).


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
import netCDF4
from netCDF4 import Dataset
import fnmatch
import glob
import numpy as np
import os
import pandas as pd


#------------#
#-- SET WD --#
#------------#
path1 = "H:/Jason/Wildfire/future_climate_data_for_model"
os.chdir(path1)	#set directory


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### Create netcdfs from results array
## @param outname (e.g. 'filename.nc')
## @param data_source (where the array data for input into the nc file comes from)
def make_ncdf(outname, data_source):
	try:
		new_nc = Dataset(outname, 'w', format = 'NETCDF4_CLASSIC')
		new_nc.createDimension('lon',  len(nc1.variables['longitude']))
		new_nc.createDimension('lat',  len(nc1.variables['latitude']))
		new_nc.createDimension('time', len(data_source))
		longitudes = new_nc.createVariable('longitude', np.float32, ('lon',))
		latitudes  = new_nc.createVariable('latitude',  np.float32, ('lat',))
		time       = new_nc.createVariable('time',      np.int,     ('time',))
		data       = new_nc.createVariable('data',      np.float32, ('time', 'lat', 'lon'))
		new_nc.variables['latitude'][:]  = nc1.variables['latitude'][:]
		new_nc.variables['longitude'][:] = nc1.variables['longitude'][:]
		new_nc.variables['data'][:]      = data_source
	finally:
		# remove file from memory
		new_nc.close()

#### finds index of nearest value in array
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx


#-------------------#
#-- PROCESS NCDFs --#
#-------------------#
# This loop calculates the spring/summer mean tmax or total ppt in each year.  
# It does this for each emissions scenario (RCP4.5, RCP8.5) within each time 
# period (historical, mid-century, late-century), for each of the 28 GCMs 
# within each period-scenario.

## create list of all netcdf files
import_list = []
for root, dirnames, filenames in os.walk('L:/LOCA_data_NM'):
    for filename in fnmatch.filter(filenames, '*.nc'):
        import_list.append(os.path.join(root, filename))

## separate the file paths into categories so they can be merged into period/RCP specific NCDFs
parameter = 'tmax'   # 'tmax' or 'precipitation'

master_list = [i for i in import_list if parameter in i]
hist_list   = [i for i in master_list if "historical_1986-2005" in i]
mid45_list  = [i for i in master_list if "mid-century_2040-2059" in i][:28]
mid85_list  = [i for i in master_list if "mid-century_2040-2059" in i][28:]
late45_list = [i for i in master_list if "late-century_2080-2099" in i][:28]
late85_list = [i for i in master_list if "late-century_2080-2099" in i][28:]
list_a = [hist_list, mid45_list, mid85_list, late45_list, late85_list]

## for each time period...
for lst in list_a:
	
	## create empty array for results from each year
	means = np.zeros((20, 91, 96), dtype = float)

    ## define date brackets for slicing by season of year
	for j in range(20):
		ll = 60 + 365 * j    #upper limit day (start of spring)
		ul = 244 + 365 * j   #lower limit day (end of summer)

        ## create empty array for a given year's results for all models
		annual = np.zeros((len(lst), 91, 96), dtype = float)

		## for each model...
		for index, path in enumerate(lst):
			
			## open file, slice by date range, take mean, append to 'annual' array
			nc1 = Dataset(path, 'r')
			data = nc1.variables['data'][ll:ul, :, :]

			if parameter == 'tmax':
				annual[index, :, :] = np.mean(data, axis = 0)   # mean if tmax
			else:
				annual[index, :, :] = np.sum(data, axis = 0)   # sum if precipitation

		## take annual mean of all models in 'annual' array; append to 'means'
		means[j, :, :] = np.mean(annual, axis = 0)

	## create netcdfs of results
	pieces = lst[0].rsplit('\\', )
	outname = pieces[1] + '_spring_summer_' + pieces[2] + '_' + pieces[3] + '.nc'
	make_ncdf(outname = outname, 
		     data_source = means)


#---------------------------------------#
#-- CALCULATE ANNUAL MODEL COVARIATES --#
#---------------------------------------#
# (1) samples netcdfs using the coordinates of the stations used to build the fire model;
# (2) takes statewide median;
# (3) writes results to file
coords = pd.read_csv('H:/Jason/Wildfire/future_climate_data_for_model/stations_used_for_model.csv')
import_list = glob.glob('H:/Jason/Wildfire/future_climate_data_for_model/*.nc')

for path in import_list:
	nc1 = Dataset(path, 'r')
	data = nc1.variables['data'][:, :, :]
	lats = nc1.variables['latitude'][:]
	lons = nc1.variables['longitude'][:]

	Xs = coords['LONGITUDE'] + 360
	Ys = coords['LATITUDE']
	Ys = [find_nearest(array = lats, value = Ys[i]) for i in range(len(Ys))]
	Xs = [find_nearest(array = lons, value = Xs[i]) for i in range(len(Xs))]

	subset  = data[:, Ys, Xs]
	## statewide median
	medians = np.median(subset, axis = [1])

	pieces = path.rsplit('\\', )
	outname = pieces[1] + '.csv'
	np.savetxt(outname, medians, delimiter=",")
