#------------#
#-- AUTHOR --#
#------------#
# Jason Schatz
# Created:  12/06/16
# Last modified:  12/07/16


#-----------------#
#-- DESCRIPTION --#
#-----------------#
# This script aggregates NLDAS-derived climate data by county and 
# New Mexico small area.  Specifically, it samples each climate 
# parameter (e.g. days per year over 90F) by census tract then 
# uses those tract values to calculate a population-weighted 
# mean for each county.  For example, if County A consista of 
# Tract 1 (population 1000) and Tract 2 (population 2000), the 
# Tract 2 spatial average would be weighted twice as much as the 
# Tract 1 average when calculating County A's mean.
#
# NOTE: the number of days over various threshold is typically
# higher for Tair than for HI and AT.  This may seem like an 
# error, but it's not. NM just has very low humidity, so AT and 
# HI are very often lower than Tair.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
import datetime as dt
import glob
import math
import netCDF4
from netCDF4 import Dataset
import numpy as np
import os
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import pandas as pd
import rasterstats
from rasterstats import zonal_stats, point_query
import re
from StringIO import StringIO


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### Convert a netcdf to a raster file
## @' param step:  step in loop
## @' param outname:  name of output (e.g. 'temp_raster.tif')
## @' param template_raster:  existing raster for grabbing spatial info
def ncdf2raster(step, outname, template_raster):
	try:
		nrows,ncols = np.shape(data1[step, :, :])
		xres = (xmax - xmin) / float(ncols)
		yres = (ymax - ymin) / float(nrows)
#		geotransform = (xmin, xres, 0, ymin, 0, -yres)
		output_raster = gdal.GetDriverByName('GTiff').Create(outname,ncols, nrows, 1, gdal.GDT_Float32)
#		output_raster.SetGeoTransform(geotransform)
		output_raster.SetGeoTransform(template_raster.GetGeoTransform())
		output_raster.SetProjection(template_raster.GetProjection())
		inverted_data = data1[step,:,:][::-1,:]
		output_raster.GetRasterBand(1).WriteArray(inverted_data)
	finally:
		del(output_raster)


#-----------------#
#-- import data --#
#-----------------#
path1="L:/NLDAS/Annual_totals_by_county/results"
os.chdir(path1)	#set directory
import_list = glob.glob('L:/NLDAS/Annual_summaries/*.nc')


#----------------------------------------#
#-- prep data for population weighting --#
#----------------------------------------#
df1 = pd.read_csv('L:/NLDAS/Annual_totals_by_county/PERMANENT/tracts_by_geography_POP.csv')
grouped = df1.groupby(['COUNTYFP'], as_index = False)
countyPOP = grouped['POP'].aggregate(np.sum)
grouped = df1.groupby(['SMALL_AREA'], as_index = False)
small_area_pop = grouped['POP'].aggregate(np.sum)
df1 = pd.merge(df1,countyPOP, how = 'left', 
	           left_on = ['COUNTYFP'], right_on = ['COUNTYFP'])   
df1 = pd.merge(df1,small_area_pop, how = 'left', 
	           left_on = ['SMALL_AREA'], right_on = ['SMALL_AREA'])
df1['county_weight'] = df1['POP_x'] / df1['POP_y']
df1['sa_weight'] = df1['POP_x'] / df1['POP']

# final column names for results dataframes
CTY_names = ['COUNTYFP']
CTY_names.extend(range(1981,2017))
SA_names = ['SMALL_AREA']
SA_names.extend(range(1981,2017))


#----------------------------------------#
#-- calculate pop-weighted areal means --#
#-- for NM counties and small areas    --#
#----------------------------------------#
for f in xrange(len(import_list)):

	## import netcdf
	nc_file = import_list[f]
	fh    = Dataset(nc_file, mode = 'r')
	lons  = fh.variables['longitude'][:]
	lats  = fh.variables['latitude'][:]
	dates = fh.variables['time'][:]
	data1 = fh.variables['data'][:]
	fh.close()

	## set raster parameters
	inDs = gdal.Open("L:/NLDAS/Annual_totals_by_county/PERMANENT/CRS_raster.tif")
	xmin,ymin,xmax,ymax = [lons.min(), lats.min(), lons.max(), lats.max()]

	## create empty list for results
	county_results = []
	sa_results = []

    ## for each year...
	for i in xrange(len(dates)):

		## transform netcdf into raster
		ncdf2raster(i, 'temp_raster.tif', inDs)
		
		## sample raster by census tract
		stats = zonal_stats('H:/Jason/GIS/Census/geometries/NM_tracts.shp', 
			                'L:/NLDAS/Annual_totals_by_county/results/temp_raster.tif', 
			                band = 1, all_touched = True, geojson_out = True)
		means  = [stats[j]['properties']['mean'] for j in range(len(stats))]
		tracts = [stats[j]['properties']['GEOID'] for j in range(len(stats))]
		result = pd.DataFrame(tracts, columns = ['TRACT'])
		result['means'] = means
		result = result.convert_objects(convert_numeric = True)
		df2 = pd.merge(df1, result, how = 'left', left_on = ['TRACT'], right_on = ['TRACT'])

		## calculate population-weighted temeprature metric for each NM county and small area
		df2['wtd_mean_cnty'] = df2['means'] * df2['county_weight']
		df2['wtd_mean_sa'] = df2['means'] * df2['sa_weight']
		grouped = df2.groupby(['COUNTYFP'], as_index = False)
		county_means = grouped['wtd_mean_cnty'].aggregate(np.sum)
		grouped = df2.groupby(['SMALL_AREA'], as_index = False)
		sa_means = grouped['wtd_mean_sa'].aggregate(np.sum)
		county_results.append(county_means)
		sa_results.append(sa_means)

	## merge results lists into final file, set column names, write to file
	county_final = reduce(lambda left,right: pd.merge(left,right, on = 'COUNTYFP'), county_results)
	sa_final = reduce(lambda left,right: pd.merge(left,right, on = 'SMALL_AREA'), sa_results)
	county_final.columns = CTY_names
	sa_final.columns = SA_names

	CTY_out = 'county' + '_' + import_list[f][26:] + '.csv'
	SA_out  = 'small_area' + '_' + import_list[f][26:] + '.csv'

	county_final.to_csv(CTY_out, index = False)
	sa_final.to_csv(SA_out, index = False)


#-------------------------------#
#-- merge CSVs to single file --#
#-------------------------------#
#### merge list of individual .csv files to one big file
## @' param which_list:  list of .csv file paths
## @' param id_var:  id variable for counties/small areas (COUNTYFP or SMALL_AREA)
## @' param outname:  output name for merged csv
def merge_csvs(which_list, id_var, outname):
	list_a = []
	for i in range(len(which_list)):
		item1 = pd.read_csv(which_list[i])
		item1 = pd.melt(item1, id_vars=[id_var])
		item1.columns = [id_var, 'year', 'ndays']
		item1['metric'] = which_list[i].rsplit('_', 4)[1] + '_' + \
		                  which_list[i].rsplit('_', 3)[3].replace('.nc.csv', '')

		list_a.append(item1)
	result = pd.concat(list_a)
	result.to_csv(outname, index = False)

## path lists for counties/small areas
path1="L:/NLDAS/NCDMs"
os.chdir(path1)	#set directory
import_list = glob.glob('L:/NLDAS/Annual_totals_by_county/results/*.csv')
county_list = [s for s in import_list if 'county_' in s]
sa_list     = [s for s in import_list if 'small_area_' in s]

## merge county and small area CSVs into individual files
merge_csvs(county_list, 'COUNTYFP', 'NM_county_temperature_NCDMs_1981-2015.csv')
merge_csvs(sa_list, 'SMALL_AREA', 'NM_small_area_temperature_NCDMs_1981-2015.csv')
