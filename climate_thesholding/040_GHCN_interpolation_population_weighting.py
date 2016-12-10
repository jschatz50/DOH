#------------#
#-- AUTHOR --#
#------------#
# Jason Schatz
# Created:  12/08/16
# Last modified:  12/08/16


#-----------------#
#-- DESCRIPTION --#
#-----------------#
# This script calculates temperature NCDMs using GHCN station data
# for New Mexico. Specifically, it samples each climate 
# parameter (e.g. days per year over 90F) by census tract then 
# uses those tract values to calculate a population-weighted 
# mean for each county.  For example, if County A consista of 
# Tract 1 (population 1000) and Tract 2 (population 2000), the 
# Tract 2 spatial average would be weighted twice as much as the 
# Tract 1 average when calculating County A's mean.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
from functools import partial
import glob
import numpy as np
import os
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import pandas as pd
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import rasterstats
from rasterstats import zonal_stats, point_query
import sklearn
from sklearn import linear_model
from sklearn.linear_model import LinearRegression


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### runs regression kriging on a set of variables, grids, etc. XY coords must be
#### named 'LONGITUDE' and 'LATITUDE', respectively.
## @' param data:        raw data
## @' param grid:        interpolation grid containing same predictors/column names as data
## @' param predictors:  list of predictors...e.g. ['pred1', 'pred2', 'pred3']
## @' param response:    response variable (column of data being predicted)
def regression_krige(data, grid, predictors, response):
	## fit linear regression
	X = data[predictors]
	y = data[[response]]
	lm1 = LinearRegression()
	lm1.fit(X, y)
	## use lm to predict on interpolation grid
	gridx = grid.LONGITUDE.unique()
	gridy = grid.LATITUDE.unique()
	grid_predictions = LinearRegression.predict(lm1, grid[predictors])
	grid_predictions.shape = (gridy.shape[0], gridx.shape[0])
	## krige station residuals
	data['residuals'] = y - LinearRegression.predict(lm1, X)
	OK = OrdinaryKriging(data[['LONGITUDE']], 
		                 data[['LATITUDE']],
		                 data[['residuals']], 
		                 variogram_model = 'linear',
	                     verbose = False, 
	                     enable_plotting = False)
	z, ss = OK.execute('grid', gridx, gridy)
	## add the regression result grid to the kriged residuals grid
	rk = grid_predictions + z
	return rk

#### generate raster from array
## @' param data:        array values
## @' param grid:        contains array coordinates (named 'LATITUDE' and 'LONGITUDE')
## @' param outpath:     path and name of exported file
def create_raster(data, grid, outpath):
	gridx = grid.LONGITUDE.unique()
	gridy = grid.LATITUDE.unique()
	xmin, ymin, xmax, ymax = [gridx.min(), gridy.min(), gridx.max(), gridy.max()]
	nrows, ncols = np.shape(rk)
	xres = (xmax - xmin) / float(ncols)
	yres = (ymax - ymin) / float(nrows)
	output_raster = gdal.GetDriverByName('GTiff').Create(outpath,
                    ncols, nrows, 1, gdal.GDT_Float32)
	#geotransform = (xmin, xres, 0, ymax, 0, -yres)   
	#output_raster.SetGeoTransform(geotransform)  
	output_raster.SetGeoTransform(template_raster.GetGeoTransform())
	output_raster.SetProjection(template_raster.GetProjection())
	output_raster.GetRasterBand(1).WriteArray(data)
	del(output_raster)


#---------------#
#-- read data --#
#---------------#
df1 = pd.read_csv('H:/Jason/Climate/past_climate/GHCND/nm_ghcnd_data_1892-2016.csv')
station_traits = pd.read_csv('H:/Jason/Climate/past_climate/GHCND/nm_ghcn_station_traits.csv')
station_traits = station_traits[['SID', 'LATITUDE', 'LONGITUDE', 'ELEVATION']]
grid = pd.read_csv('H:/Jason/Climate/past_climate/GHCND/interpolations/interpolation_grid.csv')
grid.columns = ['LATITUDE', 'LONGITUDE', 'ELEVATION']
template_raster = gdal.Open("H:/Jason/Climate/temperature_NCDMs/from_GHCN/CRS_template.tif")


#------------------------------------#
#-- prep population weighting data --#
#------------------------------------#
df2 = pd.read_csv('L:/NLDAS/Annual_totals_by_county/PERMANENT/tracts_by_geography_POP.csv')
grouped = df2.groupby(['COUNTYFP'], as_index = False)
countyPOP = grouped['POP'].aggregate(np.sum)
grouped = df2.groupby(['SMALL_AREA'], as_index = False)
small_area_pop = grouped['POP'].aggregate(np.sum)
df2 = pd.merge(df2,countyPOP, how = 'left', 
	           left_on = ['COUNTYFP'], right_on = ['COUNTYFP'])   
df2 = pd.merge(df2,small_area_pop, how = 'left', 
	           left_on = ['SMALL_AREA'], right_on = ['SMALL_AREA'])
df2['county_weight'] = df2['POP_x'] / df2['POP_y']
df2['sa_weight'] = df2['POP_x'] / df2['POP']

# final column names for results dataframes
CTY_names = ['COUNTYFP']
CTY_names.extend(range(1981, 2016))
SA_names = ['SMALL_AREA']
SA_names.extend(range(1981, 2016))


#------------------------------------#
#-- compute thresholds/percentiles --#
#------------------------------------#
## Calculate temperature percentiles (May-Sep, 1981-2010). Stations
## must have 75% data coverage for entire period to be included.
subset = df1[(df1['YEAR'] >= 1981) & (df1['YEAR'] <= 2010) & 
	         (df1['MONTH'] >= 5) & (df1['MONTH'] <= 9)]
subset = subset[['SID', 'TMAX']]
counts = subset.groupby(['SID']).agg(['count'])
counts['percent'] = counts['TMAX'] / float(np.max(counts['TMAX']))
SIDs = counts.loc[counts['percent'] >= 0.75].index.values.tolist()
subset = subset[subset['SID'].isin(SIDs)]
subset = subset.dropna()

percentiles = [90, 95, 98]
pers = {'P {}'.format(q): partial(np.percentile, q = q) for q in percentiles}
thresholds = subset.groupby('SID').agg({'TMAX': pers})
thresholds.reset_index(inplace = True)  
thresholds.columns = thresholds.columns.droplevel()
thresholds.columns = ['SID', 'p95th', 'p98th', 'p90th']

## final thresholds dataframe
thresholds['d90F']  =  90
thresholds['d95F']  =  95
thresholds['d100F'] = 100
thresholds['d105F'] = 105


#-----------------------#
#-- threshold by year --#
#-----------------------#
## calculate number of days per year each station spends above each threshold
threshold_list = ['d90F', 'd95F', 'd100F', 'd105F', 'p90th', 'p95th', 'p98th']
years = range(1981, 2016)
list_a = []
for year in years:
	subset = df1[(df1.YEAR == year)][['SID', 'TMAX']]
	counts = subset.groupby(['SID']).agg(['count'])
	counts['percent'] = counts['TMAX'] / float(np.max(counts['TMAX']))
	SIDs = counts.loc[counts['percent'] >= 0.9].index.values.tolist()
	subset = subset[subset['SID'].isin(SIDs)]

	## threshold 
	subset = pd.merge(subset, thresholds, how = 'left', 
	                  left_on = ['SID'], right_on = ['SID'])
	for thresh in threshold_list:
		subset[thresh] = subset['TMAX'] >= subset[thresh]

	## sum occurrences per year
	totals = subset.groupby(['SID']).agg(['sum'])
	totals.reset_index(inplace = True)  
	totals.columns = ['SID', 'TMAX', 'p95th', 'p98th', 'p90th', 
	                  'd90F', 'd95F', 'd100F', 'd105F']
	totals['YEAR'] = year
	list_a.append(totals)

result = pd.concat(list_a)
result.drop('TMAX', axis = 1, inplace = True)
result = pd.merge(result, station_traits, how = 'left', left_on = ['SID'], right_on = ['SID'])


#---------------------------------------#
#-- regression kriging and population --#
#-- weighted sampling of outputs      --#
#---------------------------------------#
## kriges thresholding results by year, samples by census tract, uses tract values
## to incorporate population weighting into areal averages for counties/small areas.

for thresh in threshold_list:
	## create empty list for results
	county_results = []
	sa_results = []	
	
	for year in years:
		subset = result[(result.YEAR == year)]

		## regression kriging
		rk = regression_krige(data = subset, 
		  	                  grid = grid, 
		                      predictors = ['LATITUDE', 'ELEVATION'], 
		                      response = thresh)
		rk[rk < 0] = 0   #replace 'negative days' values with 0

		## create raster
		create_raster(data = rk, 
			          grid = grid, 
			          outpath = 'H:/Jason/Climate/temperature_NCDMs/from_GHCN/myraster.tif')

		## calculate population-weighted areal means for NM counties and small areas
		# sample raster by census tract
		stats = zonal_stats('H:/Jason/GIS/Census/geometries/NM_tracts.shp', 
			                'H:/Jason/Climate/temperature_NCDMs/from_GHCN/myraster.tif', 
					         band = 1, all_touched = True, geojson_out = True)
		means  = [stats[j]['properties']['mean'] for j in range(len(stats))]
		tracts = [stats[j]['properties']['GEOID'] for j in range(len(stats))]
		result2 = pd.DataFrame(tracts, columns = ['TRACT'])
		result2['means'] = means
		result2 = result2.convert_objects(convert_numeric = True)
		df3 = pd.merge(df2, result2, how = 'left', left_on = ['TRACT'], right_on = ['TRACT'])

		# calculate population-weighted temeprature metric for each NM county and small area
		df3['wtd_mean_cnty'] = df3['means'] * df3['county_weight']
		df3['wtd_mean_sa'] = df3['means'] * df3['sa_weight']
		grouped = df3.groupby(['COUNTYFP'], as_index = False)
		county_means = grouped['wtd_mean_cnty'].aggregate(np.sum)
		grouped = df3.groupby(['SMALL_AREA'], as_index = False)
		sa_means = grouped['wtd_mean_sa'].aggregate(np.sum)
		county_results.append(county_means)
		sa_results.append(sa_means)

	## merge results lists into final file, set column names, write to file
	county_final = reduce(lambda left, right: pd.merge(left, right, on = 'COUNTYFP'), county_results)
	sa_final = reduce(lambda left, right: pd.merge(left, right, on = 'SMALL_AREA'), sa_results)
	county_final.columns = CTY_names
	sa_final.columns = SA_names

	CTY_out = 'H:/Jason/Climate/temperature_NCDMs/from_GHCN/results/' + \
	          'county' + '_' + thresh + '.csv'
	SA_out  = 'H:/Jason/Climate/temperature_NCDMs/from_GHCN/results/' + \
	          'small_area' + '_' + thresh + '.csv'

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
		item1['metric'] = which_list[i].rsplit('_', 5)[3].replace('.csv', '')   #3 for counties, 4 for SAs
		list_a.append(item1)
	result = pd.concat(list_a)
	result.to_csv(outname, index = False)

## path lists for counties/small areas
path1 = "H:/Jason/Climate/temperature_NCDMs/from_GHCN/results"
os.chdir(path1)	#set directory
import_list = glob.glob('H:/Jason/Climate/temperature_NCDMs/from_GHCN/results/individual/*.csv')
county_list = [s for s in import_list if 'county' in s]
sa_list     = [s for s in import_list if 'small_area' in s]

merge_csvs(county_list, 'COUNTYFP', 'NM_county_temperature_NCDMs_1981-2015.csv')
merge_csvs(sa_list, 'SMALL_AREA', 'NM_small_area_temperature_NCDMs_1981-2015.csv')
