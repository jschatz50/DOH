#------------#
#-- author --#
#------------#
# Jason Schatz
# Created:  12.30.2016
# Modified: 01.26.2017 - now eliminates duplicate addresses before geocoding
# Modified: 02.10.2017 - added pause if daily limit is exceeded
#                      - cleaned up docstrings
# Modified: 03.22.2017 - simplified code


#-----------------#
#-- description --#
#-----------------#
# NM address geocoder.  This script:
#
#  (1) Separates addresses into street addresses and PO Boxes.
#  (2) Attempts to geocode the street addresses.
#  (3) For addresses that are not geocoded, ensures that the zip matches the city (people
#      entering the wrong city is one of the most frequent and easily corrected errors).
#      If original city does not match zip, geocoding is reattempted with correct city.
#  (4) For all successful geocoded addresses, assigns them to the appropriate census 
#      tract, county, zip, city, doh health region, and NM small area.
#  (5) For all unsuccessfully geocoded addresses, each location is assigned (based on zip)
#      to the appropriate tract/county/zip/etc, to the extent possible.  Some zips are
#      entirely contained within a census tract, small area, county, etc.  So in these
#      cases, it is possible to assign some geographic specificity to the location,
#      assuming that the zip code is entered correctly
#  (6) Outputs fully and partially geocoded addresses
#
# Addresses should be in the following format:
#  street        city        state   zip     address
#  123 Any St    SANTA FE    NM      87505   123 Any St Santa Fe NM 87505
#
# Script uses Google's API, which is free up to 2500 request/day. If billing is enabled,
# this script will run much faster (50 requests/second rather than 1 per second) and 
# would cost $0.50 per 1000 requests above 2500/day.  So 10,000 requests costs ~$3.75.
# 
# https://developers.google.com/maps/documentation/geocoding/intro#Types


#----------------------#
#-- import libraries --#
#----------------------#
import geocoder
import geopandas as gpd
import jellyfish
import pandas as pd
import numpy as np
from shapely.geometry import Point
import sklearn.cluster
import time


#----------------------#
#-- define functions --#
#----------------------#
def prep_addresses(addresses):
	'''Sort street addresses from PO boxes and preps data for subsequent functions

	Args:
		addresses: a dataframe with columns: 'street' 'city' 'state' 'zip' 'address'

	Returns:
		dataframe of formatted street addresses
	'''

	## to uppercase (for consistency/finding duplicates)
	addresses['address'] = addresses['address'].astype(str).str.upper()

	## remove po boxes
	addresses = addresses[~addresses['address'].str.contains('PO BOX|P.O. BOX|P O BOX|P. O. BOX', na=False)]

	## remove common problematic characters from addresses
	addresses['address'] = addresses['address'].str.replace('TRLR', '')
	addresses['address'] = addresses['address'].str.replace('SP', '')
	addresses['address'] = addresses['address'].str.replace('NUM', '')
	return addresses


def geocode(addresses, originals):
	'''Geocode addresses using Google Geocoding API

	Args:
		addresses: dataframe output from 'prep_addresses' function
		originals: original address dataframe

	Returns:
		dataframe of geocoded addresses
	'''

	## prep data dictionary
	d = {'original_street':[], 'original_zip':[], 'original_address':[],
	     'lat':[], 'lon':[], 'county':[], 'city':[], 
	     'zip':[], 'address':[], 'status':[], 
	     'confidence':[], 'quality':[], 'accuracy':[]}
	
	## initialize geocoding status
	status = 'OK'

	## reduce to unique addresses
	addresses = addresses.drop_duplicates(['address'], keep='first')
	addresses.reset_index(drop=True, inplace=True)

	## geocode each unique address
	for idx in range(len(addresses)):
		g = geocoder.google(addresses['address'][idx])
		
		status = g.status
		if status == 'OVER_QUERY_LIMIT':
			print 'over query limit; pausing for 1 hr'
			time.sleep(60*60)

		try:
			d['lat'].append(g.latlng[0])
			d['lon'].append(g.latlng[1])
			d['county'].append(g.county)
			d['city'].append(g.city)
			d['zip'].append(g.postal)
			d['address'].append(g.address)
			d['status'].append(g.status)
			d['confidence'].append(g.confidence)
			d['quality'].append(g.quality)
			d['accuracy'].append(g.accuracy)
			d['original_street'].append(addresses['street'][idx])
			d['original_address'].append(addresses['address'][idx])
			d['original_zip'].append(addresses['zip'][idx])
		except:
			pass

	geocoded = pd.DataFrame.from_dict(d)
	unsuccessful = [i for i in addresses['address'] if i not in geocoded['original_address'].tolist()]
	unsuccessful = addresses[addresses['address'].isin(unsuccessful)]
	return geocoded, unsuccessful


def geocode_retry(geocoded, unsuccessful):
	'''Find the city that corresponds with a given zip code and 
	re-geocode if original city was incorrect (this is common)

	Args:
		geocoded: dataframe of successfully geocoded addresses 
			      (from 'geocode' function)
		unsuccessful: dataframe with zip codes and cities of street 
		              addresses that were not successfullly geocoded
		              in 'geocode' function

	Returns:
		dataframe of any newly successfully geocoded addresses
	'''

	## return unique zips and geocode corresponding city
	zip_list = unsuccessful['zip'].tolist()
	unique_zips = list(set(zip_list))
	city_list = [geocoder.google(zipcode).city for zipcode in unique_zips]
	zips_cities = pd.DataFrame({'zip': unique_zips, 'city2': city_list})

	## return addresses where original city does not match zip
	unsuccessful = pd.merge(unsuccessful, zips_cities, how='left', on=['zip'])
	unsuccessful['street'] = unsuccessful['street'].str.upper()
	unsuccessful['city'] = unsuccessful['city'].str.upper()
	unsuccessful['city2'] = unsuccessful['city2'].str.upper()
	unsuccessful = unsuccessful.drop_duplicates(keep='first')
	retry_set = unsuccessful.loc[unsuccessful['city'] != unsuccessful['city2']]
	retry_set['address2'] = retry_set['street'] + ' ' + retry_set['city2'] + ' NM ' + retry_set['zip'].map(str)
	retry_set.rename(columns = {'address': 'address1', 'address2': 'address'}, inplace=True) 

	## reattempt geocoding
	goods, bads = geocode(retry_set, unsuccessful)

	## replace new address with original address
	matches = pd.merge(goods, unsuccessful, how='left', left_on=['original_street'], right_on=['street'])
	matches = pd.concat([matches['address_x'], matches['address_y']], axis=1)
	goods = pd.merge(goods, matches, how='left', left_on=['address'], right_on=['address_x'])
	goods['original_address'] = goods['address_y']
	goods.drop(['address_x', 'address_y'], axis=1, inplace=True)

	## merge newly geocoded with fully geocoded addresses (if any)
	geocoded = geocoded.append(goods)
	geocoded.reset_index(drop=True, inplace=True)
	return geocoded


def assign_geographies(geocoded):
	'''Use spatial join to assign public health area, county, zip, small area, 
	and census tract for each geocoded location.

	Args:
		geocoded: dataframe of locations to sample (from fun 'geocode_retry';
		          colnames 'lon' and 'lat' required)

	Returns:
		dataframe of geocoded addresses joined with corresponding geometries
	'''	
	small_areas = gpd.read_file('P:/Jason/Geocoder/geographies/small_areas.shp')
	doh_regions = gpd.read_file('P:/Jason/Geocoder/geographies/doh_regions.shp')
	tracts = gpd.read_file('P:/Jason/Geocoder/geographies/tracts.shp')
    counties = gpd.read_file('P:/Jason/Geocoder/geographies/counties.shp')

	crs = {'init': u'epsg:4326'}
	geometry = [Point(xy) for xy in zip(geocoded.lon, geocoded.lat)]
	geo_df = gpd.GeoDataFrame(geocoded, crs=crs, geometry=geometry)

	geocoded['doh_region'] = gpd.sjoin(geo_df, doh_regions, how="left", op='intersects')['Reg_name']
	geocoded['county'] = gpd.sjoin(geo_df, counties, how="left", op='intersects')['NAME']
	geocoded['small_area'] = gpd.sjoin(geo_df, small_areas, how="left", op='intersects')['OBJECTID']
	geocoded['tract'] = gpd.sjoin(geo_df, tracts, how="left", op='intersects')['GEOID']
	geocoded.drop('geometry', axis=1, inplace=True)
	return geocoded


def zips_to_geographies(geocoded, originals):
	'''Use spatial join to assign public health area, county, zip, small area, 
	and census tract for each non-geocoded location (when those geographies 
	contain >=90% of a given zip code).

	Args:
		geocoded: dataframe of geocoded addresses joined with corresponding geometries
		          (from 'assign_geographies' function)
		originals: original address dataframe

	Returns:
		dataframe of non-geocoded addresses joined with corresponding geometries by zip
	'''	

	## return addresses that were not successfully geocoded
	unsuccessful = [i for i in originals['address'] if i not in geocoded['original_address'].tolist()]
	unsuccessful = originals[originals['address'].isin(unsuccessful)]

	## merge zip code/geography keys with unsuccessfully geocoded addresses
	key = pd.read_csv('P:/Jason/Geocoder/geography_keys/zip_key.csv')
	merged = pd.merge(unsuccessful, key, how='left', on=['zip'])
	return merged


def NM_geocoder(addresses, outpath):
	'''Master function (links sub-functions in sequence)

	Args:
		addresses: a dataframe with columns: 'street' 'city' 'state' 'zip' 'address'
		outpath:  desired path for output files to be saved

	Returns:
		dataframe of geocoded addresses joined with corresponding geometries
		dataframe of non-geocoded addresses joined with corresponding geometries
	'''		

	## run functions in sequence
	originals = addresses
	addresses = prep_addresses(addresses)
	geocoded, unsuccessful = geocode(addresses, originals)
	geocoded = geocode_retry(geocoded, unsuccessful)
	geocoded = assign_geographies(geocoded)
	zip_coded = zips_to_geographies(geocoded, originals)

	## ensure outpath is formatted correctly
	if outpath[-1] == '/':
		outpath = outpath
	else:
		outpath = outpath + '/'

	## write results to file
	geocoded.to_csv(outpath + 'geocoded_address_results.csv', index=False, encoding='utf-8')
	zip_coded.to_csv(outpath + 'partially_geocoded_address_results.csv', index=False, encoding='utf-8')
	return geocoded, zip_coded


#---------------#
#-- geocoding --#
#---------------#
## define parameters
addresses = pd.read_csv('P:/Jason/_Project_summaries/Geocoding/_CODE/sample_addresses.csv')
outpath = 'P:/Jason/_Project_summaries/Geocoding/_CODE'

## geocode
geocoded, zip_coded = NM_geocoder(addresses=addresses,
	                              outpath=outpath)
