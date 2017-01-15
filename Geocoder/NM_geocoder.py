#------------#
#-- AUTHOR --#
#------------#
## Jason Schatz
## Created: 12/30/2016
## Last modified: 01/04/2017


#---------------#
#-- FILE INFO --#
#---------------#
## NM address geocoder.  This script:
##
## (1) separates addresses into street addresses and PO Boxes
## (2) attempts to geocode the street addresses
## (3) for addresses that are not geocoded, ensures that the zip matches the city (people
##     entering the wrong city is one of the most frequent and easily corrected errors).
##     If original city does not match zip, geocoding is reattempted with correct city.
## (4) For all successful geocoded addresses, assigns them to the appropriate census 
##     tract, county, zip, city, doh health region, and NM small area.
## (5) For all unsuccessfully geocoded addresses, each location is assigned (based on zip)
##     to the appropriate tract/county/zip/etc, to the extent possible.  Some zips are
##     entirely contained within a census tract, small area, county, etc.  So in these
##     cases, it is possible to assign some geographic specificity to the location,
##     assuming that the zip code is entered correctly
## (6) Outputs fully and partially geocoded addresses
## (7) Uses affinity propogagtion to quantify address similarity to see if any non-geocoded 
##     addresses cluster with any geocoded addresses,
##
##
## Addresses should be in the following format:
## street        city        state        zip         address
## 123 Any St    SANTA FE    NM           87505       123 Any St Santa Fe NM 87505
##
##
## Script uses Google's API, which is free up to 2500 request/day. If billing is enabled,
## this script will run much faster (50 requests/second rather than 1 per second) and 
## would cost $0.50 per 1000 requests above 2500/day.  So 10,000 requests costs ~$3.75.
## 
## https://developers.google.com/maps/documentation/geocoding/intro#Types


#----------------------#
#-- IMPORT LIBRARIES --#
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
#-- DEFINE FUNCTIONS --#
#----------------------#

#############################
#### geocoding functions ####
#############################
#### FUNCTION 1: separate addresses into PO Boxes and street addresses
## @param addresses:  dataframe with 'address' column
def box_or_street(addresses):
	addresses = addresses['address'].tolist()
	addresses = [i.upper() for i in addresses]
	## remove unnecessary characters from addresses
	addresses = map(lambda x: x.replace('TRLR ', ''), addresses)   #remove trailer
	addresses = map(lambda x: x.replace('SP ', ''), addresses)   #remove space
	addresses = map(lambda x: x.replace('NUM ', ''), addresses)   #remove 'NUM' from address
	pobox_list = [i for i in addresses if 'BOX' in i]
	street_list = [i for i in addresses if i not in pobox_list]
	return street_list, pobox_list

#### FUNCTION 2: geocode using google API
## @param addresses: list of addresses
## @param original_addresses: dataframe of original addresses for comparison with output
# limited to 2500 request per day, though this can be 
# increased for $0.50/1000 additional requests
def geocode(address_list, original_addresses):
	d = {'original':[], 'lat':[], 'lon':[], 
	     'county':[], 'city':[], 'zip':[], 
	     'address':[], 'status':[], 'confidence':[], 
	     'quality':[], 'accuracy':[]}
	for address in address_list:
		try:
			g = geocoder.google(address)
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
			d['original'].append(address)
		except:
			pass
	df = pd.DataFrame.from_dict(d)
	street_list = [i for i in address_list if i not in df['original'].tolist()]
	unsuccessful = original_addresses[original_addresses['address'].isin(street_list)]
	return df, unsuccessful

#### FUNCTION 3: find the city that corresponds with a given zip code 
#### and re-geocode if original city was incorrect
## @param unsuccessful: dataframe with zip codes and cities listed
def geocode_retry(unsuccessful):
	zip_list = unsuccessful['zip'].tolist()
	city_list = [geocoder.google(zipcode).city for zipcode in zip_list]
	unsuccessful['city2'] = city_list
	unsuccessful['street'] = unsuccessful['street'].str.upper()
	unsuccessful['city'] = unsuccessful['city'].str.upper()
	unsuccessful['city2'] = unsuccessful['city2'].str.upper()

	retry_set = unsuccessful.loc[unsuccessful['city'] != unsuccessful['city2']]
	retry_set['address2'] = retry_set['street'] + ' ' + retry_set['city2'] + ' NM ' + retry_set['zip'].map(str)
	address_list = retry_set['address2'].tolist()

	df, unsuccessful = geocode(address_list, retry_set)
	return df, unsuccessful

#### FUNCTION 4: get public health area, county, zip, small area, tract for each location
## @param coords: dataframe of locations to sample (colnames 'lon' and 'lat')
def assign_geographies(coords):
	small_areas = gpd.read_file('P:/Jason/Geocoder/geographies/small_areas.shp')
	doh_regions = gpd.read_file('P:/Jason/Geocoder/geographies/doh_regions.shp')
	tracts = gpd.read_file('P:/Jason/Geocoder/geographies/tracts.shp')
    counties = gpd.read_file('P:/Jason/Geocoder/geographies/counties.shp')

	crs = {'init': u'epsg:4326'}
	geometry = [Point(xy) for xy in zip(coords.lon, coords.lat)]
	geo_df = gpd.GeoDataFrame(coords, crs=crs, geometry=geometry)

	coords['small_area'] = gpd.sjoin(geo_df, small_areas, how="left", op='intersects')['OBJECTID']
	coords['doh_region'] = gpd.sjoin(geo_df, doh_regions, how="left", op='intersects')['Reg_name']
	coords['tract'] = gpd.sjoin(geo_df, tracts, how="left", op='intersects')['GEOID']
	coords['county'] = gpd.sjoin(geo_df, counties, how="left", op='intersects')['NAME']
	return coords

#### FUNCTION 5: get county, census tracts, and NM small areas from zip codes (when possible)
## @param unsuccessful: unsuccessfully geocoded addresses/po boxes
def zips_to_geographies(unsuccessful):
	key = pd.read_csv('P:/Jason/Geocoder/geography_keys/zip_key.csv')
	merged = pd.merge(unsuccessful, key, how = 'left',  left_on = ['zip'], right_on = ['zip'])
	return merged

#### FUNCTION 6: master function (links sub-functions in sequence)
## @param addresses:  dataframe with 'address' column
## @param outpath:  desired path for output files to be saved
def NM_geocoder(addresses, outpath):
	## addresses to uppercase
	addresses['address'] = addresses['address'].str.upper()

	## separate street addresses from po boxes
	address_list, box_list = box_or_street(addresses)
	box_addresses = addresses[addresses['address'].isin(box_list)]

	## geocode street addresses
	geocoded1, unsuccessful1 = geocode(address_list, addresses)

	## check if city matches zip; re-attempt geocoding after correction
	geocoded2, unsuccessful2 = geocode_retry(unsuccessful1)

	## combine results into dataframes
	geocoded = pd.concat([geocoded1, geocoded2], axis = 0)
	geocoded = geocoded.set_index([range(0, len(geocoded))])
	unsuccessful = pd.concat([box_addresses, unsuccessful1, unsuccessful2], axis = 0)
	unsuccessful = unsuccessful.set_index([range(0, len(unsuccessful))])

	## assign geographies to locations
	geocoded = assign_geographies(geocoded)
	semi_geocoded = zips_to_geographies(unsuccessful)

	## write to file
	geocoded.to_csv(outpath + '/geocoded_address_results.csv', index = False, encoding='utf-8')
	semi_geocoded.to_csv(outpath + '/partially_geocoded_address_results.csv', index = False)

	return geocoded, semi_geocoded


###############################################
#### address similarity matching functions ####
###############################################
#### FUNCTION 1: find all matching indices in a list (sub)
## @param lst:  list
## @a:  value to match
def find(lst, a):
    return [i for i, x in enumerate(lst) if x == a]

#### FUNCTION 2: match addresses one by one
## @param master: geocoded addresses
## @param target: addresses that were not successfully geocoded
## @param outpath: desired path for output files to be saved
def address_match(master, target, outpath):
	master_list = master['address'].tolist()
	master_list = [i.upper() for i in master_list]
	target_list = target['address'].tolist()
	target_list = [i.upper() for i in target_list]

	results_list = []
	for address in target_list:
		list_a = [jellyfish.damerau_levenshtein_distance(unicode(address), unicode(master_address)) for master_address in master_list]
		nearest_addresses = [master_list[i] for i in find(list_a, np.min(list_a))]
		df1 = pd.DataFrame({'nearest_addresses': nearest_addresses})
		df1['master_address'] = address
		df1['damerau_levenshtein_distance'] = np.min(list_a)
		results_list.append(df1)
	final = pd.concat(results_list, axis = 0)
	final.to_csv(outpath + '/best_matches_for_non-geocoded_addresses.csv', index = False)
	return final

#### FUNCTION 3: use affinity propagation to cluster addresses by similarity
## @param addresses:  dataframe with 'address' column
def address_cluster(addresses):
	addresses = np.asarray(addresses['address'])
	similarity = -1 * np.array([[jellyfish.damerau_levenshtein_distance(unicode(w1), unicode(w2)) for w1 in addresses] for w2 in addresses])
	affprop = sklearn.cluster.AffinityPropagation(affinity='precomputed', 
		                                          damping=0.5)
	affprop.fit(similarity)

	list_a = []
	for cluster_id in np.unique(affprop.labels_):
		exemplar = addresses[affprop.cluster_centers_indices_[cluster_id]]
		cluster = addresses[np.nonzero(affprop.labels_ == cluster_id)]
		affinity = similarity[np.flatnonzero(addresses == exemplar).tolist()[0], np.flatnonzero(affprop.labels_ == cluster_id).tolist()]
		df1 = pd.DataFrame({'similar': cluster, 'affinity': affinity, 'exemplar': str(exemplar)})
		list_a.append(df1)
	return pd.concat(list_a)


#---------------#
#-- geocoding --#
#---------------#
## read data
addresses = pd.read_csv('P:/Jason/Geocoder/sample_data/2015_ED_EPHT_testset.csv')

## geocode
start_time = time.clock()
geocoded, semi_geocoded = NM_geocoder(addresses = addresses,
	                                  outpath = 'P:/Jason/Geocoder/sample_data')
print(" %s seconds" % (time.clock() - start_time))


#----------------------#
#-- address matching --#
#----------------------#
## get best matches for each semi-geocoded address
best_matches = address_match(master = geocoded, 
	                         target = semi_geocoded,
	                         outpath = 'P:/Jason/Geocoder/sample_data')


## address clusters (if using 1000s of addresses, this takes a long time, but works)
addresses = pd.read_csv('P:/Jason/Geocoder/sample_data/radon_addresses.csv')
clusters = address_cluster(addresses)


#----------------------------#
#-- check for missing data --#
#----------------------------#
df1 = pd.read_csv('P:/Jason/Geocoder/sample_data/2015_ED_EPHT_testset.csv')
df2 = pd.read_csv('P:/Jason/Geocoder/sample_data/geocoded_address_results.csv')
df3 = pd.read_csv('P:/Jason/Geocoder/sample_data/partially_geocoded_address_results.csv')

list1 = df2["original"].tolist()
list2 = df3['address'].tolist()
output_list = list1 + list2

input_list = df1['address'].tolist()

missing = [i for i in input_list if i not in output_list]
