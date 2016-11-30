#------------#
#-- AUTHOR --#
#------------#
# Jason Schatz
# Created:  11/21/16
# Last modified:  11/21/16


#-----------------#
#-- DESCRIPTION --#
#-----------------#
# Generates URLs for downloading NLDAS2 forcing data for 
# New Mexico for calculating daily min/max temperature,
# apparent temprature, and heat index at 1/8th degree 
# resolution. This is the data source used by the CDC, 
# so we are using it for consistency.


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


#-----------------------#
#-- GENERATE URL LIST --#
#-----------------------#
## because the URLs reference an "http services" dest, rather than the file itself, the URLs generated
## here were used in the Firefox extension "DownThemAll" rather than as a direct scrape. This was the
## recommendation on the gsfc website, so it seemed simplest.
##
## data source: http://disc.sci.gsfc.nasa.gov/uui/datasets/NLDAS_FORA0125_H_V002/summary?keywords=NLDAS
## date range:  (1981-01-01 to present{2016-10-30})
## spatial bounding box for New Mexico: 31.31,-109.2,37.07,-102.9
## username:  jason.schatz
## password:  0Uamdw!pw&w

start = "http://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FNLDAS%2FNLDAS_FORA0125_H.002%2F"
middle = ".002.grb&FORMAT=bmV0Q0RGLw&BBOX=31.31%2C-109.2%2C37.07%2C-102.9&LABEL=NLDAS_FORA0125_H.A"
last = ".002.2016328180428.pss.nc&SHORTNAME=NLDAS_FORA0125_H&SERVICE=SUBSET_GRIB&VERSION=1.02&LAYERS=Eg&DATASET_VERSION=002"

hour_list  = ['0000', '0100', '0200', '0300', '0400', '0500', '0600', '0700', 
              '0800', '0900', '1000', '1100', '1200', '1300', '1400', '1500',
              '1600', '1700', '1800', '1900', '2000', '2100', '2200', '2300']
dates_list = pd.date_range(pd.to_datetime("19810101", format='%Y%m%d'), periods=13087).tolist()

url_list   = []
for i in range(len(dates_list)):
	year  = dates_list[i].year
	month = dates_list[i].month
	day_of_month = dates_list[i].day
	date  = str(year) + str("%02d" % (month)) + str("%02d" % (day_of_month))
	doy   = str("%03d" % (dates_list[i].dayofyear))

	for j in range(len(hour_list)):
		hour = hour_list[j]
		url1 = start + str(year) + "%2F" + str(doy) + "%2FNLDAS_FORA0125_H.A" + str(date) + "." + str(hour) + middle + str(date) + "." + str(hour) + last
		url_list.append(url1)
