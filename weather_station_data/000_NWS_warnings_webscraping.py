#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script scrapes historical (1986-present) NWS weather advisories 
# and subsets by NM


#--------------------#
#-- OPEN LIBRARIES --#
#--------------------#
import urllib2
import urllib
from StringIO import StringIO
from zipfile import ZipFile
from urllib import urlopen
import re
import os
import shapefile


#-------------------#
#-- generate URLs --#
#-------------------#
base_url = "https://mesonet.agron.iastate.edu/pickup/wwa/"
list1    = range(1986,2017)

url_list = []
for i in range(len(list1)):
	url_list.append(base_url + str(list1[i]) + "_all.zip")


#--------------------#
#-- download/unzip --#
#--------------------#
for i in range(len(url_list)):
	try:
	    url1 = urlopen(url_list[i])
	    zipfile1 = ZipFile(StringIO(url1.read()))
	    zipfile1.extractall('h:/Jason/Climate/NWS_watches_warnings/' + str(list1[i]))
	except:
		pass


#-----------------------------#
#-- subset shapefiles by NM --#
#-----------------------------#
#get list of shapefiles for import
inpath_list = []
for root, dirs, files in os.walk("h:/Jason/Climate/NWS_watches_warnings"):
    for file in files:
        if file.endswith(".shp"):
             inpath_list.append(os.path.join(root, file))

outpath_list = []
for i in range(len(inpath_list)):
	outpath_list.append(re.sub(".shp", "-NM.shp", inpath_list[i]))

for i in range(len(inpath_list)):
	r = shapefile.Reader(inpath_list[i])
	w = shapefile.Writer(shapeType = shapefile.POLYGON)
	w.fields = list(r.fields)
	selection = [] 
	for rec in enumerate(r.records()):
		if rec[1][10].startswith("NM"):
			selection.append(rec) 
	# Add the geometry and records to the writer
	for rec in selection:
		w._shapes.append(r.shape(rec[0]))
		w.records.append(rec[1])
	w.save(outpath_list[i]) 
