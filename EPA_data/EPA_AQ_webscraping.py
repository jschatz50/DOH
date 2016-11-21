#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script downloads all EPA daily AQM files for 1990-present, unzips them,
# subsets by NM stations, writes them to pollutant specific folders, then merges
# them into statewide files for each pollutant


#----------------------#
#-- import libraries --#
#----------------------#
import urllib2
import urllib
from StringIO import StringIO
from zipfile import ZipFile
from urllib import urlopen
import pandas as pd
import re
import glob


#-------------------#
#-- generate URLs --#
#-------------------#
list1 = ['44201', '42401', '42101', '42602', '88101', '88502', '81102', 'SPEC',
		'WIND', 'TEMP', 'PRESS', 'RH_DEP', 'HAPS', 'VOCS', 'NONOxNOy', 'LEAD']
list2 = range(1990, 2017)

url_list = []
for i in range(len(list1)):
	url1 = "https://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/daily_" + list1[i]
	for j in range(len(list2)):
		url_list.append(url1 + "_" + str(list2[j]) + ".zip")


#---------------------------------#
#-- download/unzip/subset to NM --#
#---------------------------------#
for i in range(len(url_list)):
	try:
	    url1 = urlopen(url_list[i])
	    zipfile1 = ZipFile(StringIO(url1.read()))
	    zz = re.sub(url_list[i][-9:],b "", 
	    	re.sub("https://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/daily_", "", url_list[i]))
	    zipfile1.extractall('h:/Jason/EPA_AQ/' + zz)
	    a = zipfile1.namelist()[0]
	    df1 = pd.read_csv('h:/Jason/EPA_AQ/' + zz + "/" + a)
	    df1 = df1[(df1['State Name'] == 'New Mexico')]
	    df1.to_csv('h:/Jason/EPA_AQ/' + zz + "/" + a)
	    del df1
	    del zipfile1
	except:
		pass


#-----------------#
#-- merge files --#
#-----------------#
path_list = ['H:/Jason/EPA_AQ/42101/*.csv', 'H:/Jason/EPA_AQ/42401/*.csv', 'H:/Jason/EPA_AQ/42602/*.csv',
			'H:/Jason/EPA_AQ/44201/*.csv', 'H:/Jason/EPA_AQ/81102/*.csv', 'H:/Jason/EPA_AQ/88101/*.csv',
			'H:/Jason/EPA_AQ/88502/*.csv', 'H:/Jason/EPA_AQ/HAPS/*.csv', 'H:/Jason/EPA_AQ/LEAD/*.csv',
			'H:/Jason/EPA_AQ/NONOxNOy/*.csv', 'H:/Jason/EPA_AQ/PRESS/*.csv', 'H:/Jason/EPA_AQ/SPEC/*.csv',
			'H:/Jason/EPA_AQ/TEMP/*.csv', 'H:/Jason/EPA_AQ/VOCS/*.csv', 'H:/Jason/EPA_AQ/WIND/*.csv']

for i in range(len(path_list)):
	import_list = glob.glob(path_list[i])
	list_a = []
	for j in range(len(import_list)):
		list_a.append(pd.read_csv(import_list[j]))
	df = pd.concat(list_a)
	outname = re.sub("/[*].csv", "", path_list[i]) + ".csv"
	df.to_csv(outname)
