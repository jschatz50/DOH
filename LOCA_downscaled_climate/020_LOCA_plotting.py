#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 11/16/16


#---------------#
#--DESCRIPTION--#
#---------------#
# generate maps of select LOCA data


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
from pylab import *
import netCDF4
import fnmatch
import os
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import cartopy
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import gridspec


#-----------------------------#
#-- define mapping function --#
#-----------------------------#
def make_map(item, lower, upper):
	map1 = Basemap(width = 650000, height = 750000,
            resolution = 'l', projection = 'stere',
            lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)
	map1.pcolor(xi, yi, np.squeeze(item), vmin = lower, vmax = upper)
	map1.readshapefile('H:/Jason/GIS/Political_boundaries/NM_counties/NM_counties', 'NM_counties', color = 'w')
	map1.readshapefile('H:/Jason/GIS/Political_boundaries/NM_state_boundary/NM_state', 'NM_state', color = 'black', linewidth = 3)


#------------#
#-- SET UP --#
#------------#

# set directory
path1 = "L:/LOCA_figures"
os.chdir(path1)	#set directory

# create list of all files
thresholds = [90, 95, 100, 105, 110]
metrics = ["means", "maxes", "mins", "p90th", "p10th"]


#-------------------------------------------------#
#-- LOOP THROUGH FILES AND CREATE MAPS FOR EACH --#
#-------------------------------------------------#

# import file lists
for j in range(len(thresholds)):
	inpath = 'L:/LOCA_summaries_NM/d' + str(thresholds[j]) + 'F'
	import_list = []
	for root, dirnames, filenames in os.walk(inpath):
	    for filename in fnmatch.filter(filenames, '*.nc'):
	        import_list.append(os.path.join(root, filename))

	# generate multi-panel plots for each metric and threshold
	for k in range(len(metrics)):
		data = [i for i in import_list if metrics[k] in i]
		list_a = [Dataset(data[i], 'r') for i in range(len(data))]
		data1  = [list_a[i].variables['data'] for i in range(len(list_a))]
		list_b = [np.mean(data1[i], axis = 0) for i in range(len(list_a))]   # mean values for each grid cell across all models

		# set color bar upper and lower limits
		llim = 0
		ulim = int(math.ceil(np.max(list_b)))

		# create 2D array
		lons = list_a[0]['longitude'][:]
		lats = list_a[0]['latitude'][:]
		lon_0 = lons.mean()
		lat_0 = lats.mean()
		lon, lat = np.meshgrid(lons, lats)
		m = Basemap(width = 650000, height = 750000,
		            resolution = 'l',projection = 'stere',
       		     lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)
		xi, yi = m(lon, lat)

		# set figure specs
		fig = plt.figure(figsize = (14, 8)) 
		gs = gridspec.GridSpec(2, 3, width_ratios = [2, 1, 1]) 

		# create subplots
		ax1 = plt.subplot(gs[:, 0])
		title1 = "days/yr over " + str(thresholds[j]) + "F -- 1986-2005"
		ax1.set_title(title1)
		make_map(list_b[0], llim, ulim)
		cs1 = m.pcolor(xi, yi, np.squeeze(list_b[k]), vmin = llim, vmax = ulim)
		cbar1 = m.colorbar(cs1, location = 'bottom', pad = "5%")

		ax2 = plt.subplot(gs[0, 1])
		ax2.set_title("2040-2059 -- RCP 4.5")
		make_map(list_b[3], llim, ulim)

		ax3 = plt.subplot(gs[0, 2])
		ax3.set_title("2040-2059 -- RCP 8.5")
		make_map(list_b[4], llim, ulim)

		ax4 = plt.subplot(gs[1, 1])
		ax4.set_title("2080-2099 -- RCP 4.5")
		make_map(list_b[1], llim, ulim)

		ax5 = plt.subplot(gs[1, 2])
		ax5.set_title("2080-2099 -- RCP 8.5")
		make_map(list_b[2], llim, ulim)

		plt.tight_layout()
		#plt.show()
		outname = metrics[k] + "_annual_" + str(thresholds[j]) + "F_days.png"
		plt.savefig(outname, facecolor = fig.get_facecolor())
		plt.clf()
		plt.cla()
		plt.close()


#-----------------#
#-- SINGLE PLOT --#
#-----------------#

# separate path lists for each metric
metrics = ["means", "maxes", "mins", "p90th", "p10th"]
data = [i for i in import_list if metrics[1] in i]

list_a = [Dataset(data[i], 'r') for i in range(len(data))]
data1  = [list_a[j].variables['data'] for j in range(len(list_a))]
list_b = [np.mean(data1[i], axis = 0) for i in range(len(list_a))]   # mean values for each grid cell across all models

lons = list_a[0]['longitude'][:]
lats = list_a[0]['latitude'][:]

# parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

# use meshgrid to create 2D array
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

# plot data
make_map(0, llim, ulim)

# colorbar
cs = map1.pcolor(xi, yi, np.squeeze(list_b[0]), vmin = llim, vmax = ulim)
cbar = map1.colorbar(cs, location='bottom', pad="10%")

# title
plt.title('Days over 90F per year')
plt.show()
