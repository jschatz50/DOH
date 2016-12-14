#------------#
#-- AUTHOR --#
#------------#
# Jason Schatz
# Created:  12/09/16
# Last modified:  12/09/16


#-----------------#
#-- DESCRIPTION --#
#-----------------#
# This script calculates merges tabular data with shapefiles
# to map temperature NCDM results.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(tigris)
library(rgdal)
library(sp)
library(leaflet)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
#### makes a leaflet map
## @' param data:  shapefile containing geometry and data
## @' param factor:  factor that you're mapping
## @' param name:  descriptive name of factor you're mapping
make_map  <- function(data, factor, geo_name, name){
    popup <- paste0("NAME: ", unlist(data.frame(data[ , geo_name])), 
                    "<br>", paste(name, ": ", sep = ""), 
                     unlist(round(data.frame(data[ , factor]), 2)))
    pal   <- colorNumeric(palette = "YlGnBu",
             domain = data.frame(data[ , factor]))
    map1  <- leaflet() %>%
             addProviderTiles("CartoDB.Positron") %>%
             addPolygons(data = data, 
                         fillColor = ~pal(data@data[factor]), 
                         color = "#b2aeae",
                         fillOpacity = 0.7, 
                         weight = 1, 
                         smoothFactor = 0.2,
                         popup = popup) %>%
             addLegend(pal = pal, 
                       values = unlist(data.frame(data[ , factor])),
                       position = "bottomright",
                       title = name,
                       labFormat = labelFormat(suffix = ""))   #labelFormat(suffix = "")
return(map1)
}


#---------------#
#-- PREP DATA --#
#---------------#
yr <- 2000
metric <- 'Tair_d90F'

counties <- readOGR("H:/Jason/GIS/Census/geometries", "NM_counties")   #open geometry
counties@data$COUNTYFP <- as.integer(as.character(counties@data$COUNTYFP))
#c_data  <- read.csv('L:/NLDAS/NCDMs/NM_county_temperature_NCDMs_1981-2015.csv', header = T)   #open tabular data
c_data   <- read.csv('H:/Jason/Climate/temperature_NCDMs/from_GHCN/results/NM_county_temperature_NCDMs_1981-2015.csv', header = T)   #open tabular data
c_data   <- c_data[(c_data$year == yr & c_data$metric == metric),]   #subset tabular data
counties <- geo_join(counties, c_data, "COUNTYFP", "COUNTYFP")   #merge geometry and tabular data

small_areas <- readOGR("H:/Jason/GIS/Political_boundaries/NM_small_areas", "Sarea134Detailed")
small_areas@data$SMALL_AREA <- as.integer(as.character(small_areas@data$OBJECTID))
s_data     <- read.csv('L:/NLDAS/NCDMs/NM_small_area_temperature_NCDMs_1981-2015.csv', header = T)
#s_data      <- read.csv('H:/Jason/Climate/temperature_NCDMs/from_GHCN/results/NM_small_area_temperature_NCDMs_1981-2015.csv', header = T)
s_data      <- s_data[(s_data$year == yr & s_data$metric == metric), ]
small_areas <- geo_join(small_areas, s_data, "SMALL_AREA", "SMALL_AREA")


#---------------#
#-- MAKE MAPS --#
#---------------#
make_map(data = small_areas, 
         factor = 'ndays',
         geo_name = 'SMALL_AREA',
         name = paste(metric, yr, sep = ' '))
