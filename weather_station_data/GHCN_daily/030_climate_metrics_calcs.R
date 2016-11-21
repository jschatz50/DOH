#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script calculates desired climate metrics (e.g. 'average' number of days over 90F)


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(reshape2)
library(sp)


#---------------#
#-- read data --#
#---------------#
stations <- read.csv("H:/Jason/Climate/ghcn_station_list.csv",header <- T)
stations <- stations[, c('STATION_ID', 'LATITUDE', 'LONGITUDE', 'ELEVATION')]
tmaxes   <- read.csv("H:/Jason/Climate/GHCND/cleaned_tmaxes.csv", header = T)
tmins    <- read.csv("H:/Jason/Climate/GHCND/cleaned_tmins.csv", header = T)


#------------------#
#-- calculations --#
#------------------#

#--define starting and ending years--#
start_year <- 1981
end_year   <- 2010

#--prep data--#
tmaxes <- subset(tmaxes, YEAR >= start_year & YEAR <= end_year)
tmins  <- subset(tmins, YEAR >= start_year & YEAR <= end_year)
tmaxes <- melt(tmaxes, id = c("YEAR", "MONTH", "DAY"))
tmins  <- melt(tmins,  id = c("YEAR", "MONTH", "DAY"))

per_cov_max <- with(tmaxes, aggregate(x = value, by = list(variable), FUN = function(x) {sum(!is.na(x)) / 10957}))   #% coverage for Tmax
per_cov_min <- with(tmins,  aggregate(x = value, by = list(variable), FUN = function(x) {sum(!is.na(x)) / 10957}))   #% coverage for Tmin

completeness  <- 0.75	##according to http://journals.ametsoc.org/doi/full/10.1175/BAMS-D-11-00197.1, to compute NCDC normals, at least 10 of 30 years of each month must be complete)
tmax_stations <- per_cov_max[(per_cov_max$x >= completeness), 1]
tmin_stations <- per_cov_min[(per_cov_min$x >= completeness), 1]

tmaxes_good <- tmaxes[tmaxes$variable %in% tmax_stations,]
tmins_good  <- tmins[tmins$variable %in% tmin_stations,]

#--calculate temperature thresholds--#
#Tmax (hot days)
max_thresholds <- with(tmaxes_good, aggregate(x = value, by = list(variable), 
		               FUN   = function(x){
		               	  c(d90F  <- (sum(x >= 90,  na.rm = T) / 30) * sum(!is.na(x)) / 10957,
					        d95F  <- (sum(x >= 95,  na.rm = T) / 30) * sum(!is.na(x)) / 10957,
					        d100F <- (sum(x >= 100, na.rm = T) / 30) * sum(!is.na(x)) / 10957)
		               }))
max_thresholds <- cbind(max_thresholds$Group.1, data.frame(max_thresholds$x))
names(max_thresholds) <- c("STATION_ID", "d90F", "d95F", "d100F")
max_thresholds <- merge(stations, max_thresholds, all.y = T)

#Tmin (hot nights)
min_thresholds <- with(tmins_good, aggregate(x = value, by = list(variable), 
		               FUN = function(x){n80F = (sum(x >= 80, na.rm = T) / 30) * sum(!is.na(x)) / 10957}))
min_thresholds <- cbind(min_thresholds$Group.1, data.frame(min_thresholds$x))
names(min_thresholds) <- c("STATION_ID", "n80F")

#--merge files and write results to file--#
final <- merge(max_thresholds, min_thresholds, all.x = T)
write.csv(final, "NM_high_temperature_thresholds_1981-2010.csv", row.names=F)
