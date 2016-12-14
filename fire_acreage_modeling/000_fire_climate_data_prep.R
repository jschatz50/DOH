#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# created:       11/10/16
# last modified: 12/13/16


#---------------#
#--DESCRIPTION--#
#---------------#
# This script:
# (1) calculates 1981-2010 climate normals for precipitation and temperature
#     for different seasons of the year.
# (2) using those normals, calculates seasonal and annual deviations from 
#     normal for temperature and precipitation at MET stations across New 
#     Mexico for 1980-2015.
# (2) merges climate and fire data for analysis.


#--------------------#
#-- OPEN LIBRARIES --#
#--------------------#
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(randomForest)
library(rpart.plot)
library(vita)
library(MASS)
library(data.table)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


#--------------------#
#-- PREP FIRE DATA --#
#--------------------#
fire <- read.csv("H:/Jason/Wildfire/NM_fires_1980-2015.csv", header = T)
fire <- with(fire, aggregate(x = TOTALACRES, by = list(YEAR), FUN = sum))
names(fire) <- c("YEAR", "ACRES")


#-----------------------#
#-- PREP CLIMATE DATA --#
#-----------------------#
# Calculate the median statewide deviation from 1981-2010 normals for precip 
# and daily max temperature from ~63 stations in NM (these 63 stations have
# 95% data coverage for this period).  This section returns both seasonal 
# (winter, spring, summer, fall) and annual statewide deviations from normal 
# for each year.

## open climate data
clim <- read.csv("H:/Jason/Climate/past_climate/GHCND/nm_ghcnd_data_1892-2016.csv", header = T)
clim <- subset(clim, c(YEAR >= 1980 & YEAR <= 2015))
clim$TAVG = (clim$TMIN + clim$TMAX) / 2

## pull stations with good coverage
goods <- with(clim, aggregate(x = TMAX, by = list(SID), FUN = length))
names(goods) = c("SID", "TMAX")
goods$percent = goods$TMAX / 13149 * 100
goods = goods[goods$percent >= 95, 'SID']   #stations w/>=95% data coverage
clim = clim[clim$SID %in% goods, ]

## calculate seasonal values by year
# associate month 12 with the following year (for seasonal continuity)
clim$YEAR <- with(clim, ifelse(MONTH == 12, YEAR + 1, YEAR))

clim$season <- with(clim, ifelse(c(MONTH == 12 | MONTH == 1 | MONTH == 2), "winter",
                          ifelse(c(MONTH ==  3 | MONTH == 4 | MONTH == 5), "spring",
                          ifelse(c(MONTH ==  6 | MONTH == 7 | MONTH == 8), "summer",
                          "fall"))))

## filter out incomplete seasons
t_seasonal <- aggregate(data = clim, TMAX ~ YEAR + SID + season, mean)
p_seasonal <- aggregate(data = clim, PPT  ~ YEAR + SID + season, sum)
lengths   <- aggregate(data = clim, cbind(TMAX, PPT) ~ YEAR + SID + season, length)
names(lengths)[4:5] <- c("TCNT", "PCNT")
df_list <- list(t_seasonal, p_seasonal, lengths)
seasonal <- Reduce(function(...) merge(..., id = c("SID", "YEAR", "season"), all=T), df_list)

seasonal$TMAX2 <- ifelse(seasonal$TCNT <= 75, seasonal$TMAX2 <- NA, seasonal$TMAX2 <- seasonal$TMAX)
seasonal$TMAX  <- seasonal$TMAX2; seasonal$TMAX2 <- NULL
seasonal$PPT2  <- ifelse(seasonal$PCNT <= 75, seasonal$PPT2 <- NA, seasonal$PPT2 <- seasonal$PPT)
seasonal$PPT   <- seasonal$PPT2; seasonal$PPT2 <- NULL
seasonal$TCNT  <- NULL
seasonal$PCNT  <- NULL

## calculate 1981-2010 seasonal normals at each MET station
TMAX_normals <- aggregate(data = subset(seasonal, c(YEAR >= 1981 & YEAR <= 2010)),
                          TMAX ~ SID + season, 
                          FUN = mean)
PPT_normals  <- aggregate(data = subset(seasonal, c(YEAR >= 1981 & YEAR <= 2010)),
                          PPT ~ SID + season, 
                          FUN = sum)
normals <- merge(TMAX_normals, PPT_normals, all = T, id = c("SID", "season"))
normals <- with(normals, normals[order(SID, season), ])
names(normals)[3:4] <- c("TMAX_norm", "PPT_norm")
normals$PPT_norm <- normals$PPT_norm / 30   #to get per year (sums over all 30 years, originally)

## merge observed data and normals
merged <- merge(seasonal, normals, all.x = T, id = c("SID", "season"))
merged$TMAX_diffnorm <- merged$TMAX - merged$TMAX_norm
merged$PPT_diffnorm  <- merged$PPT  - merged$PPT_norm

statewide <- aggregate(data = merged, cbind(TMAX_diffnorm, PPT_diffnorm) ~ SID + YEAR + season, mean)
statewide$PPT_diffnorm <- statewide$PPT_diffnorm * 3   #scale PPT mean up to PPT total
statewide <- aggregate(data = statewide, cbind(TMAX_diffnorm, PPT_diffnorm) ~ YEAR + season, median)
annual   <- aggregate(data = merged, cbind(TMAX_diffnorm, PPT_diffnorm) ~ YEAR, median)

statewide <- dcast(setDT(statewide), YEAR ~ season, value.var=c("TMAX_diffnorm", "PPT_diffnorm"))
statewide <- statewide[statewide$YEAR != 2016, ]
colnames(statewide)[1] = "YEAR"

## finalize climate data
final_clim = merge(annual, statewide)
final_clim$PPT_diffnorm_spring_summer = final_clim$PPT_diffnorm_spring + final_clim$PPT_diffnorm_summer
final_clim$TMAX_diffnorm_spring_summer = (final_clim$TMAX_diffnorm_spring + final_clim$TMAX_diffnorm_summer) / 2 

final = merge(fire, final_clim)
write.csv(final, 'H:/Jason/Wildfire/merged_climate_fire_data_1980-2015.csv', row.names = F)

