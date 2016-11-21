#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 11/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script merges individual GSOD station files into a single file and
# performs basic quality control.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(data.table)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
F2C <- function(x){(x-32) * 5/9}   #degF to degC
C2F <- function(x){x * 9/5 + 32}   #degC to degF
RHfromTDP <- function(T,DP){100 * 
	                       (exp((17.625 * DP) / 
	                       (243.04 + DP)) / 
	                       exp((17.625 * T) / 
	                       (243.04 + T)))
                           }

import_csvs <- function(x){
   # Opens multiple CSVs
   # 
   # Args:
   #   x: path containing CSVs
   #
   # Returns:
   #   filenames:   list of import paths
   #   names:       filenames only
   #   import_list: list containing CSVs
   c(filenames   <- list.files(path = x, full = T, pattern = "*.csv"),
     names       <- list.files(path = x, full = F, pattern = "*.csv"),
     import_list <- lapply(filenames, read.csv, header=T))
}


#----------------------------#
#-- merge files by station --#
#----------------------------#
# path1 <- "H:/Jason/Climate/GSOD/CSVs"
# import_csvs(path1)
# names <- gsub(".csv", "", names)
# station.levs <- levels(factor(substr(names, 1, 12)))
# stations     <- substr(names, 1, 12)
#
# for (i in 1:length(station.levs)){
#    tryCatch({
#       stns    <- which(stations == station.levs[i])
#       df.list <- import.list[stns]
#       merged  <- rbindlist(df.list, fill = T)
#       outname <- paste("H:/Jason/Climate/GSOD/merged/", station.levs[i], ".csv", sep = "")
#       write.csv(merged, outname, row.names = F)
#    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }


#-----------------------------------------------------#
#-- on second thought, merge everything to one file --#
#-----------------------------------------------------#
path1 <- "H:/Jason/Climate/GSOD/CSVs"
import_csvs(path1)
names <- gsub(".csv","",names)
merged <- rbindlist(import.list,fill=T)
names(merged) <- c("STN", "WBAN", "DATE", "TAVG", "TCNT", "DP", "DPCNT", "SLP", "SLPCNT", 
	               "STP", "STPCNT", "VISIB", "VISIBCNT", "WINDSPD", "WINDCNT", "MAXWIND",
	               "MAXGUST", "TMAX", "TMIN", "PPT", "SNOWDEPTH", "FRSHTT")
merged$SID <- paste(merged$STN, sprintf("%05d", merged$WBAN), sep = "-")

##replace missing vals with NA (all different, for some reason)
merged$DP[merged$DP == 99.9] <- NA
merged$SLP[merged$SLP == 9999.9] <- NA
merged$STP[merged$STP == 999.9] <- NA
merged$VISIB[merged$VISIB == 99.9] <- NA
merged$WINDSPD[merged$WINDSPD == 99.9] <- NA
merged$MAXWIND[merged$MAXWIND == 99.9] <- NA
merged$MAXGUST[merged$MAXGUST == 999.9] <- NA
merged$TMAX[merged$TMAX == 999.9] <- NA
merged$TMIN[merged$TMIN == 99.9] <- NA
merged$SNOWDEPTH[merged$SNOWDEPTH == 999.9] <- NA
merged$PPT <- as.numeric(substr(as.character(merged$PPT), 1, 5))
merged$PPT[merged$PPT == 99.99] <- NA

##covert knots to mph
merged$WINDSPD <- merged$WINDSPD * 1.15077944802
merged$MAXWIND <- merged$MAXWIND * 1.15077944802
merged$MAXGUST <- merged$MAXGUST * 1.15077944802

##filter by number of daily obs (at least 20 of 24 hours for most)
merged[(merged$TCNT     < 20), "TAVG"]    <- NA
merged[(merged$TCNT     < 20), "TMAX"]    <- NA
merged[(merged$TCNT     < 20), "TMIN"]    <- NA
merged[(merged$DPCNT    < 20), "DP"]      <- NA
merged[(merged$SLPCNT   < 12), "SLP"]     <- NA
merged[(merged$STPCNT   < 12), "STP"]     <- NA
merged[(merged$VISIBCNT < 12), "VISIB"]   <- NA
merged[(merged$WINDCNT  < 20), "WINDSPD"] <- NA

##filter 
merged[(merged$TMIN >   90), 'TMIN']  <- NA
merged[(merged$TMAX < (-20)), 'TMAX'] <- NA

##calculate mean daily RH
merged$RH <- RHfromTDP(T = F2C(merged$TAVG), DP = F2C(merged$DP))   #units should be in degC
merged[(merged$RH>100), 'RH'] <- NA
#merged$RH = 100 * (exp((17.625*F2C(merged$DP))/(243.04+F2C(merged$DP)))/exp((17.625*F2C(merged$TAVG))/(243.04+F2C(merged$TAVG))))

##write results
write.csv(merged, "NM_GSOD_data_1941-2016.csv", row.names = F)

