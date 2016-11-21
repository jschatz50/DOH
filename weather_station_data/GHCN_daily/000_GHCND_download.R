#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script downloads GHCN daily weather files for NM and formats them for futher use


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(reshape2)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


#---------------------------#
#-- GENERATE LIST OF URLs --#
#---------------------------#
stations <- read.csv("H:/Jason/Climate/station_list.csv",header=T)

url_list <- list()
for (i in 1:nrow(stations)){
   url_list[[i]] <- paste("http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/",stations$STATION_ID[i], ".dly", sep="")
}

urls <- unlist(url_list)


#-------------------#
#-- download data --#
#-------------------#
for (i in 1:length(urls)){
   tryCatch({
      URL <- urls[i] 
      dest <- paste("H:/Jason/Climate/GHCND/", gsub("http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/", "", urls[i]), sep="")
      download.file(URL, dest, quiet=FALSE)
      Sys.sleep(1)	#so their server doesn't kick us off
   }, error <- function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#-------------------#
#-- process files --#
#-------------------#
files = list.files(path="H:/Jason/Climate/GHCND/downloads", full=T)
names = list.files(path="H:/Jason/Climate/GHCND/downloads")
names = gsub(".dly", "", names)
column.widths <- c(11, 4, 2, 4, rep(c(5, 1, 1, 1),31))

for (i in 1:length(files)){
   tryCatch({
      data <- read.fwf(paste(files[i], sep = ""), column.widths)
      data[data == -9999] <- NA
      data[ ,c(seq(6, 128, 4), seq(7, 128, 4), seq(8, 128, 4))] = NULL
      names(data) = c("SID", "YEAR", "MONTH", "VAR", seq(1, (ncol(data) - 4), 1))
      data = melt(data, id=c("SID", "YEAR", "MONTH", "VAR")); names(data)[5] = "DAY"
      data = dcast(data, SID + YEAR + MONTH + DAY ~ VAR)
      data = data[rowSums(is.na(data[c("TMAX", "TMIN")])) != 2, ]
      outname = paste("H:/Jason/Climate/GHCND/CSVs/", names[i], ".csv", sep="")
      write.csv(data, outname, row.names=F)
   }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#--------------------------------------#
#-- which files didn't come through? --#
#--------------------------------------#
names2 = list.files(path = "H:/Jason/Climate/GHCND/CSVs")
names2 = gsub(".csv", "", names2)

missing = names[names %not in% names2]









