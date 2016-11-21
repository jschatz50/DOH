#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script downloads global summary of the day (GSOD) weather files,
# outputs the files as CSVs, and checks to see which downloads failed.


#------------------#
#--OPEN LIBRARIES--#
#------------------#
library(RCurl)


#--------------------#
#--DEFINE FUNCTIONS--#
#--------------------#
#returns objects that are not in another object
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


#---------------------------#
#-- generate list of URLs --#
#---------------------------#
stations <- read.csv("H:/Jason/Climate/GSOD/nm_gsod_stations.csv",header=T)

url_list <- list()
for (i in 1:nrow(stations)){
   code   <- paste(stations$USAF[i], sprintf("%05d", stations$WBAN[i]), sep = "-")
   years  <- seq(stations$start_year[i], stations$end_year[i], 1)
   list_a <- list()
   for (j in 1:length(years)){
      list_a[[j]] <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/gsod/", years[j], "/", code, "-", years[j], ".op.gz", sep = "")
   }
   url_list[[i]] <- list_a
}
urls <- do.call(c, unlist(url_list, recursive = FALSE))

#-------------------#
#-- download data --#
#-------------------#
for (i in 1:length(urls)){
   tryCatch({
      URL  <- urls[i] 
      dest <- paste("H:/Jason/Climate/GSOD/GZs/", substr(urls[i], 44, 66), sep = "")
      download.file(URL, dest, quiet = FALSE)
      Sys.sleep(1)
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#-------------------#
#-- process files --#
#-------------------#
files <- list.files(path = "H:/Jason/Climate/GSOD/GZs", full = T)
names <- list.files(path = "H:/Jason/Climate/GSOD/GZs")
names <- gsub(".op.gz", "", names)
column.widths <- c(6,1,5,2,8,4,4,1,2,4,4,1,2,2,6,1,2,3,5,1,2,3,4,1,2,3,4,1,2,3,4,2,5,3,5,4,4,2,6,1,5,2,6)

for (i in 1:length(files)){
   tryCatch({
      data <- read.fwf(paste(files[i], sep = ""), column.widths, skip = 1)
      data[ ,c(seq(2, 42, 2))] <- NULL
      names(data) <- c("STN---", "WBAN", "YEARMODA", "TEMP", "TCNT", "DEWP", "DPCNT", "SLP",
                      "SLPCNT", "STP", "STPCNT", "VISIB", "VISIBCNT", "WDSP", "WNDCNT",
                      "MXSPD", "GUST", "MAX", "MIN", "PRCP", "SNDP", "FRSHTT")
      outname <- paste("H:/Jason/Climate/GSOD/CSVs/", names[i], ".csv", sep = "")
      write.csv(data, outname, row.names = F)
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#--------------------------------------#
#-- which files didn't come through? --#
#--------------------------------------#
names   <- list.files(path = "H:/Jason/Climate/GSOD/CSVs")
names   <- gsub(".csv", "", names)
urls2   <- substr(urls, 44, 60)
missing <- urls2[urls2 %not in% names]
getURL("ftp://ftp.ncdc.noaa.gov/pub/data/gsod/",verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)


#-------------------------------#
#-- re-download missing files --#
#-------------------------------#
url_list2 <- list()
for (i in 1:length(missing)){
   year1 <- substr(missing[i], 14, 18)
   url_list2[[i]] <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/gsod/", year1, "/", missing[i], ".op.gz", sep = "")
}
urls2 <- unlist(url_list2)

for (i in 1:length(urls2)){
   tryCatch({
      URL  <- urls2[i] 
      dest <- paste("H:/Jason/Climate/GSOD/GZs/", substr(urls2[i], 44, 66), sep = "")
      Sys.sleep(1)
      download.file(URL, dest, quiet = FALSE)
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.csv(urls2,"remaining_urls.csv",row.names=F)
