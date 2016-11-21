#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script finds which stations are within ~1 mile of one another, which will allow
# us to merge stations that have stayed in the same location and changed station ID #s
# or only moved a short distance


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(spatstat)
library(sp)
library(rgeos)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
import_csvs = function(x){
   # Opens multiple CSVs
   # 
   # Args:
   #   x: path containing CSVs
   #
   # Returns:
   #   filenames:   list of import paths
   #   names:       filenames only
   #   import_list: list containing CSVs
   c(filenames   = list.files(path = x, full = T, pattern = "*.csv"),
     names       = list.files(path = x, full = F, pattern = "*.csv"),
     import_list = lapply(filenames, read_csv))
}

#move file from one directory to another
file.move <- function(from, to) {
   todir <- dirname(to)
   if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive = TRUE)
   file.rename(from = from,  to = to)
}


#------------------------#
#-- NEIGHBOR DISTANCES --#
#------------------------#
df1 <- read.csv(file.choose(), header = T)	#station locations
coordinates(df1) <- cbind(df1$LONGITUDE, df1$LATITUDE)

d <- gDistance(df1, byid=T)
d[d > 0.017] <- NA	#~within a mile
a <- unique(which(colSums(!is.na(d)) > 1))
df2 <- data.frame(df1)
df3 <- df2[a, ]
write.csv(df3, "matched_stations.csv", row.names = F)


#-----------------------------------------#
#--MOVE MATCHED STATIONS TO SEPARATE DIR--#
#-----------------------------------------#
matches = df3$STATION_ID

path1 = "H:/Jason/Climate/GHCND/CSVs"
import_csvs(path1)

names = gsub(".csv","",names)
names2 = names[names %in% df3$STATION_ID]
names2 = paste(names2,".csv",sep="")

for (i in 1:length(names2)){
   file.move(from = paste("H:/Jason/Climate/GHCND/CSVs/", names2[i], sep=""),
             to = paste("H:/Jason/Climate/GHCND/CSVs/individual_stations/", names2[i], sep=""))
}


#-------------------------------------------#
#-- MERGE THE MATCHED STATIONS (SLOW WAY) --#
#-------------------------------------------#
name1   <- "USC00298192.csv"
name2   <- "USC00299820.csv"
#names3 <- "USC00292324.csv"

df1 <- read.csv(name1, header = T)
df2 <- read.csv(name2, header = T)

df1 <- df1[,c("SID","YEAR","MONTH","DAY","TMAX","TMIN","PRCP")]
names(df1)[7] = "PPT"
df2 <- df2[,c("SID","YEAR","MONTH","DAY","TMAX","TMIN","PRCP")]
names(df2)[7] = "PPT"

final = rbind(df1, df2)
write.csv(final, paste(name1, name2, sep = "-"), row.names = F)


