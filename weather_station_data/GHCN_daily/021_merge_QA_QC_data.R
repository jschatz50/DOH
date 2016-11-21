#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script emrges the QA/QCd data


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(tidyverse)


#----------------------#
#-- READ/MERGE/WRITE --#
#----------------------#
tmins  <- read_csv("H:/Jason/Climate/GHCND/cleaned_tmins.csv")
tmaxes <- read_csv("H:/Jason/Climate/GHCND/cleaned_tmaxes.csv")
ppt    <- read_csv("H:/Jason/Climate/GHCND/cleaned_ppt_trange.csv")

tmins$row  <- NULL
tmaxes$row <- NULL

tmins  <- melt(tmins, id = c("YEAR", "MONTH", "DAY"))
tmaxes <- melt(tmaxes, id = c("YEAR", "MONTH", "DAY"))
names(tmins)  <- c("YEAR", "MONTH", "DAY", "SID", "TMIN")
names(tmaxes) <- c("YEAR", "MONTH", "DAY", "SID", "TMAX")

temps <-merge(tmins, tmaxes, all = T, by = c("YEAR", "MONTH", "DAY", "SID"))
all   <- merge(temps, ppt, all = T, by = c("SID", "YEAR", "MONTH", "DAY"))
final <- all[apply(all[,5:8], 1, function(x) sum(!is.na(x))>0),]

final$TMIN   <- as.numeric(final$TMIN)
final$TMAX   <- as.numeric(final$TMAX)
final$TRANGE <- final$TMAX - final$TMIN
final$TRANGE <- NULL

write.csv(final, "ghcnd_data_1892-2016.csv", row.names = F)