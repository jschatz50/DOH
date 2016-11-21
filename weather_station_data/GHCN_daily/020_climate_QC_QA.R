#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script performs QA/QC on the GHCND data.
# First, impossible values are removed.
# Second, for each month/sensor, locations whose data diverge significantly from
# statewide trends are eliminated
# Third, for each sensor, dates where data diverge significantly from the average
# (over time) relationship with statewide weather are eliminated.  See README for more.


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(data.table)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
C2F <- function(x){x*9/5+32}  #degC to degF
mm2in <- function(x){x/25.4}  #mm to inches
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


#--------------------------------------#
#-- REDUCE DFS TO RELEVANT VARIABLES --#
#--------------------------------------#
path1 <- "H:/Jason/Climate/GHCND/CSVs/individual_stations"
import_csvs(path1)

list_a <- list()
for (i in 1:length(import_list)){
   tryCatch({
      df1 <- import_list[[i]]
      df1 <- df1[ ,c("SID", "YEAR", "MONTH", "DAY", "TMAX", "TMIN", "PPT")]
	  names(df1)[7] <- "PPT"
	  list_a[[i]] <- data.frame(cbind(names[i],nrow(df1)/365))
	  write.csv(df1, paste("H:/Jason/Climate/GHCND/CSVs/processed/", names[i], sep = ""), row.names = F)
   }, error = function(e){cat("ERROR :",conditionMessage(e), "/n")})
}

num_years <- rbindlist(list_a)
names(num_years) <- c("file", "years")
write.csv(num_years, "number_of_years_per_station2.csv", row.names = F)


#---------------------------#
#-- PERCENT DATA COVERAGE --#
#---------------------------#
path1 <- "H:/Jason/Climate/GHCND/CSVs/processed"
import_csvs(path1)

list_a <- list()
for (i in 1:length(import_list)){
   df1        <- import_list[[i]]
   first_yr   <- min(df1$YEAR)
   last_yr    <- max(df1$YEAR)
   n_years    <- length(levels(factor(df1$YEAR)))
   possible_days <- n_years * 365.25
   tmax_pcov  <- round((sum(!is.na(df1$TMAX)) / possible_days) * 100, 1)
   tmin_pcov  <- round((sum(!is.na(df1$TMIN)) / possible_days) * 100, 1)
   ppt_pcov   <- round((sum(!is.na(df1$PPT)) / possible_days) * 100, 1)
   name       <- gsub(".csv", "", names[i])
   list_a[[i]] <- data.frame(SID = name, 
   							first_yr = first_yr, 
   							last_yr = last_yr, 
   							n_years = n_years, 
   							tmax_pcov = tmax_pcov, 
   							tmin_pcov = tmin_pcov, 
   							ppt_pcov =ppt_pcov)
}
final <- rbindlist(list_a)
write.csv(final, "percent_coverage_by_station.csv", row.names = F)


#---------------------------------------------------------------------#
#--					       			QA/QC	   		       		    --#
#-- remove impossible, poorly correlated, and highly unusual values --#
#---------------------------------------------------------------------#

df    <- read.csv("H:/Jason/Climate/GHCND/CSVs/percent_coverage_by_station.csv", header = T)   #percent coverage by station
df    <- subset(df, keep.=='y')
names <- paste(unlist(df$SID), ".csv", sep = "")

list_a <- list()
for (i in 1:length(names)){
   list_a[[i]] = read.csv(paste("H:/Jason/Climate/GHCND/CSVs/processed/", names[i], sep = ""), header = T)
}

final <- rbindlist(list_a)

#--covert to english units--#
final$TMAX   <- C2F(final$TMAX / 10)
final$TMIN   <- C2F(final$TMIN / 10)
final$PPT    <- mm2in(final$PPT / 10)
final$TRANGE <- abs(final$TMAX - final$TMIN)

#--remove impossible values--#
final[final$PPT>12, 'PPT'] = NA
final[(final$TMAX > 122 | final$TMAX < (-30)), 'TMAX']   <- NA
final[(final$TMIN > 100 | final$TMIN < (-50)), 'TMIN']   <- NA
final[(final$TRANGE < 1 | final$TRANGE > 70),  'TRANGE'] <- NA

#--for each month, correlate TMIN/TMAX with daily median of all sensors--#
#--Tmin--#

tmins <- data.frame(dcast(final, YEAR + MONTH + DAY ~ SID, value.var = "TMIN", fun.aggregate = mean))
tmins$median <- apply(tmins[ , 4:146], 1, median, na.rm = T)
tmins$year_month <- paste(tmins$YEAR, tmins$MONTH, sep = "_")

years <- levels(factor(tmins$year_month))
list_b <- list()
for (i in 1:length(years)){
   df1 <-= subset(tmins, year_month == years[i])
      for (j in 4:146){
	     tryCatch({
            corr=cor(df1[ , j], df1$median, use = "complete.obs")
            df1[,j] <- if (corr < 0.4) {df1[ , j] = NA} else {df1[ , j] = df1[ , j]}
         }, error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
      }
   list_b[[i]] <- df1
}

result <- data.frame(rbindlist(list_b))
result$row <- seq(1, nrow(result), 1)

fordata <- result[ , c("YEAR", "MONTH", "DAY", "row")]

#for all years, remove residuals >20F from relationship with median
for (i in 4:146){
   a <- result[complete.cases(result[,i]), c(i, 150)]
   a[which(abs(lm(result[,i] ~ result$median)$residuals) > 20), 1] <- NA
   fordata <- merge(fordata, a, all.x = T)
}

write.csv(fordata,"cleaned_tmins.csv",row.names=F)


#--Tmax--#

tmaxes            <- data.frame(dcast(final, YEAR + MONTH + DAY ~ SID, value.var = "TMAX", fun.aggregate = mean))
tmaxes$median     <- apply(tmaxes[,4:146], 1, median, na.rm = T)
tmaxes$year_month <- paste(tmaxes$YEAR, tmaxes$MONTH, sep = "_")

years  <- levels(factor(tmaxes$year_month))
list_b <- list()
for (i in 1:length(years)){
   df1 <- subset(tmaxes, year_month == years[i])
   for (j in 4:146){
      tryCatch({
         corr <- cor(df1[,j],df1$median,use = "complete.obs")
         df1[,j] <- if (corr < 0.4) {df1[,j] = NA} else {df1[,j] = df1[,j]}
      }, error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
   }
   list_b[[i]] <- df1
}

result <- data.frame(rbindlist(list_b))
result$row <- seq(1,nrow(result),1)
fordata <- result[,c("YEAR","MONTH","DAY","row")]

#for all years, remove residuals >20F from relationship with median
for (i in 4:146){
	a <- result[complete.cases(result[,i]), c(i, 150)]
	a[which(abs(lm(result[,i] ~ result$median)$residuals) > 20),1] <- NA
	fordata <- merge(fordata,a, all.x = T)
}

write.csv(fordata, "cleaned_tmaxes.csv", row.names = F)

####
x = data.frame(final)[,c('SID', 'YEAR', 'MONTH', 'DAY', 'PPT', 'TRANGE')]
write.csv(x, "cleaned_ppt_trange.csv", row.names=F)

