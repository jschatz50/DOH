#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 11/10/16


#---------------#
#--DESCRIPTION--#
#---------------#
# This script loads fire history data for New Mexico from 1980-2015, merges it with climate data
# from the same time period, and test for correlations between burned acres and weather.  The
# objective is to find a model that can be applied to future climate projections as a first order
# estimate of how climate change may affect fires in New Mexico.


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
# This script calculates the statewide median deviation from 1981-2010 normals for
# precip and daily average temperature from ~63 stations (those with >95% data coverage.
# It returns both seasonal (winter, spring, summer, fall) and annual statewide deviations 
# from normal for each year.

clim <- read.csv("H:/Jason/Climate/past_climate/GHCND/nm_ghcnd_data_1892-2016.csv", header = T)
clim <- subset(clim, c(YEAR >= 1980 & YEAR <= 2015))
clim$TAVG = (clim$TMIN + clim$TMAX) / 2

#--pull stations with good coverage--#
goods <- with(clim, aggregate(x = TAVG, by = list(SID), FUN = length))
names(goods) = c("SID", "TAVG")
goods$percent = goods$TAVG / 13149 * 100
goods = goods[goods$percent >= 95, 'SID']   #stations w/>=95% data coverage
clim = clim[clim$SID %in% goods, ]

#--calculate seasonal values by year--#
#--associate month 12 with the following year (for seasonal continuity)--#
clim$YEAR <- with(clim, ifelse(MONTH == 12, YEAR + 1, YEAR))

clim$season <- with(clim, ifelse(c(MONTH == 12 | MONTH == 1 | MONTH == 2), "winter",
                          ifelse(c(MONTH ==  3 | MONTH == 4 | MONTH == 5), "spring",
                          ifelse(c(MONTH ==  6 | MONTH == 7 | MONTH == 8), "summer",
                          "fall"))))

#filter out incomplete seasons
t_seasonal <- aggregate(data = clim, TAVG ~ YEAR + SID + season, mean)
p_seasonal <- aggregate(data = clim, PPT  ~ YEAR + SID + season, sum)
lengths   <- aggregate(data = clim, cbind(TAVG, PPT) ~ YEAR + SID + season, length)
names(lengths)[4:5] <- c("TCNT", "PCNT")
df_list <- list(t_seasonal, p_seasonal, lengths)
seasonal <- Reduce(function(...) merge(..., id = c("SID", "YEAR", "season"), all=T), df_list)

seasonal$TAVG2 <- ifelse(seasonal$TCNT <= 75, seasonal$TAVG2 <- NA, seasonal$TAVG2 <- seasonal$TAVG)
seasonal$TAVG  <- seasonal$TAVG2; seasonal$TAVG2 <- NULL
seasonal$PPT2  <- ifelse(seasonal$PCNT <= 75, seasonal$PPT2 <- NA, seasonal$PPT2 <- seasonal$PPT)
seasonal$PPT   <- seasonal$PPT2; seasonal$PPT2 <- NULL
seasonal$TCNT  <- NULL
seasonal$PCNT  <- NULL

#--calculate 1981-2010 seasonal normals--#
TAVG_normals <- aggregate(data = subset(seasonal, c(YEAR >= 1981 & YEAR <= 2010)),
                          TAVG ~ SID + season, 
                          FUN = mean)
PPT_normals  <- aggregate(data = subset(seasonal, c(YEAR >= 1981 & YEAR <= 2010)),
                          PPT ~ SID + season, 
                          FUN = sum)
normals <- merge(TAVG_normals, PPT_normals, all = T, id = c("SID", "season"))
normals <- with(normals, normals[order(SID, season), ])
names(normals)[3:4] <- c("TAVG_norm", "PPT_norm")
normals$PPT_norm <- normals$PPT_norm / 30   #to get per year (sums over all 30 years, originally)

#--merge observed and normals--#
merged <- merge(seasonal, normals, all.x = T, id = c("SID", "season"))
merged$TAVG_diffnorm <- merged$TAVG - merged$TAVG_norm
merged$PPT_diffnorm  <- merged$PPT  - merged$PPT_norm

statewide <- aggregate(data = merged, cbind(TAVG_diffnorm, PPT_diffnorm) ~ SID + YEAR + season, mean)
statewide$PPT_diffnorm <- statewide$PPT_diffnorm * 3   #scale PPT mean up to PPT total
statewide <- aggregate(data = statewide, cbind(TAVG_diffnorm, PPT_diffnorm) ~ YEAR + season, median)
annual   <- aggregate(data = merged, cbind(TAVG_diffnorm, PPT_diffnorm) ~ YEAR, median)

require(data.table) 
statewide <- dcast(setDT(statewide), YEAR ~ season, value.var=c("TAVG_diffnorm", "PPT_diffnorm"))
statewide <- statewide[statewide$YEAR != 2016, ]
colnames(statewide)[1] = "YEAR"

final_clim = merge(annual, statewide)
final_clim$PPT_diffnorm_spring_summer = final_clim$PPT_diffnorm_spring + final_clim$PPT_diffnorm_summer
final_clim$TAVG_diffnorm_spring_summer = (final_clim$TAVG_diffnorm_spring + final_clim$TAVG_diffnorm_summer) / 2 


#------------------#
#-- LINEAR MODEL --#
#------------------#
final = merge(fire, final_clim)
#final = final[-32,]

# first differencing
df1 = data.frame(matrix(nrow = (nrow(final)-1), ncol = 13))
for (i in 2:14){
   df1[,(i-1)] <- diff(final[,i])
}
names(df1) = names(final)[2:14]

#--test predictive ability of model--#
n = 10000

results1 = data.frame(matrix(nrow = n, ncol = 2))
for (i in 1:n){
   training <- df1[sample(nrow(df1), 25), ]
   test     <- df1[df1$ACRES %not in% training$ACRES, ]

   rlm1 <- rlm(ACRES ~ PPT_diffnorm_spring_summer + TAVG_diffnorm_spring_summer, training)

   predicted  <- predict(rlm1, test)
   observed <- test$ACRES
   results1[i,2] <- cor(predicted, observed)^2
   results1[i,1] <- sqrt(mean(predicted-observed)^2)
}
names(results1) = c("RMSE","RSQUARED")
summary(results1)

#--fit full model--#
summary(rlm1 <- rlm(ACRES ~ PPT_diffnorm_spring_summer + TAVG_diffnorm_spring_summer, df1))
weights1 <- data.frame(year = final$YEAR[-1], resid = rlm1$resid, weight = rlm1$w)
weights2 <- weights1[order(rlm1$w), ]

plot(predict(rlm1, df1), type = 'l', col = 'red', ylim = range(predict(rlm1,df1), df1$ACRES))
lines(df1$ACRES)

cor(predict(rlm1, df1), df1$ACRES)^2


#-------------------#
#-- RANDOM FOREST --#
#-------------------#
final = merge(fire, final_clim)
final = final[-32,]   #excluding 2011, which is an enormous outlier

# first differencing
df1 = data.frame(matrix(nrow = (nrow(final)-1), ncol = 13))
for (i in 2:14){
   df1[,(i-1)] <- diff(final[,i])
}
names(df1) = names(final)[2:14]

# number of combinations of years for training/test data 
n = 1000

results = data.frame(matrix(nrow = n, ncol = 3))
for (i in 1:n){
   training <- df1[sample(nrow(df1), 25), ]
   test     <- df1[df1$ACRES %not in% training$ACRES, ]
   formula1 <- formula(ACRES ~ PPT_diffnorm_spring_summer + TAVG_diffnorm_spring_summer)
   rf <- randomForest(formula1, data = training, importance = TRUE, ntree = 1000, mtry = 1)
   pred <- predict(rf)
   obs  <- training$ACRES
   results[i,1] <- cor(pred, obs)^2   
   predicted <- predict(rf, test)
   observed  <- test$ACRES
   results[i,2] <- cor(predicted,observed)^2
   results[i,3] <- sqrt(mean(predicted-observed)^2)
}

names(results) <- c("training_rsq", "test_rsq", "test_rmse")
(mean(results$training_rsq))
(mean(results$test_rsq))
(mean(results$test_rmse))

#  Xs <- training[,c("PPT_diffnorm_summer", "PPT_diffnorm_spring", 
#                    "TAVG_diffnorm_spring")]
#  Ys <- training[,c("ACRES")]
#  pimp.varImp.reg <- PIMP(Xs, Ys, rf, S=100, ncores=3)
#  (pimp.t.reg <- PimpTest(pimp.varImp.reg))
#  varImpPlot(rf)
partialPlot(rf, training, PPT_diffnorm_spring_summer)
plot(observed,type='l',ylim=range(predicted,observed),xaxt='n')
     axis(1, at=1:nrow(test), labels=test$MONTH)
lines(predicted,col='red')

(rf <- randomForest(formula1, data = df1, importance = TRUE, ntree = 1000, mtry = 1))
pred <- predict(rf)
obs  <- df1$ACRES
plot(pred,type='l')
lines(obs,col='red')



#---------------#
#-- leftovers --#
#---------------#

#--detrend data using linear time trend (detrended residuals + 1981-2010 mean)--#

final = merge(fire, final_clim)
df1 = data.frame(matrix(nrow = nrow(final), ncol = 13))   #empty dataframe
for (i in 2:14){
   lm1 <- lm(final[ , i] ~ final[ , 1])
   df1[,(i-1)] <- lm1$residuals + mean(final[2:31, i])
}
names(df1) <- names(final)[2:14]


#--from first differences to original values--#

x        <- data.frame(rnorm(1:100)); names(x) = "x"
x_diff   <- data.frame(diff(x[ , 1])); names(x_diff) = "x"
init     <- x[1,]   #first value...could do with any value

list_a = list()
list_a[[1]] = init
for (i in 1:nrow(x_diff)){
   list_a[[i+1]] <- list_a[[i]] + x_diff$x[i]
}

rev   <- unlist(list_a); names(res) <- NULL
x$rev <- rev


#--ols regression with log response variable--#

lm1 <- lm(log(ACRES) ~ PPT_diffnorm_spring_summer + TAVG_diffnorm_spring_summer, df1); summary(lm1)
100 * (exp(lm1$coefficients[2]) -1)   #percent change in arithmetic mean per unit change in PPT
100 * (exp(lm1$coefficients[3]) -1)   #percent change in arithmetic mean per unit change in TAVG
toy1 <- data.frame(PPT_diffnorm_spring_summer  = c(-5,-2,-1,0,0,1,2,5),
          TAVG_diffnorm_spring_summer = c(5,2,1,0,0,-1,-2,-5))
predict(rlm1, toy1)



