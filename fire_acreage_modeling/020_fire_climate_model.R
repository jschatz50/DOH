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
# (1) loads climate and fire history data for New Mexico from 1980-2015;
# (2) tests various models of the relationship between burned 
#     acres and weather;
# (3) uses those models to predict future burn acreage given
#     downscaled climate data (LOCA) for New Mexico.


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
library(scales)
library(ggplot2)
library(data.table)


#----------------------#
#-- DEFINE FUNCTIONS --#
#----------------------#
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


#---------------#
#-- READ DATA --#
#---------------#
historical <- read.csv("H:/Jason/Wildfire/merged_climate_fire_data_1980-2015.csv", header = T)
future = read.csv('H:/Jason/Wildfire/future_climate_data_for_model/statewide_medians/results.csv', header=T)


#-------------------------#
#-- ROBUST LINEAR MODEL --#
#-------------------------#
## detrend data
df1 = data.frame(matrix(nrow = nrow(historical), ncol = 13))   #empty dataframe
for (i in 2:14){
   tryCatch({
      lm1 <- lm(historical[ , i] ~ historical[ , 1])
      df1[,(i-1)] <- lm1$residuals + mean(historical[2:31, i])
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
names(df1) <- names(historical)[2:14]

## test predictive ability of model (10,000 random draws of years
## for the training/test set, since there are only 36 years)
n = 10000
results1 = data.frame(matrix(nrow = n, ncol = 2))
for (i in 1:n){
   training <- df1[sample(nrow(df1), 25), ]
   test     <- df1[df1$ACRES %not in% training$ACRES, ]

   rlm1 <- rlm(ACRES ~ PPT_diffnorm_spring_summer + TMAX_diffnorm_spring_summer, training)

   predicted  <- predict(rlm1, test)
   observed <- test$ACRES
   results1[i,2] <- cor(predicted, observed)^2
   results1[i,1] <- sqrt(mean(predicted-observed)^2)
}
names(results1) = c("RMSE","RSQUARED")
summary(results1)

## fit full model
summary(rlm1 <- rlm(ACRES ~ PPT_diffnorm_spring_summer + TMAX_diffnorm_spring_summer, df1))
weights1 <- data.frame(year = historical$YEAR[], resid = rlm1$resid, weight = rlm1$w)
weights2 <- weights1[order(rlm1$w), ]   # weights assigned to each observation

plot(predict(rlm1, df1), type = 'l', col = 'red', ylim = range(predict(rlm1,df1), df1$ACRES))
lines(df1$ACRES)
cor(predict(rlm1, df1), df1$ACRES)^2


#----------------------------------------------------#
#-- predict fire acreage using climate projections --#
#----------------------------------------------------#
predictions = predict(rlm1, future)
future$predicted_acreage = predictions
write.csv(future, 'acreage_predictions.csv', row.names = F)

## boxplots
future$PERIOD = factor(future$PERIOD, c("1986 to 2005", "2040 to 2059", "2080 to 2099"))
ggplot(aes(y = predicted_acreage, x = PERIOD, fill = SCENARIO), data = future) + 
           geom_boxplot() + 
           theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 panel.border = element_rect(colour = "black", fill=NA, size=2)) +
           scale_y_continuous(labels = comma) + 
           ylab('predicted acres burned per year')
           

#-------------------#
#-- RANDOM FOREST --#
#-------------------#
historical = historical[-32,]   #excluding 2011, which is an enormous outlier

## detrend
df1 = data.frame(matrix(nrow = nrow(historical), ncol = 13))   #empty dataframe
for (i in 2:14){
   tryCatch({
      lm1 <- lm(historical[ , i] ~ historical[ , 1])
      df1[,(i-1)] <- lm1$residuals + mean(historical[2:31, i])
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
names(df1) <- names(historical)[2:14]

## number of combinations of years for training/test data 
n = 1000

results = data.frame(matrix(nrow = n, ncol = 3))
for (i in 1:n){
   training <- df1[sample(nrow(df1), 25), ]
   test     <- df1[df1$ACRES %not in% training$ACRES, ]
   formula1 <- formula(ACRES ~ PPT_diffnorm_spring_summer + TMAX_diffnorm_spring_summer)
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
partialPlot(rf, training, TMAX_diffnorm_spring_summer)
plot(observed, type = 'l', ylim = range(predicted,observed), xaxt = 'n')
     axis(1, at = 1:nrow(test), labels = test$MONTH)
lines(predicted, col = 'red')

(rf <- randomForest(formula1, data = df1, importance = TRUE, ntree = 1000, mtry = 1))
pred <- predict(rf)
obs  <- df1$ACRES
plot(pred, type = 'l')
lines(obs, col = 'red')


############################################
############################################
#---------------#
#-- leftovers --#
#---------------#
## first differencing
df1 = data.frame(matrix(nrow = (nrow(historical)-1), ncol = 13))
for (i in 2:14){
   df1[,(i-1)] <- diff(historical[,i])
}
names(df1) = names(historical)[2:14]

## from first differences to original values
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

## ols regression with log response variable

lm1 <- lm(log(ACRES) ~ PPT_diffnorm_spring_summer + TAVG_diffnorm_spring_summer, df1); summary(lm1)
100 * (exp(lm1$coefficients[2]) -1)   #percent change in arithmetic mean per unit change in PPT
100 * (exp(lm1$coefficients[3]) -1)   #percent change in arithmetic mean per unit change in TAVG
toy1 <- data.frame(PPT_diffnorm_spring_summer  = c(-5,-2,-1,0,0,1,2,5),
          TAVG_diffnorm_spring_summer = c(5,2,1,0,0,-1,-2,-5))
predict(rlm1, toy1)