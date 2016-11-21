#----------#
#--AUTHOR--#
#----------#
# Jason Schatz 
# 10/16


#--------------------#
#--FILE DESCRIPTION--#
#--------------------#
# This script uses three methods for interpolation/spatial prediction:
# (1) geoR (universal kriging)
# (2) gstat (regression kriging...equivalent to geoR but a bit faster/more flexible)
# (3) Spatial prediction with INLA


#----------------------#
#-- IMPORT LIBRARIES --#
#----------------------#
library(lattice)
library(sp)
library(gstat)
library(geoR)
library(raster)
library(reshape2)
library(rgdal)
library(maptools)
library(scales)
library(maps)
library(mapdata)
library(SDMTools)
library(INLA)

#///////////////////////////////////////#
#///////////////////////////////////////#
#/// (1) universal kriging with geoR ///#
#///////////////////////////////////////#
#///////////////////////////////////////#

#---------------#
#-- prep data --#
#---------------#
data <- read.csv("H:/Jason/Climate/GHCND/interpolations/NM_high_temperature_thresholds_1981-2010.csv", header = T)
shp1 <- readOGR("H:/Jason/GIS/Political_boundaries/NM_state_boundary", "NM_state")

locations  <- data
SID        <- as.data.frame(locations[,1])
names(SID) <- "SID"
locations  <- locations[,2:3]
locations  <- SpatialPointsDataFrame(locations,SID)

grid <- read.csv("H:/Jason/Climate/GHCND/interpolations/interpolation_grid.csv", header = T)
coordinates(grid) = ~LON + LAT
gridded(grid) = T
grid <- as(grid, "SpatialPixelsDataFrame") # to full grid

geoRdata <- as.data.frame(data)
geodata  <- as.geodata(geoRdata, 2:3, 5)

#--set covariates on data--#
LAT  <- geoRdata$LAT
LON  <- geoRdata$LON
ELEV <- geoRdata$ELEV

#--set covariates on grid--#
grids <- as.data.frame(grid)
grid2 <- data.frame(grids$LON,grids$LAT)
names(grid2) <- c("LON", "LAT")
LAT2  <- grids$LAT
LON2  <- grids$LON
ELEV2 <- grids$ELEV

#----------------------------------------------#
#-- fit variogram and select variogram model --#
#----------------------------------------------#
vario <- variog(geodata, estimator.type = "modulus", max = 20000)
plot(vario)

#--set sill and range--#
sill  <- 1700
range <- 6

#--evaluate maximum likelihood fits of exponential and spherical models--#
ml.1  <- likfit(geodata, trend = ~LAT + LON + ELEV, ini.cov.pars = c(sill, range),
                cov.model = "exponential")
ml.10 <- likfit(geodata, trend = ~LAT + LON + ELEV, ini.cov.pars = c(sill, range),
                cov.model="spherical")

c(ml.1$AIC,  ml.1$BIC)
c(ml.10$AIC, ml.10$BIC)

#--perform univerals kriging--#
kc <- krige.control(trend.d = ~LAT  + LON  + ELEV,
                    trend.l = ~LAT2 + LON2 + ELEV2, obj.model = ml.10)
uk <- krige.conv(geodata, loc = grid2, krige = kc)

#--pull out results--#
a <-uk$predict
uk.df <- data.frame(a,grid2$LON, grid2$LAT)
names(uk.df) <- c("UK", "X", "Y")
coordinates(uk.df) <- ~X + Y
gridded(uk.df) = T
uk.spdf <- as(uk.df, "SpatialPixelsDataFrame")
uk.r <- raster(uk.spdf)
values(uk.r)[values(uk.r) < 0] <- 0	#set negative vals to 0

#--map results--#
plot(uk.r)
plot(shp1, add = T)

#--write raster to file--#
writeRaster(uk.r, "days_over_90F", "GTiff", overwrite = T)

#----------------------#
#-- cross validation --#
#----------------------#
xvalid1 <- xvalid(geodata, model = ml.10)
print(1 - var(xvalid1$error) / var(xvalid1$data))
cor(xvalid1$dat, xvalid1$predicted)	                   #correlation
sqrt(sum(xvalid1$error ^ 2) / length(xvalid1$error))   #rmse




#/////////////////////////////////////////#
#/////////////////////////////////////////#
#/// (2) regression kriging with gstat ///#
#/////////////////////////////////////////#
#/////////////////////////////////////////#

#---------------#
#-- prep data --#
#---------------#
data <- read.csv("H:/Jason/Climate/GHCND/interpolations/NM_high_temperature_thresholds_1981-2010.csv", header = T)
coordinates(data) = ~ LON + LAT

grid <-read.csv("H:/Jason/Climate/GHCND/interpolations/interpolation_grid.csv", header = T)
coordinates(grid) = ~ LON + LAT
gridded(grid) = T
grid <- as(grid, "SpatialPixelsDataFrame") # to full grid

df <- data.frame(matrix(NA, nrow = 2, ncol = 68))
rownames(df) <- c("RMSE", "cor")


#--------------------#
#-- define formula --#
#--------------------#
formula1 <- as.formula(d90F ~ LAT + LON + ELEV)


#--------------------#
#-- fit variograms --#
#--------------------#
null.vgm <- vgm(var(data$d90F), "Sph", 6, nugget = 80)
vgm_r    <- fit.variogram(variogram(formula1, data = data), model = null.vgm)
plot(variogram(formula1, data = data), vgm_r)


#-----------------#
#-- run kriging --#
#-----------------#
lm1 <- lm(formula1, data = data)
lm_rk <- krige(formula1, locations = data, newdata = grid)
ok_rk <- krige(residuals(lm1) ~ 1, locations = data, newdata = grid, model = vgm_r)
lm_rk$var1.rk <- lm_rk$var1.pred + ok_rk$var1.pred

spplot(lm_rk, "var1.rk")
r2 <- raster(lm_rk[1])
#writeRaster(r1,"days_over_90","GTiff",overwrite=T)


#----------------------#
#-- cross-validation --#
#----------------------#
rk.cv <- krige.cv(formula = formula1, locations = data, model = vgm_r, nfold = 10)
rmse  <- sqrt(sum(rk.cv$residual ^ 2) / length(rk.cv$residual))
cor   <- cor(rk.cv$observed, rk.cv$observed - rk.cv$residual)
print(result <- rbind(rmse, cor))

#--------------------------------#
#-- geostatistical simulations --#
#--------------------------------#
##takes a long time to run##
#rk.sim<-krige(formula1,locations=data,newdata=grid,model=vgm_r,nsim=2,debug.level=-1)
#spplot(rk.sim)




#////////////////////////////////////#
#////////////////////////////////////#
#/// spatial prediction with INLA ///#
#////////////////////////////////////#
#////////////////////////////////////#
#spde: stochastic partial differential equations
#inla: integrated nested laplace approximations 

#---------------#
#-- prep data --#
#---------------#
data <- read.csv("H:/Jason/Climate/GHCND/interpolations/NM_high_temperature_thresholds_1981-2010.csv",header = T)
data$LAT1 <- data$LAT
coordinates(data) = ~ LON + LAT

grid <- read.csv("H:/Jason/Climate/GHCND/interpolations/interpolation_grid.csv", header = T)
grid$LAT1 = grid$LAT
coordinates(grid) = ~ LON + LAT
gridded(grid) = T
grid <- as(grid, "SpatialPixelsDataFrame") # to full grid

shp1 <- readOGR("H:/Jason/GIS/Political_boundaries/NM_state_boundary", "NM_state")
borders1 <- data.frame(t(bbox(shp1)))


#--------------#
#-- set defs --#
#--------------#
mesh1 <- inla.mesh.2d(loc = data, loc.domain = borders1, max.edge = 0.5); plot(mesh1)
spde1 <- inla.spde2.matern(mesh = mesh1, alpha = 1.5)


#--------------------------#
#-- SPDE w/no covariates --#
#--------------------------#
formula <- y ~ -1 + Intercept + f(spatial.field, model = spde)
A.pred  <- inla.spde.make.A(mesh = mesh1)
stack.pred <- inla.stack(data = list(y = NA),
	                     A = list(A.pred),
                         effects = list(c(s.index, list(Intercept = 1))),
                         tag = 'pred')
join.stack.noelev <- inla.stack(stack.est.noelev,stack.pred)
output <- inla(formula,
               data = inla.stack.data(join.stack.noelev, spde = spde1),
               family = "gaussian",
               control.predictor = list(A = inla.stack.A(join.stack.noelev), compute = TRUE))

index.pred     <- inla.stack.index(join.stack.noelev,"pred")$data
post.mean.pred <- output$summary.linear.predictor[index.pred,"mean"]
post.sd.pred   <- output$summary.linear.predictor[index.pred,"sd"]
proj.grid      <- inla.mesh.projector(mesh1, 
	                                  xlim = range(grid@coords[,1]), 
	                                  ylim=range(grid@coords[,2]),
	                                  dims=c(grid@grid@cells.dim[1], grid@grid@cells.dim[2]))
post.mean.pred.grid <- inla.mesh.project(proj.grid, post.mean.pred)
post.sd.pre.grid    <- inla.mesh.project(proj.grid, post.sd.pred)

r1 <- raster(t(post.mean.pred.grid))
r1 <- flip(r1, 2)
values(r1)[values(r1) < 0] <- 0	#set negative vals to 0
extent(r1) <- bbox(grid)
plot(r1)
plot(shp1, add=T)


#-----------------------#
#-- SPDE w/covariates --#
#-----------------------#
formula <- y ~ -1 + Intercept + ELEV + LAT1 + f(spatial.field, model = spde)
A.pred  <- inla.spde.make.A(mesh = mesh1, loc = grid@coords)
A.est   <- inla.spde.make.A(mesh = mesh1, loc = data@coords)
s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde1$n.spde)

stack.est <- inla.stack(tag = 'est',
                        data = list(y = NA),
                        A = list(A.est, 1),
                        effects = list(c(s.index, list(Intercept = 1)),
                        list(ELEV = data@data$ELEV, LAT1 = data@data$LAT1)))
stack.pred <- inla.stack(tag='pred',
                         data=list(y=NA),
                         A=list(A.pred, 1),
                         effects=list(c(s.index,list(Intercept=1)),
                         list(ELEV=grid@data$ELEV, LAT1 = grid@data$LAT1)))
join.stack <- inla.stack(stack.est,stack.pred)

##takes a long time to run##
output <- inla(formula, 
	           num.threads = 3,
	           data = inla.stack.data(join.stack, spde = spde1),
	           family = "gaussian",
	           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE))
index.pred     <- inla.stack.index(join.stack, "pred")$data
post.mean.pred <- output$summary.linear.predictor[index.pred, "mean"]
post.sd.pred   <- output$summary.linear.predictor[index.pred, "sd"]
proj.grid      <- inla.mesh.projector(mesh1,
	                                  xlim = range(grid@coords[,1]), 
	                                  ylim=range(grid@coords[,2]), 
	                                  dims=c(grid@grid@cells.dim[1], grid@grid@cells.dim[2]))
post.mean.pred.grid <- inla.mesh.project(proj.grid, post.mean.pred)
post.sd.pre.grid    <- inla.mesh.project(proj.grid, post.sd.pred)
