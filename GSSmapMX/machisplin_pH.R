#libraries
library(spatialEco)
library(raster)
library(rgdal)
#read dataset in a shapefile 
#shape <- readOGR(dsn=getwd(), layer="edaf_puntos_sii")
shape <- readOGR(dsn='/home/mguevara/Downloads/mapa_salinidad_mx/inputs/datos_inegi/', layer="edaf_puntos_sii")
proj4string(shape)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#MexiCE limit 
#CEuntry limit from the global administrative areas project\
#https://gadm.org/
lim <- getData('GADM', country='MEX', level=2)
#lim <- readRDS('gadm36_MEX_2_sp.rds')
lim <- spTransform(lim, CRS(projection(shape)))
#visualize the points
#and histograms
hist(shape$CE)
#log transformed
hist(log1p(shape$CE))

#generate a year only CElum
year <- strsplit(as.character(shape$FECHA), '/')
y <- numeric()
for (i in 1:length(year)){
y[i] <- year[[i]][3]} 
shape$year <- as.numeric(y)
#select only those years after 1998 
shape@data <- shape@data[shape@data$year > 1998,]
shape@data <- na.omit(shape@data)

#algorthms for quantitative pedology
library(aqp)
sp4=shape@data
sp4$IDPROF <- paste0("IDPROF_", sp4$COORD_Y, "_", sp4$COORD_X)
#generate a soil profile CEllection object
depths(sp4) <- IDPROF  ~ LIM_SUP + LIM_INF
site(sp4) <- ~ COORD_X + COORD_Y
coordinates(sp4) <- ~ COORD_X + COORD_Y
proj4string(sp4)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
library(GSIF)
try(CE<- mpspline(sp4, 'CE', d = t(c(0,30,100,150))))
try(PH <- mpspline(sp4, 'PH', d = t(c(0,30,100,150))))
dat <- data.frame(id = sp4@site$IDPROF,
X = sp4@sp@coords[,1],
Y = sp4@sp@coords[,2],
CE030 = CE$var.std[,1],
PH030 = PH$var.std[,1],
CE30100 = CE$var.std[,2],
PH30100 = PH$var.std[,2])

hist(dat$PH030)
summary(dat$PH030)
soil1=subset(dat,!is.na(dat$PH030))
coordinates(soil1) <- ~ X+Y
proj4string(soil1)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')

predictors_b <- readRDS('predictors_mx_soil_map.rds')
soil1 <- spTransform(soil1, CRS(projection(predictors_b)))
bubble(soil1,"PH030", main="Harmonized Electrical Conductivity (0-30 cm)")

crs(predictors_b); crs(soil1)


predictors.ov=over(soil1, predictors_b)
soil1@data <- cbind(soil1@data, soil1@coords ,data.frame(predictors.ov) )

#soil1a=soil1@data[,c('X', 'Y', "PH030", names(predictors_b))]
soil1a=soil1@data[,c('X', 'Y', "PH30100", names(predictors_b))]

soil1a <- na.omit(as.data.frame(soil1a))
soil1a <- soil1a[, colSums(soil1a != 0) > 0]


lim_g <- spTransform(lim, CRS(projection(predictors_b)))
#Mydata <- data.frame(long=soil1a$X, lat=soil1a$Y, Tran=soil1a$PH030) 
Mydata <- data.frame(long=soil1a$X, lat=soil1a$Y, Tran=soil1a$PH30100) 
raster_stack <- stack(predictors_b)
#raster_stack <- aggregate(raster_stack, 5, mean)
raster_stack <- mask(raster_stack, lim_g)
remove(predictors_b)
library(MACHISPLIN)
interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=1, tps=TRUE)
#interp.rast <- readRDS('machinspline_CE_model_result.rds')
library(car)

(pred_sc <- interp.rast[[1]]$final)
pred_sc_2 <- raster.transformation(pred_sc, min(na.omit(soil1$PH030)),max(na.omit(soil1$PH030)),  trans = 'stretch')
#
library(wesanderson)
plot(pred_sc_2, col= viridis::viridis(10), zlim=c(0,1), legend=FALSE)

legend("topright", legend = c("pH"), fill = viridis::viridis(10)))

saveRDS(interp.rast, file='PH_30100_machisplin_model_object.rds')
writeRaster(stack(pred_sc, pred_sc_2), "PH_30100_machisplin_prediction_prediction_scaled.tif")

#saveRDS(interp.rast, file='PH_030_machisplin_model_object.rds')
#writeRaster(stack(pred_sc, pred_sc_2), "PH_030_machisplin_prediction_prediction_scaled.tif")
