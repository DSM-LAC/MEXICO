#libraries
library(raster)
library(rgdal)
#read dataset in a shapefile 
shape <- readOGR(dsn = ".", layer = "edaf_puntos_sii")
proj4string(shape)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#Mexico limit 
#country limit from the global administrative areas project\
#https://gadm.org/
lim <- readRDS('gadm36_MEX_2_sp.rds')
lim <- spTransform(lim, CRS(projection(shape)))
#visualize the points

#Figure 1

x <- shape
#saturate for improving visualization of variability below 30% 
x$CO[x$CO>30] <- 30
#x <- x[is.na(x@data$CO)==FALSE,]
bubble(x, "CO",
        panel=function(...) {
          sp.polygons(lim, fill='white')
          sp:::panel.bubble(...)
        }) 
#and histograms
hist(shape$CO)
#log transformed
hist(log1p(shape$CO))

#generate a year only colum
year <- strsplit(as.character(shape$FECHA), '/')
y <- numeric()
for (i in 1:length(year)){
y[i] <- year[[i]][3]} 
shape$year <- as.numeric(y)
#select only those years after 1998 
shape@data <- shape@data[shape@data$year > 1998,]
shape@data <- na.omit(shape@data)

#bulk density funciton (Saini, 1996)
estimateBD <- function(SOC, method="Saini1996"){
OM <- SOC * 1.724
if(method=="Saini_1996"){BD <- 1.62 - 0.06 * OM}
return(BD)
}
shape@data$BLD <- estimateBD(shape@data$CO, method="Saini_1996")

#coarse fragments data
shape$CRF <- as.numeric(shape$PEDREG)
shape$CRF[shape$CRF==1] <-  0
shape$CRF[shape$CRF==2] <-  20
shape$CRF[shape$CRF==3] <-  40
shape$CRF[shape$CRF==4] <-  60
shape$CRF[shape$CRF==5] <-  80

#fit continuous functions (soil depth relationships)

#algorthms for quantitative pedology
library(aqp)
sp4=shape@data
sp4$IDPROF <- paste0("IDPROF_", sp4$COORD_Y, "_", sp4$COORD_X)
#generate a soil profile collection object
depths(sp4) <- IDPROF  ~ LIM_SUP + LIM_INF
site(sp4) <- ~ COORD_X + COORD_Y
coordinates(sp4) <- ~ COORD_X + COORD_Y
proj4string(sp4)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
library(GSIF)
try(SOC <- mpspline(sp4, 'CO', d = t(c(0,100))))
try(BLD <- mpspline(sp4, 'BLD', d = t(c(0,100))))
try(CRF <- mpspline(sp4, 'CRF', d = t(c(0,100))))
dat <- data.frame(id = sp4@site$IDPROF,
X = sp4@sp@coords[,1],
Y = sp4@sp@coords[,2],
SOC = SOC$var.std[,1],
BLD = BLD$var.std[,1],
CRF = CRF$var.std[,1])
head(dat)

#Figure 2

agg <- slab(sp4, fm= ~ SOC+BLD+CRF)
library(lattice)

xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
             xlab='median bounded by 25th and 75th percentiles',
             lower=agg$p.q25, upper=agg$p.q75, ylim=c(105,-2),
             panel=panel.depth_function,
             alpha=0.25, sync.colors=TRUE,
             par.settings=list(superpose.line=list(col=c('gray'), lwd=2)),
             prepanel=prepanel.depth_function,
             cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
             layout=c(1, 3), strip=strip.custom(bg=grey(0.8)),
             scales=list(x=list(tick.number=4, alternating=3, relation='free'))
             )

#Remove zero values

dat$SOC[dat$SOC==0] <- NA
dat <- na.omit(dat)

#Calculate SOC stocks
OCSKGM <- OCSKGM(ORCDRC = dat$SOC*10, BLD = dat$BLD*1000,
CRFVOL = dat$CRF, HSIZE = 100)
dat$OCSKGM <- OCSKGM
dat$meaERROR <- attr(OCSKGM,"measurementError")
dat <- dat[dat$OCSKGM>0,]
summary(dat)

#covariates files by state
(lis <- list.files(pattern='COV.rds'))
 (lis <- lis[-c(8,11, 12, 34)])
#empty data frame
m <- readRDS(paste0( lis[1]))
	m1 <- raster::as.data.frame(m[1,][-1,])
	m1$SOC <- numeric()
#prepare for extraction
d1 <- dat
coordinates(d1) <- ~X+Y
proj4string(d1)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
d1 <- spTransform(d1, CRS(projection(m)))
#loop for extraction
for(i in 1:length(lis)){
m <- readRDS(lis[i])
 proj4string(m) <- CRS( "+proj=longlat +ellps=GRS80 +no_defs")
ov <- over(d1, m)
	  ov$x <-raster::as.data.frame(d1)$X
	  ov$y <- raster::as.data.frame(d1)$Y
	  ov$SOC <- d1$OCSKGM 
	  ov <- ov[complete.cases(ov),]
	  m1 <- rbind(m1, ov)
	print(i)
}

##MODELING
Train <- m1[ trainIndex,]
Test  <- m1[-trainIndex,]

##relax line

library(caret)
library(raster)
library(doMC)
library(doParallel)
#run model in paralel
#Test <- readRDS("TEST_SOCsaini.rds" )
#Train <- readRDS("TRAIN_SOCsaini.rds")
set.seed(3456)
cl <- makeCluster(detectCores(), type='PSOCK')
	registerDoParallel(cl)
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5,saveDetails=TRUE )
	(RFE <- rfe(Train[,1:12], log1p(Train[,15]), sizes=c(1:6), rfeControl=control) )
stopCluster(cl = cl)

saveRDS(RFE, file='model_validation_modelingSainiRFE.rds')
saveRDS(Train, file='Train_modelingSainiRFE.rds')
saveRDS(Test, file='Test_modelingSainiRFE.rds')

##now with all data
set.seed(3456)
cl <- makeCluster(detectCores(), type='PSOCK')
	registerDoParallel(cl)
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5,saveDetails=TRUE )
	(RFE <- rfe(m1[,1:12], log1p(m1[,15]), sizes=c(1:6), rfeControl=control) )
stopCluster(cl = cl)

saveRDS(RFE, file='model_allData_modelingSainiRFE.rds')

##PREDICTION

library(foreach)
library(doSNOW)
library(caret)
library(rgdal)
library(raster)
#lis <- list.files(pattern='COV.rds')
#fit <- readRDS('modelRFE.rds')
fit <- readRDS('model_allData_modelingSainiRFE.rds')

(lis1 <- list.files(pattern='COV.rds'))

 (lis <- lis1[-c(8,11, 12, 15, 36,39, 40)])



NA2median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))

for(i in 30:length(lis)){
#i=36
print(lis[i])
testing <- readRDS(lis[i])
testing@data[] <- lapply(testing@data, NA2median)


#testing <- testing[1:(dim(testing)[1]/2),]
#testing <- testing[(dim(testing)[1]/2):(dim(testing)[1]),]
#gridded(testing) <- TRUE
#saveRDS(testing, file='Sonora1_2_COV.rds')

edo <- unlist(strsplit(lis[i], 'COV'))[1]
library(foreach)
library(doSNOW)
cl<-makeCluster(6) #change the 4 to your number of CPU cores
registerDoSNOW(cl)  
num_splits<-15
split_testing<-sort(rank(1:nrow(testing@data))%%4)
predictions<-foreach(i=unique(split_testing),
.combine=c,.packages=c("randomForest")) %dopar% {
as.numeric(predict(fit,newdata=testing@data[split_testing==i,]))

}
stopCluster(cl)

testing <- cbind(data.frame(testing@coords), SOC=expm1(predictions) , edo=edo)
write.csv(testing, file=paste0(edo, '_TODASsaini_soc90mRFE.csv'))
print(edo) 
}

#make rasters

library(rgdal)
library(raster)
require(data.table)
filenames <- list.files(path=getwd(), pattern='saini_soc90mRFE.csv', full.names=TRUE)
filenames
raslis <- list()
for (i in 1:length(filenames)){
d <- fread(filenames[i])
coordinates(d) <- ~x+y 
gridded(d) <- TRUE
r <- raster(d['SOC'])
writeRaster(r,
	file=paste0('tiles/',levels(as.factor(d$edo)), '_TODASsaini_soc90mRFE.tif'), 			
	overwrite=TRUE)
raslis[[i]] <- r
print(levels(as.factor(d$edo)))
}

require(data.table)
filenames <- list.files(path=getwd(), pattern='TODASsaini_soc90mRFE.csv', full.names=TRUE)
filenames
DAT <- rbindlist(lapply(filenames, fread))

library(dplyr)
library(lubridate)

DAT %>% 
  group_by(edo = as.factor(edo)) %>% 
  summarise(SOCsum = sum(SOC), SOCmin = min(SOC),
SOCmax = max(SOC), SOCsd = sd(SOC))     ->  DAT

x <-DAT
x %>% 
  group_by(edo = as.factor(as.character(edo))) %>% 
  summarise(SOCsum = sum(SOC), 
SOCmin = min(SOC),
SOCmax= max(SOC),
SOCsd = sd(SOC),
SOCmedian = median(SOC),
SOCmean = mean(SOC))     ->  xx

dat %>% group_by(edo) %>% tally() -> n

coordinates(testing) <- ~ x+y
gridded(testing) <- TRUE
#plot(raster(testing['prediction']))

writeRaster(raster(testing['prediction']), file=paste0(edo, '_TODAShoney_soc90mRFE.tif'), overwrite=TRUE)#change
print(edo)	
}

##Chihuahua, Coahuila, and Sonora working in two pieces

###Make mosaic of rasters 

ListRasters <- function(list_names) {
  raster_list <- list() # initialise the list of rasters
   for (i in 1:(length(list_names))){ 
    grd_name <- list_names[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name)
   }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

wgs84.tif.list <- list.files(path=getwd(), pattern=glob2rx("*.tif"), full.names=T,recursive=F)

list_names <- NULL
for (i in 1:length(wgs84.tif.list)) {
  list_names <- c(list_names, wgs84.tif.list[i])
}

raster.list <-sapply(list_names, FUN = ListRasters)

raster.list <- raslis
raster.list$fun <- mean

names(raster.list) <- NULL

raster.list$fun <- mean
mos <- do.call(mosaic, raster.list)
plot(mos)
#r <- mos
##


###calculate stocks
library(raster)
library(rgdal)
r <- raster("soc_predicted_SAINIbd.tif")

alt <- getData('alt', country='MEX')
slope <- terrain(alt, opt='slope')
 aspect <- terrain(alt, opt='aspect')
  hill <- hillShade(slope, aspect, 40, 270)

 
plot(hill, col=grey(0:100/100), legend=FALSE, main='', cex.axis=1.5)
plot(raster('soc_predicted_SAINIbd.tif'), col=rainbow(25, alpha=0.35), add=TRUE)


  l <- readRDS("gadm36_MEX_2_sp.rds" )
 projection(r) <- projection(l)

results <- data.frame()

lev <- levels(as.factor(l$NAME_1))

for (i in 1:length(lev)){
#i=1
r1 <- crop(r, l[l$NAME_1 == lev[i],])
r1 <- mask(r1, l[l$NAME_1 == lev[i],])

pixelsum <- cellStats(r1, sum)
 kgSOC <- (90*90) * pixelsum
 pgSOC <- kgSOC*1e-12
 #Petagrams of SOC at 1m depth across Mexico 
print(pgSOC)

results[i, 1] <- lev[i]
results[i, 2] <- pgSOC
results[i, 3] <- cellStats(r1, mean)
results[i, 4] <- cellStats(r1, min)
results[i, 5] <- cellStats(r1, max)
results[i, 6] <- cellStats(r1, sd)
results[i, 7] <- length(Which(r1, cells=TRUE))
print(lev[i])
}

names(results) <- c('estado', 'COS', 'mean', 'min', 'max', 'sd', 'npixeles')

#plot
set.seed(42) 
r <- read.csv('reporteCOSestado.csv')


r$qnt<- cut(r$COS , breaks=quantile(r$COS),
                                    labels=1:4, include.lowest=TRUE)

library(ggrepel)
 library(ggplot2)
 ggplot(r) + 
     geom_point(aes(log1p(npixeles), log1p(COS)), size = 5, color = 'grey') + 
     geom_label_repel( aes(log1p(npixeles), log1p(COS),  fill = factor(qnt), 
         label = estado), 
         fontface = 'bold', 
         color = 'white', 
         box.padding = unit(0.35, "lines"), 
         point.padding = unit(0.5, "lines") ) + 
     theme_classic(base_size = 20)



