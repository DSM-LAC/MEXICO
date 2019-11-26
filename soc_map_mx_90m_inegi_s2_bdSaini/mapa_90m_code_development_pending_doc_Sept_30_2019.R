
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
#lista de covariables
lis <- readRDS('covariateList.rds')
#visualize the points
#figure 1
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
#bulk density funciton
estimateBD <- function(SOC, method="Saini1996"){
OM <- SOC * 1.724
if(method=="Saini1996"){BD <- 1.62 - 0.06 * OM}
if(method=="Drew1973"){BD <- 1 / (0.6268 + 0.0361 * OM)}
if(method=="Jeffrey1979"){BD <- 1.482 - 0.6786 * (log(OM))}
if(method=="Grigal1989"){BD <- 0.669 + 0.941 * exp(1)^(-0.06 * OM)}
if(method=="Adams1973"){BD <- 100 / (OM /0.244 + (100 - OM)/2.65)}
if(method=="Honeyset_Ratkowsky1989"){BD <- 1/(0.564 + 0.0556 * OM)}
return(BD)
}
shape@data$BLD <- estimateBD(shape@data$CO, method="Saini1996")
#coarse fragments data
shape$CRF <- as.numeric(shape$PEDREG)
shape$CRF[shape$CRF==1] <-  0
shape$CRF[shape$CRF==2] <-  20
shape$CRF[shape$CRF==3] <-  40
shape$CRF[shape$CRF==4] <-  60
shape$CRF[shape$CRF==5] <-  80

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


agg <- slab(sp4, fm= ~ CO + BLD + CRF)
library(lattice)

xyplot(top ~ p.q50 | variable, data=agg, ylab='Profundidad (cm)',
             xlab='Valor medio de la variable dentro de los 25th y 75th percentiles',
             lower=agg$p.q25, upper=agg$p.q75, ylim=c(105,-2),
             panel=panel.depth_function,
             alpha=0.25, sync.colors=TRUE,
             par.settings=list(superpose.line=list(col=c('brown'), lwd=2)),
             prepanel=prepanel.depth_function,
             cf=agg$contributing_fraction, cf.col='brown', cf.interval=5, 
             layout=c(3, 1), strip=strip.custom(bg=grey(0.8)),
             scales=list(x=list(tick.number=4, cex=1.5,alternating=3, relation='free'), y=list(cex=1.5))
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

#lis is a covariates files list names by state

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

library(SuperLearner)
x <-m1[,1:12]
x$response <-  m1[,15] 

fit_unique<- CV.SuperLearner(
				     log1p(x[,13]),
                                     x[,1:12],	
                                     V=5,
                                     SL.library=list("SL.xgboost", "SL.ranger",
                                                     "SL.ksvm", "SL.kernelKnn" ,"SL.bayesglm"))

plot(fit_unique, Y=log1p(m1[,15])) + theme_bw(base_size=20)


control <- trainControl(method="repeatedcv", number=5, repeats=5)
set.seed(7)
modelRanger <- train(log1p(response)~., data=x, method="ranger", trControl=control)
# train the GBM model
set.seed(7)
modelGbm <- train(log1p(response)~., data=x, method="gbm", trControl=control, verbose=FALSE)
# train the SVM model
set.seed(7)
modelSvm <- train(log1p(response)~., data=x, method="svmRadial", trControl=control)
#Bayes GLM
set.seed(7)
modelBayes <- train(log1p(response)~., data=x, method="bayesglm", trControl=control)
#Bayes ELM
#set.seed(7)
#modelElm <- train(response~., data=x, method="elm", trControl=control)
#keras model
set.seed(7)
modelKknn<- train(log1p(response)~., data=x, method="kknn", trControl=control)
# collect resamples
results <- resamples(list(RAN=modelRanger, GBM=modelGbm, KNN=modelKknn,SVM=modelSvm, BYS=modelBayes))
# summarize the distributions
summary(results)
# boxplots of results
bwplot(results)
# dot plots of results
dotplot(results)

fit_cv<- CV.SuperLearner(
				     log1p(x[,126]),
                                     x[,1:125],	
                                     V=5,
                                     SL.library=list("SL.xgboost", "SL.ranger",
                                                     "SL.ksvm", "SL.kernelKnn" ,"SL.bayesglm"))

plot(fit_cv)

fit_all<- SuperLearner(
				     log1p(x[,126]),
                                     x[,1:125],	
                                    
                                     SL.library=list("SL.xgboost", "SL.ranger",
                                                     "SL.ksvm", "SL.kernelKnn" ,"SL.bayesglm"))

#predict ... 


library(caret)
trainIndex <- createDataPartition(m1$SOC, p = .70, 
                                  list = FALSE, 
                                  times = 1)

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





#saveRDS(RFE, file='modelRFEv1DREW.rds')
#lis <- list.files(pattern='COV.rds')
#RFE <- readRDS('modelRFEv1.rds')
#lisras <- list()
#library(caret)
#library(raster)
#beginCluster()
#for(i in 2:length(lis)){
#m <- stack(readRDS(lis[i]))
#library(randomForest)
#p <- clusterR(m, predict, args=list(model=RFE))
#edo <- unlist(strsplit(lis[i], 'COV'))[1]
#writeRaster(p, file=paste0(edo, '_soc90mRFE.tif'))
#lisras[[i]] <- p
#print(i)
#}
#endCluster()

library(foreach)
library(doSNOW)
library(caret)
library(rgdal)
library(raster)
#lis <- list.files(pattern='COV.rds')
#fit <- readRDS('modelRFE.rds')
fit <- readRDS('model_allData_modelingSainiRFE.rds')
#(lis1 <- list.files(pattern='COV.rds'))
#(lis <- lis1[-c(8,11, 12, 15, 36,39, 40)])
NA2median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
for(i in 1:length(lis)){
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



writeRaster(raster(testing['prediction']), file=paste0(edo, '_TODASsaini_soc90mRFE.tif'), overwrite=TRUE)#change
print(edo)	
}

##Chihuahua and Sonora working in two pieces

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

plot(r, col=colorRampPalette(c("azure","darkgray" ,"gold4", "orange",'firebrick1'))(255), cex.axis=1.5)
plot(sh, col='black', border=NA, add=TRUE)
plot(sh2, col='blue', border=NA, add=TRUE)
map('world', add=TRUE)
scalebar(below="kilometers", type='bar', lonlat=TRUE, divs=2, cex=1.3)
plot(res, col=colorRampPalette(c("gray100","gray" , "orange",'firebrick1'))(255), cex.axis=1.5)
plot(sh, col='black', border=NA, add=TRUE)
plot(sh2, col='blue', border=NA, add=TRUE)
map('world', add=TRUE)
scalebar(below="kilometers", type='bar', lonlat=TRUE, divs=2, cex=1.3)

plot(res2, col=colorRampPalette(c("gray100","gray" , "orange",'firebrick1'))(255), cex.axis=1.5, legend.only=TRUE)

plot(r2, col=colorRampPalette(c("azure","darkgray" ,"gold4", "orange",'firebrick1'))(255), cex.axis=1.5, legend.only=TRUE)

sh <- shapefile('agebur15gw.shp')
sh2 <- shapefile('C_A_Etapa_I_a_VI_2016actual.shp')

e <- extent(hill)
e <- as(e, 'SpatialPolygons')
proj4string(e) <- CRS(projection(sh))
sh2 <- spTransform(sh2, CRS(projection(sh)))
 
plot(hill, col=grey(0:100/100), legend=FALSE, main='', cex.axis=1.5)
plot(raster('soc_predicted_SAINIbd.tif'), col=rainbow(25, alpha=0.35), add=TRUE)
plot(sh[e,], col='black', border=NA, add=TRUE)
plot(sh2[e,], col='blue', border=NA, add=TRUE)
maps::map('world', add=TRUE)
scalebar(500, below="kilometers", type='bar', lonlat=TRUE, divs=5, xy=click())

err <- raster("RFE_residuals_bd_saini.tif")
err <- expm1(err)
err <- pred/err
xxx[xxx>100] <- 100
extent()
plot(hill, col=grey(0:100/100), legend=FALSE, main='', cex.axis=1.5)
plot(xxx,  col=rev(heat.colors(100, alpha=0.50)), add=TRUE)
plot(sh[e,], col='black', border=NA, add=TRUE)
plot(sh2[e,], col='blue', border=NA, add=TRUE)
maps::map('world', add=TRUE)
scalebar(500, below="kilometers", type='bar', lonlat=TRUE, divs=5, xy=click())



##crops visualize smaller areas
library(raster)
library(rgeos)
points <- SpatialPoints(cbind(-100, 25))
pbuf <- gBuffer(points, widt=1)
#plot(pbuf, add=TRUE, col='green')
#alt <- getData('alt', country='MEX')
#hill <- hillShade(slope, aspect, 40, 270)

cr <- mask(crop(raster('soc_predicted_SAINIbd.tif'), pbuf), pbuf)
cr <- trim(cr)
cr2 <- mask(crop(raster("RFE_residuals_bd_saini.tif"), pbuf), pbuf)
cr2 <- trim(cr2)
cr3 <- mask(crop(hill, pbuf), pbuf)
cr3 <- trim(cr3)
proj4string(pbuf) <- CRS(projection(sh))

plot(cr3, col=grey(0:100/100), legend=FALSE, main='', cex.axis=1.5)
plot(cr2,  col=rev(heat.colors(100, alpha=0.50)), add=TRUE)#col=rainbow(25, alpha=0.35)
plot(sh[pbuf,], col='black', border=NA, add=TRUE)
plot(sh2[pbuf,], col='blue', border=NA, add=TRUE)
plot(pbuf, add=TRUE, cex=5)
#maps::map('world', add=TRUE)
#scalebar(50, below="kilometers", type='bar', lonlat=TRUE, divs=5, xy=click())


plot(cr3, col=grey(0:100/100), legend=FALSE, main='', cex.axis=1.5)
plot(cr,  col=rainbow(25, alpha=0.35), add=TRUE)
plot(sh[pbuf,], col='black', border=NA, add=TRUE)
plot(sh2[pbuf,], col='blue', border=NA, add=TRUE)
#maps::map('world', add=TRUE)
plot(pbuf, add=TRUE, cex=5)
#scalebar(50, below="kilometers", type='bar', lonlat=TRUE, divs=5, xy=click())



