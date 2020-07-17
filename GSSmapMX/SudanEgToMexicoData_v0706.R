library(sp); library(foreign); library(rgdal); library(car);library(carData); library(maptools)
library(spacetime); library(gstat); library(automap);library(randomForest);library(fitdistrplus);
library(e1071); library(caret); library(raster); library(soilassessment); library(soiltexture);
library(GSIF); library(aqp); library(plyr); library(Hmisc); library(corrplot); library(factoextra)
library(spup); library(purrr); library(lattice);library(ncf);library(npsurv); library(lsei);
library(nnet); library(class); library(mda); library(RColorBrewer); library(vcd); library(grid); 
library(neuralnet);library(readxl); library(psych);library(qrnn); library(dplyr)

###
{
  predictors=readGDAL("inputs/predictors_fao/dem.tif");predictors$band1=ifelse(is.na(predictors$band1),mean(!is.na(predictors$band1)),predictors$band1); summary(predictors$band1)
names(predictors)[1] <- 'dem'
  predictors$max_temp=readGDAL("inputs/predictors_fao/Max_temperature.tif")$band1;predictors$max_temp=ifelse(is.na(predictors$max_temp),mean(!is.na(predictors$max_temp)),predictors$max_temp); summary(predictors$max_temp)
  predictors$mean_temp=readGDAL("inputs/predictors_fao/Mean_temperature.tif")$band1;predictors$mean_temp=ifelse(is.na(predictors$mean_temp),mean(!is.na(predictors$mean_temp)),predictors$mean_temp); summary(predictors$mean_temp)
  predictors$lcover=readGDAL("inputs/predictors_fao/lcover.tif")$band1;predictors$lcover=ifelse(is.na(predictors$lcover),mean(!is.na(predictors$lcover)),predictors$lcover); summary(predictors$lcover)
  predictors$min_temp=readGDAL("inputs/predictors_fao/Min_temperature.tif")$band1;predictors$min_temp=ifelse(is.na(predictors$min_temp),mean(!is.na(predictors$min_temp)),predictors$min_temp); summary(predictors$min_temp)
  predictors$precipitation=readGDAL("inputs/predictors_fao/Precipitation.tif")$band1;predictors$precipitation=ifelse(is.na(predictors$precipitation),mean(!is.na(predictors$precipitation)),predictors$precipitation); summary(predictors$precipitation)
  predictors$swir1=readGDAL("inputs/predictors_fao/Band6.tif")$band1;predictors$swir1=ifelse(is.na(predictors$swir1),mean(!is.na(predictors$swir1)),predictors$swir1); summary(predictors$swir1)
  predictors$swir2=readGDAL("inputs/predictors_fao/Band7.tif")$band1;predictors$swir2=ifelse(is.na(predictors$swir2),mean(!is.na(predictors$swir2)),predictors$swir2); summary(predictors$swir2)
  predictors$BBlue=readGDAL("inputs/predictors_fao/Band3.tif")$band1;predictors$BBlue=ifelse(is.na(predictors$BBlue),mean(!is.na(predictors$BBlue)),predictors$BBlue); summary(predictors$BBlue)
  predictors$BGreen=readGDAL("inputs/predictors_fao/Band4.tif")$band1;predictors$BGreen=ifelse(is.na(predictors$BGreen),mean(!is.na(predictors$BGreen)),predictors$BGreen); summary(predictors$BGreen)
  predictors$BRed=readGDAL("inputs/predictors_fao/Band1.tif")$band1;predictors$BRed=ifelse(is.na(predictors$BRed),mean(!is.na(predictors$BRed)),predictors$BRed); summary(predictors$BRed)
  predictors$BIRed=readGDAL("inputs/predictors_fao/Band2.tif")$band1;predictors$BIRed=ifelse(is.na(predictors$BIRed),mean(!is.na(predictors$BIRed)),predictors$BIRed); summary(predictors$BIRed)
}



# predictors$dem=predictors$band1
# predictors$band1=NULL

###### PREDICTOR DATA PREPARATION #####
#Step 4: Assess the distribution of the predictors and transform where necessary 
#Check and remove NA if available in any layer

summary(predictors)
  {
predictors$precipitation=ifelse((predictors$precipitation)<10,10+((predictors$precipitation)),predictors$precipitation)
 predictors$BBlue=ifelse((predictors$BBlue)<0,-1*((predictors$BBlue)),predictors$BBlue)
 predictors$BGreen=ifelse((predictors$BGreen)<0,-1*((predictors$BGreen)),predictors$BGreen)
 predictors$BRed=ifelse((predictors$BRed)<0,-1*((predictors$BRed)),predictors$BRed)
 predictors$BIRed=ifelse((predictors$BIRed)<0,-1*((predictors$BIRed)),predictors$BIRed)
 predictors$swir1=ifelse((predictors$swir1)<0,-1*((predictors$swir1)),predictors$swir1)
 predictors$swir2=ifelse((predictors$swir2)<0,-1*((predictors$swir2)),predictors$swir2)
predictors$BIRed=predictors$BIRed+0.001; predictors$BRed=predictors$BRed+0.001; predictors$swir1=predictors$swir1+0.001
 }
# #Step 5: Derive the remote sensing indices of salinity
predictors$SI1=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI1");summary(predictors$SI1)
predictors$SI2=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI2");summary(predictors$SI2)
predictors$SI3=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI3");summary(predictors$SI3)
predictors$SI4=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI4");summary(predictors$SI4)
predictors$SI5=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI5");summary(predictors$SI5)
predictors$SI6=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI6");summary(predictors$SI6)
predictors$SAVI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SAVI");summary(predictors$SAVI)
predictors$VSSI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"VSSI");summary(predictors$VSSI)
predictors$NDSI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"NDSI");summary(predictors$NDSI)
predictors$NDVI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"NDVI");summary(predictors$NDVI)
predictors$SR=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SR");summary(predictors$SR)
predictors$CRSI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"CRSI");summary(predictors$CRSI)
predictors$BI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"BI");summary(predictors$BI)

#Remove NA if available
predictors$CRSI=ifelse(is.na(predictors$CRSI),mean(!is.na(predictors$CRSI)),predictors$CRSI)
predictors$SI4=ifelse(is.na(predictors$SI4),mean(!is.na(predictors$SI4)),predictors$SI4)

##Check the distribution of the image indices
{
hist(predictors@data[,"NDVI"],cex=1)
summary(predictors$NDVI)
predictors$NDVI=sqrt(predictors$NDVI)
predictors$SI4=ifelse((predictors$SI4)<0,-1*sqrt(abs(predictors$SI4)),sqrt(predictors$SI4))
hist(predictors$SI5)
predictors$NDVI=ifelse(is.na(predictors$NDVI),mean(!is.na(predictors$NDVI)),predictors$NDVI)
summary(predictors$NDVI)
}


##Step 6: Convert the predictors to dataframe and perform PCA,
predicters=predictors@data[,c("SI1","SI2","SI3","SI4","SI5","SI6","SAVI","VSSI","NDSI","NDVI","SR","CRSI", "BI")]

##Step 7: Check and plot correlation between predictors and perform PCA analysis 
soil.cor=cor(predicters)
corrplot(soil.cor,method="number",number.cex = 0.7)

pca<-prcomp(predicters[], scale=TRUE)
fviz_eig(pca) #Plot the predictor importance to know number of PCs to include in the list of predictors
summary(pca)
fviz_pca_var(pca, col.var = "contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

#Check for variable contribution
{
  var_coord_func <- function(loadings, comp.sdev){loadings*comp.sdev}
  loadings <- pca$rotation;sdev <- pca$sdev
  var.coord <- t(apply(loadings, 2, var_coord_func, sdev)) 
  var.cos2 <- var.coord^2;  comp.cos2 <- apply(var.cos2, 2, sum)
  contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
  var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
  head(var.contrib[, 1:11])
}
var.pca=as.data.frame(var.contrib[, 1:11])
write.csv(var.pca,file="pca_variable_contribution.csv")

#Step 8: Convert the PCs to grid maps and attach their corresponding values to the soil dataset
Pred.pcs<-predict(pca,predicters[])
predictors@data$PCA1=Pred.pcs[,1]
predictors@data$PCA2=Pred.pcs[,2]
predictors@data$PCA3=Pred.pcs[,3]
predictors@data$PCA4=Pred.pcs[,4]
predictors@data$PCA5=Pred.pcs[,5]
predictors@data$PCA6=Pred.pcs[,6]
predictors@data$PCA7=Pred.pcs[,7]
predictors@data$PCA8=Pred.pcs[,8]


####



#libraries
library(raster)
library(rgdal)
#read dataset in a shapefile 
shape <- readOGR(dsn='inputs/datos_inegi/', layer="edaf_puntos_sii")
proj4string(shape)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#MexiCE limit 
#CEuntry limit from the global administrative areas project\
#https://gadm.org/
lim <- readRDS('gadm36_MEX_2_sp.rds')
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

hist(dat$CE030)
summary(dat$CE030)
soil1=subset(dat,!is.na(dat$CE030))
coordinates(soil1) <- ~ X+Y
proj4string(soil1)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')

predictors_b <- readRDS('predictors_mx_soil_map.rds')

soil1 <- spTransform(soil1, CRS(projection(predictors_b)))


bubble(soil1,"CE030", main="Harmonized Electrical Conductivity (0-30 cm)")

soil1$dummy=(soil1$CE030)+0.001# add "+0.001" if minimum in summary(soil1$CE030) is zero
hist(soil1$dummy, main="Frequency distribution (before transformation)", xlab="Harmonized EC (dS/m)")
soil1$Tran=(soil1$dummy^(as.numeric(car::powerTransform(soil1$dummy, family ="bcPower")["lambda"]))-1)/(as.numeric(car::powerTransform(soil1$dummy, family ="bcPower")["lambda"]))
hist(soil1$Tran, main="Frequency distribution (after transformation)",xlab="Harmonized EC (dS/m)")
crs(predictors); crs(soil1)


soil1$dummy=(soil1$CE30100)+0.001# add "+0.001" if minimum in summary(soil1$CE030) is zero
hist(soil1$dummy, main="Frequency distribution (before transformation)", xlab="Harmonized EC (dS/m)")
soil1$Tran2=(soil1$dummy^(as.numeric(car::powerTransform(soil1$dummy, family ="bcPower")["lambda"]))-1)/(as.numeric(car::powerTransform(soil1$dummy, family ="bcPower")["lambda"]))
hist(soil1$Tran2, main="Frequency distribution (after transformation)",xlab="Harmonized EC (dS/m)")
crs(predictors_b); crs(soil1)


#####
dummyRaster <- function(rast){
  rast <- as.factor(rast)
  result <- list()
  for(i in 1:length(levels(rast)[[1]][[1]])){
    result[[i]] <- rast == levels(rast)[[1]][[1]][i]
    names(result[[i]]) <- paste0(names(rast), 
                                 levels(rast)[[1]][[1]][i])
  }
  return(stack(result))
}

soil_map <- shapefile('inputs/mapa_suelos/eda251mcw.shp')


ref <- projectRaster(raster(predictors), crs=projection(soil_map))
soil_map$des <- as.factor(soil_map$DESCRIPCIO)
soilmap <- rasterize(soil_map, ref, "des")
writeRaster(soilmap, file='soil_map_rasterized_DESCRIPCIO.tif')

library(raster)
soilmap <- raster('soil_map_rasterized_DESCRIPCIO.tif')
soilmap_dummy <- dummyRaster(soilmap)

saveRDS(predictors, file='predictors_mx.rds')

predictors <- readRDS('predictors_mx.rds')

predictors_b <- as(stack(stack(predictors), projectRaster(soilmap_dummy,raster(predictors))), 'SpatialGridDataFrame'); saveRDS(predictors, file='predictors_mx_soil_map.rds')

predictors.ov=over(soil1, predictors_b)
soil1@data <- cbind(soil1@data, data.frame(predictors.ov) )

soil1a=soil1@data[,c("Tran", names(predictors_b))]
soil1a <- na.omit(as.data.frame(soil1a))

 soil1a <- soil1a[, colSums(soil1a != 0) > 0]
#regmodelSuit(soil1a,Tran , dem , max_temp , mean_temp , lcover , min_temp , precipitation , swir1 , swir2 , BBlue , BGreen , BRed , BIRed , SI1 , SI2 , SI3 , SI4 , SI5 , SI6 , SAVI , VSSI , NDSI , NDVI , SR , CRSI , BI , PCA1 , PCA2 , PCA3 , PCA4 , PCA5 , PCA6 , PCA7 , PCA8)

#Tran + dem + max_temp + mean_temp + lcover + min_temp + precipitation + swir1 + swir2 + BBlue + BGreen + BRed + BIRed + SI1 + SI2 + SI3 + SI4 + SI5 + SI6 + SAVI + VSSI + NDSI + NDVI + SR + CRSI + BI + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8

control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=2)
# run the RFE algorithm
results <- rfe(soil1a[1], soil1a[-1], sizes=c(1:8), rfeControl=control)
# summarize the results
print(results)


library(SuperLearner)
fit_cv<- CV.SuperLearner(
				     soil1a[,1],
                                     soil1a[-1],	
                                     V=5,
                                     SL.library=list("SL.xgboost", "SL.ranger",
                                                     "SL.ksvm", "SL.kernelKnn" ,"SL.bayesglm"))


ranger.ec=train(Tran~(dem + max_temp + mean_temp + lcover + min_temp + precipitation + swir1 + swir2 + BBlue + BGreen + BRed + BIRed + SI1 + SI2 + SI3 + SI4 + SI5 + SI6 + SAVI + VSSI + NDSI + NDVI + SR + CRSI + BI + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8),  data = soil1a,  method = "ranger", trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  
  # Model performance metrics
data.frame(
  RMSE = RMSE,
  Rsquare = R_square
)
  
}





















