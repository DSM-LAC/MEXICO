

#load libraries
library(foreach)
library(doSNOW)
library(caret)
library(rgdal)
library(raster)
#covariate files list
(lis1 <- list.files(pattern='COV.rds'))
(lis <- lis1[-c(8,11, 12, 15, 36,39, 40)])

#r <- raster("soc_predicted_SAINIbd.tif")
#lim <- readRDS( "gadm36_MEX_2_sp.rds")

#import the models

model <- readRDS( file='model_validation_modelingSainiRFE.rds')

#import the training and testing dataset
train <- readRDS(file='Train_modelingSainiRFE.rds')
test <- readRDS(file='Test_modelingSainiRFE.rds')
#predict the models to the test dataset
library(randomForest)
test1<- cbind(test[,1:14], res = abs(log1p(test$SOC)-predict(model, test)))
#run a model for the residuals
set.seed(3456)
cl <- makeCluster(detectCores(), type='PSOCK')
	registerDoParallel(cl)
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5,saveDetails=TRUE )
	(fit <- rfe(test1[,1:12], log1p(test1[,15]), sizes=c(1:6), rfeControl=control) )
stopCluster(cl = cl)
#define median values for NAs in the covariates file
NA2median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
#loop a statewise prediction
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
write.csv(testing, file=paste0(edo, '_TODASsaini_soc90mRFE_RESIDUALS.csv'))
print(edo) 
}



library(rgdal)
library(raster)
require(data.table)
filenames <- list.files(path=getwd(), pattern='TODASsaini_soc90mRFE_RESIDUALS.csv', full.names=TRUE)
filenames
raslis <- list()
for (i in 1:length(filenames)){
d <- fread(filenames[i])
coordinates(d) <- ~x+y 
gridded(d) <- TRUE
r <- raster(d['SOC'])
writeRaster(r,
	file=paste0('tiles/',levels(as.factor(d$edo)), '_TODASsaini_soc90mRFE_RESIDUALS.tif'), 			
	overwrite=TRUE)
raslis[[i]] <- r
print(levels(as.factor(d$edo)))
}

raslis$fun <- mean
names(raslis) <- NULL
raslis$fun <- mean
mos <- do.call(mosaic, raslis)


