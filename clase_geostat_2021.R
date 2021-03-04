#' ---
#' title: "Variografia de la respiracion de suelos ..."
#' author: "Ovalle, Ramos, Guevara, 2021"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---



#'Este codigo sirve para la clase de geostat 2021
#'y el tema es variogramas.
#'El objetivo de este codigo es dejar evidencia reproducible de resultados de analisis de variogramas
#'con datos de respiracion de suelos a diversas escalas de trabajo, con diferentes configuraciones de 
#'datos disponibles y a lo largo de diversos escenarios naturales.  
#'Feb, 2021

#+ <p>&nbsp;</p>
#+ <p>&nbsp;</p>

#' Metodologia:
#' Bases de datos


#+ <p>&nbsp;</p>
#+ <p>&nbsp;</p>



#'pregunta a R donde estamos
#+ , attr.source='.numberLines'
getwd()
#'definir directorio de trabajo
#+, attr.source='.numberLines startFrom="2"'
setwd("/Users/marioguevara/Downloads/clase")
#'librarias necesarias

#'install.packages(c('raster', 'automap', 'rgdal'))
#+, attr.source='.numberLines startFrom="3"'
library(raster)
library(rgdal)
library(automap)
library(gstat)
#'pregunta a R como funciona getData
#'?getData
#'descarga limite de pais
#+, attr.source='.numberLines startFrom="7"'
lim <- getData('GADM', country='COL', level=2)
#'library for global base maps
library(maps)
#'show the map
#'map('world')

#'que tipo de objeto es este?

class(lim)
#'leer la base de datos de trabajo

srdb <- read.csv('https://raw.githubusercontent.com/bpbond/srdb/master/srdb-data.csv')
class(srdb)
str(srdb)
#ver nombres
 names(srdb)
 #selecciona variables
  soilCN <-  srdb[c('Longitude' , 'Latitude', 'Soil_CN')]
#'resumen de la base de datos
summary(soilCN)
#'selecciona variables especificas de otra base de datos
 #'plot(soilCN$Longitude, soilCN$Latitude)
 #'nombres de las variables en base de datos 
 plot(srdb$Longitude, srdb$Latitude)
#'sobrepongan puntos
#'points(soilCN$Longitude, soilCN$Latitude, col='blue')
#'quitamos valores no asignados
soilCN <- na.omit(soilCN)
 #' conocer la distribucion estadistica de los datos
hist(soilCN$Soil_CN)
#'vemos en logaritmo
hist(log1p(soilCN$Soil_CN))
#'agregamos una columna con los datos transformados
soilCN$Soil_CN_log <- log1p(soilCN$Soil_CN)
#'resumen estadistico
summary(soilCN)
#desviacion standard
apply(soilCN, 2, sd)
 #'define un objeto espacial
 coordinates(soilCN) <- ~ Longitude + Latitude
 #'define CRS
 proj4string(soilCN) <- CRS( '+proj=longlat +datum=WGS84')
  library(automap)
 #+ results=FALSE, warnings=FALSE
variogram <- autofitVariogram(Soil_CN ~ 1, soilCN)

plot(Soil_CN)

names(variogram)
 df <- data.frame(variogram$exp_var)
 
 #varianza minima
 #nugget sill radio
 #tipo de modelo
 #correlacion
  soilCN_aoi_1 <- crop(soilCN, drawExtent())
 soilCN_aoi_1 <- crop(soilCN, drawExtent())

 
dim(soilCN)
cor(df$dist, df$gamma)
variogram$sserr


tratamiento <- c('logTransformed')

#'visualize the variogram
plot(variogram)

x <- soilCN
x <- spTransform(x, CRS('+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs '))
au = autoKrige(Soil_CN ~ , x)

#get digital elevation model with both with 

#get environmental covariates
#https://essd.copernicus.org/preprints/essd-2020-264/
covariates <- readRDS(url('https://www.hydroshare.org/resource/9f981ae4e68b4f529cdd7a5c9013e27e/data/contents/prediction_factors_15km.rds'))

x_pred <- stack(au$krige_output)
predicted <-    projectRaster(x_pred, crs=CRS(projection(soilCN)))

#'PREPARACION DE DATOS PARA VARIOGRAMA

#'DEFINIR PROYECCION GEOGRAFICA
#'reproyectar si es necesario 
#'UNIDADES METRICAS 

#'?autofitVariogram

#'library(automap)
