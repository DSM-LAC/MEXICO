#pregunta a R donde estamos
getwd()
#definir directorio de trabajo
setwd("/Users/marioguevara/Downloads/clase")
#librarias necesarias
install.packages(c('raster', 'automap', 'rgdal'))
library(raster)
library(rgdal)
library(automap)
#pregunta a R como funciona getData
#?getData
#descarga limite de pais
lim <- getData('GADM', country='MEX', level=2)
#que tipo de objeto es este?
class(lim)
#leer la base de datos de trabajo
srdb <- read.csv('https://raw.githubusercontent.com/bpbond/srdb/master/srdb-data.csv')
#nombres de las variables en base de datos 
 names(srdb)
 #selecciona variables especificas de otra base de datos
soilCN <-  srdb[c('Longitude' , 'Latitude', 'Soil_CN')]
#quitamos valores no asignados
soilCN <- na.omit(soilCN)
 #define un objeto espacial
 coordinates(soilCN) <- ~ Longitude + Latitude
variogram = autofitVariogram(Soil_CN ~1, soilCN)

#### PREPARACION DE DATOS PARA VARIOGRAMA

DEFINIR PROYECCION GEOGRAFICA
reproyectar si es necesario 
UNIDADES METRICAS 

#?autofitVariogram

library(automap)







