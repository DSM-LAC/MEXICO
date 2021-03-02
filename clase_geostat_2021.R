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
setwd("/Users/marioguevara/Downloads/clase")
#'librarias necesarias

#'install.packages(c('raster', 'automap', 'rgdal'))
library(raster)
library(rgdal)
library(automap)
library(gstat)
#'pregunta a R como funciona getData
#'?getData
#'descarga limite de pais

lim <- getData('GADM', country='COL', level=2)
#'library for global base maps
library(maps)
#show the map
map('world')

#'que tipo de objeto es este?

class(lim)
#'leer la base de datos de trabajo

srdb <- read.csv('https://raw.githubusercontent.com/bpbond/srdb/master/srdb-data.csv')
class(srdb)
str(srdb)
str(srdb)
#'nombres de las variables en base de datos 
plot(srdb$Longitude, srdb$Latitude)
#'sobrepongan puntos
points(soilCN$Longitude, soilCN$Latitude, col='blue')
 names(srdb)
  soilCN <-  srdb[c('Longitude' , 'Latitude', 'Soil_CN')]
#'resumen de la base de datos
summary(soilCN)
#'selecciona variables especificas de otra base de datos
 #'plot(soilCN$Longitude, soilCN$Latitude)
 
 #' conocer la distribucion estadistica de los datos
hist(soilCN$Soil_CN)

hist(log1p(soilCN$Soil_CN))

soilCN$Soil_CN_log <- log1p(soilCN$Soil_CN)

#desviacion standard
apply(soilCN, 2, sd)
class(soilCN)
#'quitamos valores no asignados
soilCN <- na.omit(soilCN)
 #'define un objeto espacial
 coordinates(soilCN) <- ~ Longitude + Latitude
 #'define CRS
 proj4string(soilCN) <- ~ '+proj=longlat +datum=WGS84'
  library(automap)
 #+ results=FALSE, warnings=FALSE
variogram = autofitVariogram(Soil_CN ~1, soilCN)

plot(variogram)
#'PREPARACION DE DATOS PARA VARIOGRAMA

#'DEFINIR PROYECCION GEOGRAFICA
#'reproyectar si es necesario 
#UNIDADES METRICAS 

#'?autofitVariogram

#'library(automap)



