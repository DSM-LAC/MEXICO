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
library(maps)
library(automap)

#funcion para medir tiempo

###benchmark tool


#funcion para medir tiempo
CronometroON<- function(){
      tic<-proc.time()[3]
      assign(".tic", tic, envir=baseenv())
      invisible(tic)
      }
CronometroOFF<- function(){
      tic <- get(".tic", envir=baseenv())
      toc<-proc.time()[3]-tic
      hrs<-as.integer(toc/3600)
      minu<- as.integer(((toc/3600)-hrs)*60)
      seg<- ((((((toc/3600)-hrs)*60)))-minu)*60
      time<-paste(as.character(hrs),"hrs ",as.character(minu),"min ",as.character(round(seg,digit=2)),"seg",sep="")
      return(time)
      }
CronometroON()
print('una oracion larga puede tardar mas que una oracion corta')
CronometroOFF()
CronometroON()
print('si o no?')
CronometroOFF()



#TODOS LOS INSUMOS
#'leer la base de datos de trabajo
srdb <- read.csv('https://raw.githubusercontent.com/bpbond/srdb/master/srdb-data.csv')
country <- 'FRA'

#'pregunta a R como funciona getData
#'?getData
#'descarga limite de pais
#+, attr.source='.numberLines startFrom="7"'
lim <- getData('GADM', country=country, level=2)
#'library for global base maps

#'show the map
#'map('world')

#'que tipo de objeto es este?

class(lim)
class(srdb)
str(srdb)
#ver nombres
 names(srdb)
 #selecciona variables
  RS <-  srdb[c('Longitude' , 'Latitude', 'Rs_annual')]
#'resumen de la base de datos
summary(RS)
#'selecciona variables especificas de otra base de datos
 #'plot(soilCN$Longitude, soilCN$Latitude)
 #'nombres de las variables en base de datos 
 plot(srdb$Longitude, srdb$Latitude)
#'sobrepongan puntos
points(RS$Longitude, RS$Latitude, col='blue')
#'quitamos valores no asignados
RS <- na.omit(RS)
 #' conocer la distribucion estadistica de los datos
hist(RS$Rs_annual)
#'vemos en logaritmo
hist(log1p(RS$Rs_annual))
#'agregamos una columna con los datos transformados
RS$RS_log <- log1p(RS$Rs_annual)
#'resumen estadistico
summary(RS)
#desviacion standard
apply(RS, 2, sd)
 #'define un objeto espacial
 coordinates(RS) <- ~ Longitude + Latitude
 #'define CRS
 proj4string(RS) <- CRS( '+proj=longlat +datum=WGS84')
   #+ results=FALSE, warnings=FALSE
variogram <- autofitVariogram(Rs_annual ~ 1, RS)
plot(variogram)
#nugget/sill
#0.22 
#reproyecta el limite a un sistema de coordenadas plano
lim_proj <- spTransform(lim, CRS('+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs '))

#reproyecta los puntos limite a un sistema de coordenadas plano

RSproj <- spTransform(RS, CRS('+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs '))

#recorta el limite de Francia
RS_projFRA <- RSproj [lim_proj,]
#comienza a medir tiempo
CronometroON()
#kriging
au = autoKrige(Rs_annual ~ 1. , RS_projFRA, lim_proj )
#deja de medir tiempo
CronometroOFF()
###Ejercicio 1 como cambia la estructura espacial con el area de interes y con el numero de datos? 
###Ejercicio 2 como cambia el tiempo de computo con areas de distintos numeros de datos?
###Fin










correlaciones <- c(0.5, 0.4, 0.2)
ndatos <- c(100, 87, 15)

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

class()
dim()
str()
plot()

names


x_pred <- stack(au$krige_output)
predicted <-    projectRaster(x_pred, crs=CRS(projection(soilCN)))

#'PREPARACION DE DATOS PARA VARIOGRAMA

#'DEFINIR PROYECCION GEOGRAFICA
#'reproyectar si es necesario 
#'UNIDADES METRICAS 

#'?autofitVariogram

#'library(automap)
