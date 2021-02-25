#' ---
#' title: "Clase geostat 2021 PCT"
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
#'pregunta a R como funciona getData
#'?getData
#'descarga limite de pais

lim <- getData('GADM', country='MEX', level=2)
#'que tipo de objeto es este?

class(lim)
#'leer la base de datos de trabajo

srdb <- read.csv('https://raw.githubusercontent.com/bpbond/srdb/master/srdb-data.csv')

#'nombres de las variables en base de datos 

 names(srdb)
 #'selecciona variables especificas de otra base de datos
 
 soilCN <-  srdb[c('Longitude' , 'Latitude', 'Soil_CN')]
#'quitamos valores no asignados
soilCN <- na.omit(soilCN)
 #'define un objeto espacial
 coordinates(soilCN) <- ~ Longitude + Latitude
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

#' con esta linea hacen el pdf 
#' rmarkdown::render("copy_codigo_render.R")

