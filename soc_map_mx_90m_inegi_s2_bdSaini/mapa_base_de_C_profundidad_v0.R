###Codigo de R de Josselin y Mario para generar un mapa nacional de carbono organico en el suelo, y un mapa
#nacional de profundidad efectiva de raices en una resolucion espacial de 1x1km. 
##Autores: Josselin Arias, Armando Navarrete, Mario Guevara et al., 2021




#carga librerias
library(raster)
library(rgdal)

dat <- read.csv('/Users/marioguevara/Downloads/sal/v0_carbono_anavarrete/soc_12_2020.xls - soc.csv')
#read prediction factors
covar <- readRDS('PREDICTORES_ACTUALIZADOS.rds')
plot(covar)
names(covar)

