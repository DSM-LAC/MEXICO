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



x_pred <- stack(au$krige_output)
predicted <-    projectRaster(x_pred, crs=CRS(projection(soilCN)))

#'PREPARACION DE DATOS PARA VARIOGRAMA

#'DEFINIR PROYECCION GEOGRAFICA
#'reproyectar si es necesario 
#'UNIDADES METRICAS 

#'?autofitVariogram

#'library(automap)


###
###Google earth
###



// Step 1. Create filters


//

var palettes = require('users/gena/packages:palettes');
var palette_slp = palettes.colorbrewer.RdYlGn[9];

var palette_ndvi = palettes.misc.tol_rainbow[7];
// Center the map on the area of interest.

// Point of interest
var aoi = {
longitude: -117.37,
latitude: 33.98,
};

var radius_of_interest_meters = 50000;
var aoi = ee.Geometry.Point([aoi.longitude, aoi.latitude])
.buffer(radius_of_interest_meters);

var bound = aoi
//var bound = aoi.geometry().bounds();

Map.centerObject(bound, 11);
Map.addLayer(bound, {}, 'bound');

//Construct start and end dates:
var start = ee.Date('2017-01-01');
var finish = ee.Date('2017-12-31');

// Step 2. Create predictor layers

// Predictor 1 Load WorldClim BIO Variables V1
var bioclim = ee.Image('WORLDCLIM/V1/BIO')
            .clip(bound);
            
// Seperate bands for analysis
var bio01 = bioclim.select('bio01');
var bio02 = bioclim.select('bio02');
var bio03 = bioclim.select('bio03');
var bio04 = bioclim.select('bio04');
var bio05 = bioclim.select('bio05');            
var bio06 = bioclim.select('bio06');
var bio07 = bioclim.select('bio07');
var bio08 = bioclim.select('bio08');
var bio09 = bioclim.select('bio09');
var bio10 = bioclim.select('bio10');
var bio11 = bioclim.select('bio11');
var bio12 = bioclim.select('bio12');
var bio13 = bioclim.select('bio13');
var bio14 = bioclim.select('bio14');
var bio15 = bioclim.select('bio15');
var bio16 = bioclim.select('bio16');
var bio17 = bioclim.select('bio17');
var bio18 = bioclim.select('bio18');
var bio19 = bioclim.select('bio19');

//Map.addLayer(bio01, {}, 'bio01');

// Predictor 2. Elevation and terrain modelling

// Load SRTM Digital Elevation Data 30m
var elevation = ee.Image('USGS/SRTMGL1_003')
            .clip(aoi);
            
// Load Global SRTM Topographic Diversity
var topo_d = ee.Image('CSP/ERGo/1_0/Global/SRTM_topoDiversity')
            .clip(aoi);
            
// Load Global SRTM CHILI (Continuous Heat-Insolation Load Index)
var chili = ee.Image('CSP/ERGo/1_0/Global/SRTM_CHILI')
            .clip(aoi);
            
// Apply slope algorithm to elevation.
var slope = ee.Terrain.slope(elevation);

// Get the aspect (in degrees).
var aspect = ee.Terrain.aspect(elevation);

// Calculate curvature
var curvature = ee.Terrain.slope(slope);

// Calculate hillshade
var hillshade = ee.Terrain.hillshade(elevation);

// Display all elevation related variables
//Map.addLayer(elevation, {}, 'Elevation');
//Map.addLayer(topo_d, {}, 'Topographic Diversity');
//Map.addLayer(chili, {}, 'Continuous Heat-Insolation Load Index');
Map.addLayer(slope, { palette: palette_slp}, 'Slope');
//Map.addLayer(aspect, {}, 'Aspect');
//Map.addLayer(curvature, {}, 'Curvature');

// Predictor 3. Load AG100: ASTER Global Emissivity Dataset 100-meter V003
var emissivity = ee.Image('NASA/ASTER_GED/AG100_003')
            .clip(bound);

// Seperate bands for analysis
var temperature = emissivity.select('temperature');
var temperature_sdev = emissivity.select('temperature_sdev');
var temp_c = temperature.multiply(0.01).subtract(273.15);
var temp_c_sdev = temperature_sdev.multiply(0.01).subtract(273.15);

//Map.addLayer(temp_c, {}, 'temperature');

// Predictor 4. Precipitation
var chirps= ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD')
            .filterBounds(bound)
            .filterDate(start, finish);

var chirps_mean = chirps.select('precipitation')
          .reduce('mean')          
          .rename('mean_precipitation')
          .float();  
var chirps_max = chirps.select('precipitation')
          .reduce('max')          
          .rename('max_precipitation')
          .float(); 
var chirps_min = chirps.select('precipitation')
          .reduce('min')          
          .rename('min_precipitation')
          .float();           

var precip_mean = chirps_mean.clip(aoi);
var precip_max = chirps_max.clip(aoi);
var precip_min = chirps_min.clip(aoi);

//Map.addLayer(precip_mean, {}, 'precip_mean');
Map.addLayer(precip_max, {min: -20, max: 60}, 'precip_max');
//Map.addLayer(precip_min, {}, 'precip_mean');

// Predictor 5 Landsat 8 Annual NDVI
var ndvi = ee.ImageCollection('LANDSAT/LC08/C01/T1_ANNUAL_NDVI')
            .filterBounds(bound)
            .filterDate(start, finish)
            .reduce('mean')
            .clip(bound);
 
 
var  ndviVis = {
  min: 0.0,
  max: 0.5,
  palette:
      ['3d3d3d', '080a02', '106e12', '37a930', '03ff17'],
};
            
Map.addLayer(ndvi, ndviVis, 'ndvi');      

// Predictor 6 Forest cover
var hansen = ee.Image('UMD/hansen/global_forest_change_2017_v1_5').clip(aoi);
var treeCanopyCover = hansen.select('treecover2000');
var treeCanopyCoverVis = {
  min: 0.0,
  max: 100.0,
  palette:
      ['3d3d3d', '080a02', '080a02', '080a02', '106e12', '37a930', '03ff17'],
};

Map.addLayer(treeCanopyCover, treeCanopyCoverVis, 'Tree Canopy Cover');

// Predictor 6 Global Surface Water

var gsw = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').clip(aoi);
var surfacewater = gsw.select('occurrence');
var surfacewaterVis = {
  min: 0.0,
  max: 100.0,
  palette: ['ffffff', 'ffbbbb', '0000ff'],
};

Map.addLayer(surfacewater, surfacewaterVis, 'Surface Water');

// Predictor 7 Mangrove

var mangrove = ee.ImageCollection('LANDSAT/MANGROVE_FORESTS')
            .filterBounds(bound)
            ;
var mangrovesVis = {
  min: 0,
  max: 1.0,
  palette: ['d40115'],
};

//Map.addLayer(mangrove, mangrovesVis, 'Mangroves');

// Predictor 8 MODIS Leaf Area Index
var modislai = ee.ImageCollection('MODIS/006/MCD15A3H')
            .filterBounds(bound)
            .filterDate(start, finish)
            .reduce('mean')
            .clip(bound);
var fpar = modislai.select('Fpar_mean');
var lai = modislai.select('Lai_mean');
var defaultVisualizationVis = {
  min: 0.0,
  max: 100.0,
  palette: ['e1e4b4', '999d60', '2ec409', '0a4b06'],
};

Map.addLayer(lai, defaultVisualizationVis, 'Leaf Area Index');
Map.addLayer(fpar, defaultVisualizationVis, 'Fraction of Photosynthetically Active Radiation');  

// Predictor 9 Forest Canopy Height

var gfch = ee.Image('NASA/JPL/global_forest_canopy_height_2005').clip(bound);
var forestCanopyHeight = gfch.select('1');
var forestCanopyHeightVis = {
  min: 0.0,
  max: 30.0,
  palette: [
    'ffffff', 'fcd163', '99b718', '66a000', '3e8601', '207401', '056201',
    '004c00', '011301'
  ],
};

//Map.addLayer(forestCanopyHeight, forestCanopyHeightVis, 'Forest Canopy Height');

// Predictor 10 Accessibility to cities
var acc = ee.Image('Oxford/MAP/accessibility_to_cities_2015_v1_0').clip(bound);
var accessibility = acc.select('accessibility');
var accessibilityVis = {
  min: 0.0,
  max: 41556.0,
  gamma: 4.0,
};

//Map.addLayer(accessibility, accessibilityVis, 'Accessibility');

// Predictor 11 Mean Net Primary Production


var dataset = ee.ImageCollection('MODIS/006/MOD17A3H')
             .filterBounds(bound)
             .reduce('mean')
             .clip(bound)
;
var npp = dataset.select('Npp_mean');
var nppVis = {
  min: 0.0,
  max: 19000.0,
  palette: ['bbe029', '0a9501', '074b03'],
};

Map.addLayer(npp, nppVis, 'NPP');

// Step 3. Export the images, specifying scale and region.
Export.image.toDrive({
  image: elevation,
  description: 'elevation',
  maxPixels: 1e11,
  scale: 100,
  region: bound,
  folder: 'bio_enviro'
});
