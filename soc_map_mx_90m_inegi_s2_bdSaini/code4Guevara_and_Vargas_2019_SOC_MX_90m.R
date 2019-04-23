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

