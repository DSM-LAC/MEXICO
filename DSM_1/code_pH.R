##################
###PASO 1 PREPARACIÓN DEL AMBIENTE DE TRABAJO EN R
###1.1 INSTALAR PAQUETES NECESARIOS (install.packages) SOLO UNA VEZ (RECUERDA CITAR APROPIADAMENTE CADA PAQUETE) Y CARGA SUS LIBRERIAS ASOCIADAS (library)
###1.2 DEFINE EL DIRECTORIO DE TRABAJO (setwd)
###1.3 CARGA LAS FUNCIONES QUE LOS AUTORES DEL PRESENTE EJERCICIO GENERARON PARA FACILITAR EL PROCEDIMIENTO (separar1, CronometroON/OFF, )
##################
install.packages(c("raster","Hmisc","randomForest","sp","ipred","gtools","foreign","rgdal","plyr","e1071", "gstat"))
library(randomGLM)
require(foreign)
library(sp)
library(rgdal)
library(raster)
library(gtools)
library(plyr)
library(randomForest)
library(ipred)
library(e1071)
library(gstat)
setwd("J:/USUARIOS/CONABIO/EJERCICIOS_inferencia espacial/Usobiodiv/recorte")#### set working directory / define el área de trabajo, copiar y pegar la dirección donde esten guardados los datos, notese que se usa la diagonal invertida "/"
separar1 <- function(dataframe, seed=NULL,trainsize=0.80) {   ### Divide datos disponibles para entrenar y validar ### cambiar en trainsize = 0.80 = 80%
  if (!is.null(seed)) set.seed(seed)                          ### Divide available data for train and test
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, round(trainsize*nrow(dataframe)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
  }
CronometroON<- function(){            ### CronometroOFF calculates the time spent since the last call of CronometroON
      tic<-proc.time()[3]             ### CronometroOFF calcula el tiempo transcurrido desde la última ves que CronometroON fue ejecutada
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
ListFiles<-function(ruta=getwd(), ext=""){ ### lists the file names with the especified extention(s) contained in a folder
files<-list.files(path=ruta)               ### función para listar todos los archivos de una carpeta qeu tengan la extensión especificada
files
files<-files[substr(files,start=nchar(as.character(files))-3,stop=nchar(as.character(files))) %in% ext]
return(files)
}
files<-ListFiles(ext=c(".tif", ,".img")) ### image files to be used as covariates / lee todos los archivos con extensión .tif que serán utilizados como coavariables  # podran incluirse todas las extensiones que se manejan en el paquete "raster"

##################
###PASO 2 BASES DE DATOS
###2.1 CARGA BASE DE DATOS Y GENERA LA ESTADÍSTICA DESCRIPTIVA DE LA VARIABLE OBJETIVO
###2.2 CONSTRUYE UN ESPACIO DE COVARIABLES RELACIONADO A LA BASE DE DATOS DE LA VARIABLE OBJETIVO
###2.3 DIVIDE EN DOS LA BASE DE DATOS (80% PARA ENTRENAR Y 20% PARA VALIDAR EL MODELO)
##################
##### loading field data / carga los datos de campo
data<-read.dbf("data.dbf",as.is=F)
str(data)
data<- data[which(data$NHORIZON == "1"),] ## selects non-repeating data points / selecciona puntos sin repetición # SELECCIÓN DE LOS DATOS DEL HORIZONTE 1,
data<-data.frame(X_COORD=data$X_COORD, Y_COORD=data$Y_COORD, ph=data$PH) ## create a new data frame only with the data to be used/elimina los atributos que no se utilizan y crea un nuevo data frame # EXTRAE COORDENADAS DE LA BASE DE DATOS ORIGINAL Y ESCOGE LA VARIABLE OBJETIVO (ph=data$PH)
data<-na.omit(data) ## omits incomplete records / omite los registros incompletos
###### build the histogrm of the target variable and basic statistics/ graficamos el histograma, calculamos la estadística básica de los datos a modelar
nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))
par(mar=c(3.1, 3.1, 1.1, 2.1))
boxplot(data$ph, horizontal=TRUE,  outline=TRUE, frame=F,ylim=c(2,12), col = "darkgray", main="Distribución estadística de datos de pH")
hist(data$ph,xlab = "PH", main="Histograma del pH para 467 datos", col='lightgrey')
summary(data) ## basic stats / estadística básica
sd(data)
###
ras<-raster("BalHid.tif") ### uses a raster image as base grid / carga la imagen para crear la malla base ## AQUI VA EL NOMBRA CON EXTENSIÓN DE ALGUNA DE LAS CAPAS AUXILARES (PREDICTORES) COMO REFERENCIA PARA LA EXTRACCIÓN
ras
net<-as.data.frame(ras,xy=T) ### se convierte la imagen en puntos XY
str(net) ## verifica la estructura de los datos
summary(net) ## calcula la estadistica básica de los datos, y vemos los puntos que no contienen dato (NA's)
net<-na.omit(net) ## eliminamos los NA's
net$BalHid<-NULL ## eliminamos columna
### for para la extraccion a puntos
for (i in 1:length(files))
  {
    print( files[i])
    raster<-raster(x=files[i])
    net[files[i]]<-extract(x=raster, y=net[,1:2])

    data[files[i]]<-extract(x=raster, y=data[,1:2])
  }
##### verificamos que la extracción se haya realizado
str(net)
summary(net)
str(data)
summary(data)
#### chages na data to 0.001 to avoid discontinuity in the map /
####cambia los NA's por 0.001 para evitar la dicontinuidad en la predicción dem mapa
data[is.na(data)]<-0.001
net[is.na(net)]<-0.001
#### start the modeling precess / empezamos el preceso de modelado
#### creates the train (80%) and test (20%) data sets / se crean dos nuevos conjuntos de datos el de entrenamiento (80%) y la prueba (20%)
#### a partir del original
datosNummodel <- separar1(data, seed=2) # se ejecuta la función para separar los datos
entrenamiento <- datosNummodel$trainset  # se crea la base de entranamiento
prueba <- datosNummodel$testset  # se crea la base de prueba
######## Plot to see the spatial distribution of train and test dataset /
######## grafica de la distribución espacial de los putnos de entrenamiento y prueba
plot(entrenamiento$X_COORD , entrenamiento$Y_COORD, xlab="X", ylab="Y", pch=1, cex=0.5,
main="Distribución espacial de los puntos de campo" )
par(new=TRUE)
plot(prueba$X_COORD , prueba$Y_COORD, pch=16, cex=0.5, xaxt='n', ann=FALSE, yaxt='n' )
legend(2400000,1250000,c("Entrenamiento n= 374","Prueba n= 93"), pch=c(1,16))
## verificamos que se hayan creado correctamente
str(entrenamiento)
str(prueba)

##################
###PASO 3 RANDOM FOREST
###3.1 CONSTRUYE UN MODELO PREDICTIVO EMPLEANDO VALIDACIÓN CRUZADA CON 10 FOLDS PARA IDENTIFICAR LOS PARÁMETROS MÁS IMPORTANTES (EN RANDOM FOREST mtry)
###2.2 GENERA UN MODELO PREDICTIVO CON LOS MEJORES PARAMETROS (mtry y ntree), UNA PREDICCIÓN AL ESPACIO DE COVARIABLES Y UN MAPA DE LA VARIABLE OBJETIVO
###3.3 VALIDA EL MAPA CON LOS DATOS QUE NO ENTRARON AL MODELO Y COMPARA LOS RESULTADOS DE LA VALIDACIÓN CRUZADA
##################
CronometroON()
rfTuning <- tune.randomForest(ph~., data = entrenamiento,ntree=1000, mtry=seq(from=2,to=10,by=1),
                      tunecontrol= tune.control(sampling="cross",cross=10))  ## random forest tuning / comando para encontrar los mejores parámetros del modelo (tuning) CAMBIAR cross POR EL FOLD DESEADO
t1<-CronometroOFF()
t1 # tiempo que tardo la afinación del modelo
CronometroON()
rfvc <- errorest(ph ~ ., data = entrenamiento,model = randomForest, ntree=1000,
               gamma=as.numeric(rfTuning$best.parameters[1]),
                 , estimator="cv", est.para=control.errorest(k=10,random=TRUE,predictions=TRUE)) ## random forest CV / comando para ejecutar una validación cruzada del modelo,
                                                                                               ## notese que se utilizan los parametros objeidos con la función tune.randomForest SI CAMBIO cross ARRIBA PONER EL k el mismo valor
t2<-CronometroOFF()
t2 # medicion del tiempo
# root mse / cálculo del RMSE
rmse <- rfvc$error
# correlación
corr <- cor(entrenamiento$ph, rfvc$predictions)
# print them / mostramos en pantalla los datos de validación cruzada root mean aquared error (rmse) y correlación entre observados y modelados BASE CROSS VALIDATION CRUZADA DEJANDO 10 FUERA
rmse
corr
####################################### run models / se entrenan los modelos utilizando los mejores parámetros
bestRF<-randomForest(ph~.,ntree=as.numeric(rfTuning$best.parameters[2]),mtry=as.numeric(rfTuning$best.parameters[1]),data=entrenamiento)
################################## start predictions and map creation / en esta sección se hará la predicción y creamos y guardamos el mapa
### prediction for neagtive WB / predicción para el balance hídrico negativo
cov<-data.frame(X_COORD=net$x, Y_COORD=net$y,net[,3:19]) ## crea nuevo data.frame para la predicción, debido a que se deben tener
								   ## los mismos nombres de los atributos en las bases
str(cov) ## verificamos
pred <- predict(bestRF, newdata =cov, se.fit=T) ## realizamos la predicción
###  map creation / creamos el mapa
map<-data.frame(z=pred , cov[,1:2])
str(map)
coordinates(map)=~X_COORD+Y_COORD
gridded(map)<-TRUE
proj4string(map) <- CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs") # definimos proyección del mapa
r<-raster(map)
plot(r)
#writeRaster(r,"randomForestPH.tif") ## save map / guardamos el mapa ## GUARDA EL MAPA GENERADO EN FORMATO TIF EN LA CARPETA DE TRABAJO
########### predict values in the test dataset to validate / predecimos el modelo en la base de prueba para validar
predTest <- predict(bestRF, newdata =prueba, se.fit=T) ## realizamos la predicción
prueba$RF<-predTest
str(prueba)
corTest<-cor(prueba$ph,prueba$RF)
corTest
rmseTest<-sqrt(mean((prueba$ph-prueba$RF)^2))
rmseTest
##corTest es correlación entre observados y modelados y rmseTest es el root mean squared error VALIDACIÓN EXTERNA BASADA EN LOS DATOS QUE NO ENTRARON A SER MODELADOS

##################
###PASO 4 KRIGING ORDINARIO
###4.1 EXTRAE LOS RESIDUALES DEL MODELO ANTERIOR Y ANALIZA SU ESTRUCTURA ESPACIAL
###4.2 ENTRENA UN MODELO GEOESTADISTICO SIMPLE DE LOS RESIDUALES (KRIGING ORDINARIO)
###4.3 AÑADE LOS RESIDUALES "KRIGEADOS" AL MAPA GENERADO CON RANDOM FOREST Y VALIDALO CON LA BASE DE PRUEBA
##################
########################### para el ressidual kriging
#### priemro calculamos los resudiales o erreres de la predicción
dat<-entrenamiento[,1:2]
dat$ph<-entrenamiento$ph
dat$RF<- predict(bestRF, newdata =entrenamiento, se.fit=T) ## realizamos la predicción
dat$Residual<-entrenamiento$ph-dat$RF
dat$RF<-NULL
str(dat)
######## el kriging necesita objetos espaciales por eso en estas lineas los creamos
coordinates(dat) <- ~X_COORD+Y_COORD ## convertimos a objeto espacial
dat=remove.duplicates(dat, zero = 0.0, remove.second = TRUE) ## verificamos que no existan datos duplicados
spplot(dat, zcol = "Residual",col.regions = bpy.colors(5)) ### grafico de los residuales
bubble(dat, "Residual", main = "SOIL PH residuals", col=c(1,1)) ### grafico de burbuja de los residuales

######### calculamos y ajustamos el el semi-variograma
gDepth <- gstat(formula = Residual~1, data = dat)
vgDepth <- variogram(gDepth)
plot(vgDepth, plot.nu = FALSE)
vgmDepth <- vgm(nugget=0.08, psill =0.04, range =100000, model = "Exp")
plot(vgDepth, vgmDepth)
vgmDepth <- fit.variogram(vgDepth, vgmDepth, fit.method=7)
plot(vgDepth, vgmDepth)
vgmDepth

### generamos un objeto nuevo donde haremos la predicción
mask=cov
names(mask)
### ESPACIO DE COVARIABLES (SIN COVARIABLES) PARA INDICAR EL ESPACIO DE PREDICCIÓN convertimos los objetos faltantes a espaciales
coordinates(mask) <- ~X_COORD+Y_COORD
gridded(mask)=TRUE
proj4string(mask) <- CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs")
proj4string(dat) <- CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs")

#### run the model / corremos el modelo
CronometroON()
Depthkrig <- krige(Residual~1, locations=dat, newdata=mask, model=vgmDepth)
t3<-CronometroOFF()
t3
spplot(Depthkrig, zcol="var1.pred") #### graficamos el mapa de las predicciones de kriging
########### guardaremos la predicción de kriging
###  map creation / creamos el mapa
map<-data.frame(z=Depthkrig$var1.pred , cov[,1:2])
str(map)
coordinates(map)=~X_COORD+Y_COORD
gridded(map)<-TRUE
proj4string(map) <- CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs") # definimos proyección del mapa
r1<-raster(map)
plot(r1)
#writeRaster(x=r1 ,filename="ResidualKrigingPH.tif", overwrite=TRUE) ## save map / guardamos el mapa
########## sum the residual kriging plus the random forest prediction / sumamos los residuales mas la predicción de random forest
FinalMap<-r+r1
plot(FinalMap)
#writeRaster(x=FinalMap,filename="RandomForestResidualKrigingPH.tif", overwrite=TRUE) ## save map / guardamos el mapa

########## validate FinalMap / validamos el mapa FinalMap
prueba$RF_RK<-extract(FinalMap, prueba[,1:2])
prueba2<-na.omit(prueba)
corTestRFRK<-cor(prueba2$ph,prueba2$RF_RK)
corTestRFRK
rmseTestRFRK<-sqrt(mean((prueba2$ph-prueba2$RF_RK)^2))
rmseTestRFRK     #VALIDACIÓN EXTERNA DEL MAPA DE PH GENERADO CON RANDOM FOREST MAS LAS SUMA DE SUS RESIDUALES "KRIGEADOS"
#### end / fin
#save.image("ObjRF_RKfinal2.RData") #### save all objetcs
