#Pilar Durante

####################RANGO DE 0-20cm###########
dat <- readRDS("D:/_Pitu/_Tesis2016/_Tesis_PD/Bases de datos/_4Murcia/_BD/dat_subsetPrueba.rds")

#Desecho los horizontes a mayores de 20cm de profundidad
   dat_20<-dat_subset[dat_subset$top<=20,]
   for (i in 1:nrow(dat_20)) {
     if (dat_20$top[i]<20 & dat_20$low[i]>=20) {
       dat_20$SOC_20[i]<-(((20-dat_20$top[i])/20)*dat_20$SOC[i])
     } else {
       dat_20$SOC_20[i]<-(((dat_20$low[i]-dat_20$top[i])/20)*dat_20$SOC[i])
     }
   }  
 #Profile colection y sumo todos
   dat_aqp <- dat_20
   depths(dat_aqp) <- id ~ top + low
   site(dat_aqp) <- ~ x + y + soilGroup
   coordinates(dat_aqp) <- ~ x + y
   
   dat_aqp$SOC_20t <- profileApply(dat_aqp, FUN=function(x) sum(x$SOC_20))  

   
   dat_ponderda <- data.frame(id = dat_aqp@site$id,
                            x = dat_aqp@sp@coords[,1],
                            y = dat_aqp@sp@coords[,2],
                            SOC_20 =dat_aqp@site$SOC_20t)
   write.csv(dat_ponderda, "dat_PruebaPOND.csv")
