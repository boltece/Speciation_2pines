#### install dependencies ####

install.packages("corrplot")
library(corrplot)
install.packages("ENMeval")
library(ENMeval)
install.packages("raster")
library(raster)
install.packages("sp")
library(sp)
install.packages("rgdal")
library(rgdal)
install.packages("dismo")
library(dismo)
install.packages("maptools")
library(maptools)
install.packages("rgeos")
library(rgeos)
install.packages("phyloclim")
library(phyloclim)
install.packages("rJava")
library(rJava)
install.packages("ENMTools")
library(ENMTools)
library(tibble)



#### occurrence data from GBIF ####

## sp thin to thin out occurrence data 1 point per 10km

install.packages("spThin")
library(spThin)

setwd("~/path/to/dir/location_data")


## Pinus rigida locations

Rig_loc <- read.csv("GBIF_prunded_Prigida.csv")

thin_R <- thin(Rig_loc, lat.col = "Latitude", long.col = "Longitude", spec.col = "Species",
               thin.par=10, reps=1, locs.thinned.list.return = FALSE, write.files = TRUE,
               max.files = 5, out.dir = "~/Desktop/JIB_attempt2", out.base = "thinned_Rdata",
               write.log.file = TRUE, log.file = "spatial_thin_log.txt",
               verbose = TRUE)

R_xy <- read.csv("thinned_Rdata_thin1.csv")
R_loc <- R_xy[,2:3] 
View(R_loc)
# delete (-88.16528, 38.77611) (-84.00000, 32.00000) (-76.00981,35.13000)
R_loc <- R_xy[,2:3]  
index <- which(R_loc$Longitude <(-88))
index 
R_loc <- R_loc[-index,]
View(R_loc)
index2 <-which(R_loc$Latitude == (32.00000))
index2
R_loc <- R_loc[-index2,]
index3 <-which(R_loc$Latitude == (35.13000))
index3
R_loc <- R_loc[-index3,]

Longitude <- as.numeric(R_loc[,1])
Longitude
Latitude<- as.numeric(R_loc[,2])
yz <- cbind(Longitude, Latitude)
yz




## Pinus pungens locations

Pung_loc <- read.csv("GBIF_Jetton_Ppungens.csv")

thin_P <- thin(Pung_loc, lat.col = "Latitude", long.col = "Longitude", spec.col = "Species",
               thin.par=10, reps=1, locs.thinned.list.return = FALSE, write.files = TRUE,
               max.files = 5, out.dir = "~/Desktop/JIB_attempt2", out.base = "thinned_Pdata",
               write.log.file = TRUE, log.file = "spatial_thin_log.txt",
               verbose = TRUE)


P_xy <- read.csv("thinned_Pdata_thin1.csv")
P_loc <- P_xy[,2:3] 
Longitude <- as.numeric(P_loc[,1])
Longitude
Latitude<- as.numeric(P_loc[,2])
xy <- cbind(Longitude, Latitude)
xy





#### correlation test of environmental variable ####

setwd("~/path/to/dir/WorldClim_data/wc2.1_2.5m_bio/")
bio2_present <- raster("wc2.1_2.5m_bio_2.tif")
bio10_present <- raster("wc2.1_2.5m_bio_10.tif")
bio11_present <- raster("wc2.1_2.5m_bio_11.tif")
bio15_present <- raster("wc2.1_2.5m_bio_15.tif")
bio17_present <- raster("wc2.1_2.5m_bio_17.tif")

stacked <- stack(bio2_present, bio10_present, bio11_present, bio15_present, bio17_present)
ext <- extent(-95, -60, 25, 55)
stack_cropped <- crop(stacked, ext)
bioclim <- unstack(stack_cropped)
writeRaster(bioclim[[1]],"~/path/to/dir/WorldClim_data/cropped_current_climate_XlgExtent/Bio2.asc" )
writeRaster(bioclim[[2]],"~/path/to/dir/WorldClim_data/cropped_current_climate_XlgExtent/Bio10.asc" )
writeRaster(bioclim[[3]],"~/path/to/dir/WorldClim_data/cropped_current_climate_XlgExtent/Bio11.asc" )
writeRaster(bioclim[[4]],"~/path/to/dir/WorldClim_data/cropped_current_climate_XlgExtent/Bio15.asc" )
writeRaster(bioclim[[5]],"~/path/to/dir/WorldClim_data/cropped_current_climate_XlgExtent/Bio17.asc" )

setwd("~/path/to/dir/WorldClim_data/cropped_current_climate_XlgExtent/")

bio2 <- raster("Bio2.asc")
bio10 <- raster("Bio10.asc")
bio11 <- raster("Bio11.asc")
bio15 <- raster("Bio15.asc")
bio17 <- raster("Bio17.asc")
env_stack_now <- stack( bio2, bio10, bio11, bio15, bio17)
crs(env_stacked_now) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
env_stack_now

occs.z <- cbind(xy, raster::extract(env_stack_now, xy))
write.csv(occs.z, "~/path/to/dir/Bioclim_data_Pungens_occs.csv")

occs.zR <- cbind(yz, raster::extract(env_stack_now, yz))
write.csv(occs.z, "~/path/to/dir/Bioclim_data_Rigida_occs.csv")

env_P <- occs.z[,3:7]
corr_P <- cor(env_P, use = "complete.obs", method="pearson")
corr_P_NA <- round(corr_P,2)
corr_P_NA[which(abs(corr_P_NA) < .75)] <- NA
abs(corr_P_NA)


env_R <- occs.zR[,3:7]
corr_R <- cor(env_R, use = "complete.obs", method="pearson")
corr_R_NA <- round(corr_R,2)
corr_R_NA[which(abs(corr_R_NA) < .75)] <- NA
abs(corr_R_NA)

#### Pinus pungens SDMs ####

setwd("~/path/to/dir/location_data")
P_xy <- read.csv("thinned_Pdata_thin1.csv")
P_loc <- P_xy[,2:3] 
Longitude <- as.numeric(P_loc[,1])
Longitude
Latitude<- as.numeric(P_loc[,2])
xy <- cbind(Longitude, Latitude)
xy




bb<- bbox(xy)
summary(bb)
bb.buf <- extent(bb[1]-5,bb[3]+5, bb[2]-5,bb[4]+5)
summary(bb.buf)
envs.backg <- crop(env_stack_now,bb.buf)  
plot(envs.backg[[1]], main=names(envs.backg)[1])
points(xy)
bg <- randomPoints(envs.backg[[1]], n=10000)
bg <- as.data.frame(bg)
head(bg)
plot(envs.backg[[1]], legend=FALSE)
points(bg, col='red')


block <- get.block(xy, bg, orientation = "lat_lon")
table(block$occs.grp)

evalplot.grps(pts = xy, pts.grp = block$occs.grp, envs = envs.backg) + 
  ggplot2::ggtitle("Spatial block partitions: occurrences")

evalplot.grps(pts = bg, pts.grp = block$bg.grp, envs = envs.backg) + 
  ggplot2::ggtitle("Spatial block partitions: background")




occs.z <- cbind(xy, raster::extract(env_stack_now, xy))
bg.z <- cbind(bg, raster::extract(env_stack_now, bg))

head(bg) 
colnames(bg) <- c("Longitude", "Latitude")

os <- list(abs.auc.diff = TRUE, pred.type = "cloglog", validation.bg = "partition")
os2 <- list(abs.auc.diff = TRUE, pred.type = "raw", validation.bg = "partition")
tune.args <- list(fc = c("L","LQ","LQH","H"), rm = 1:5)
ps <- list(orientation = "lat_lon")





#### NEW ENMeval SCRIPT 2.0 ####
### this is what was used for all model outputs

P_train <- ENMevaluate(occ=xy, env=env_stack_now, bg=bg, tune.args=tune.args, partitions = 'block', other.settings= os, partition.settings = ps, doClamp=TRUE, parallel = FALSE, algorithm = 'maxent.jar',overlap = TRUE)
# * Running ENMeval v2.0.0 with legacy arguments. These will be phased out in the next version.
#*** Running initial checks... ***
  
#  * Clamping predictor variable rasters...
#* Model evaluations with spatial block (4-fold) cross validation and lat_lon orientation...

#*** Running ENMeval v2.0.0 with maxent.jar v3.4.1 from dismo package v1.3.3 ***
#Calculating niche overlap for statistic D...
#|=============================================================================================================| 100%
#Calculating niche overlap for statistic I...
#|=============================================================================================================| 100%
#ENMevaluate completed in 11 minutes 32.3 seconds.


setwd("~/path/to/dir/")
write.csv(bg, "Ppungens_ENM_background_pts.csv")


pred <- P_train@predictions[[which(P_train@results$delta.AICc==0)]] #LQ_1
plot(pred)
pred

res <- eval.results(P_train)
res
opt.aicc <- res$delta.AICc==0
opt.aicc

m4.mx <- eval.models(P_train)[["fc.LQ_rm.1"]]
m4.mx # save html in Updated_script_models folder as Maxent_model_pungens_Delta0
plot(m4.mx)
eval.predictions(P_train)


## response curves

dismo::response(eval.models(P_train)[["fc.LQ_rm.1"]])
#saved as Rplot_P_train_LQ_1_pungens_response_curves
dismo::response(m4.mx, range="pa")
dismo::response(m4.mx, range="p")


write.csv(P_train@results,file="Ppungens_ENM_eval_results.csv")
eval_occs <- eval.occs(P_train)
write.csv(eval_occs, "eval_occs_P_train.csv")
eval.bg(P_train) %>% head()
eval_bg <- eval.bg(P_train)
write.csv(eval_bg, "eval_bg_P_train.csv")


#### null model ####
mod.null <- ENMnulls(P_train, mod.settings = list(fc = "LQ", rm = 1), no.iter = 100)
null.emp.results(mod.null)
null_emp_stats <- null.emp.results(mod.null)
write.csv(null_emp_stats, "null_versus_empirical_models_P_train.csv")
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
# saved as Rplot_null_versus_empirical_models_pungens_P_train.pdf

#### metadata for Pungens model ####
rmm <- eval.rmm(P_train)
rangeModelMetadata::rmmToCSV(rmm, "rmm_metadata_Pungens_P_train.csv")




#### Predictions ####

## Current climate distribution prediction
P_raw_now4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_now, args="outputformat=raw")
plot(P_raw_now4)
writeRaster(P_raw_now4, "~/Documents/TGG_SDM/Corrected_Model/Updated_script_models/pungens_raw_current.asc")
P_cloglog_now4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_now, args="outputformat=cloglog")
plot(P_cloglog_now4)
writeRaster(P_cloglog_now4, "~/Documents/TGG_SDM/Corrected_Model/Updated_script_models/pungens_cloglog_current.asc")



## mid-Holocene 
## CCSM
setwd("~/Documents/TGG_SDM/WorldClim_data/ccmidbi_2-5m/")
bio2c_mid <- raster("ccmidbi2.tif")
bio10c_mid <- raster("ccmidbi10.tif")
bio11c_mid <- raster("ccmidbi11.tif")
bio15c_mid <- raster("ccmidbi15.tif")
bio17c_mid <- raster("ccmidbi17.tif")

stacked_mid <- stack(bio2c_mid, bio10c_mid, bio11c_mid, bio15c_mid, bio17c_mid)
crs(stacked_mid) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_all <- crop(stacked, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1]) 
unstacked <- unstack(cropped_all)

setwd("~/Documents/TGG_SDM/WorldClim_data/cropped_midCCSM_climate_XlgExtent/")
writeRaster(unstacked[[1]], "~/path/to/dir/WorldClim_data/cropped_midCCSM_climate_XlgExtent/Bio2.asc")
writeRaster(unstacked[[2]], "~/path/to/dir/WorldClim_data/cropped_midCCSM_climate_XlgExtent/Bio10.asc")
writeRaster(unstacked[[3]], "~/path/to/dir/WorldClim_data/cropped_midCCSM_climate_XlgExtent/Bio11.asc")
writeRaster(unstacked[[4]], "~/path/to/dir/WorldClim_data/cropped_midCCSM_climate_XlgExtent/Bio15.asc")
writeRaster(unstacked[[5]], "~/path/to/dir/WorldClim_data/cropped_midCCSM_climate_XlgExtent/Bio17.asc")


setwd("~/Documents/TGG_SDM/WorldClim_data/cropped_midMIROC_climate_XlgExtent/")


bio2c <- raster("Bio2.asc")
bio10c <- raster("Bio10.asc")
bio11c <- raster("Bio11.asc")
bio15c <- raster("Bio15.asc")
bio17c <- raster("Bio17.asc")
env_stack_midCCSM <- stack( bio2c, bio10c, bio11c, bio15c, bio17c)
crs(env_stack_midCCSM) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


P_raw_mid4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_midCCSM, args="outputformat=raw")
plot(P_raw_mid4)
writeRaster(P_raw_mid4, "~/path/to/dir/pungens_raw_midHolocene_ccsm.asc")
P_cloglog_mid4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_midCCSM, args="outputformat=cloglog")
plot(P_cloglog_mid4)
writeRaster(P_cloglog_mid4, "~/path/to/dir/pungens_cloglog_midHolocene_ccsm.asc")

## MIROC-ESM
setwd("~/path/to/dir/WorldClim_data/mrmidbi_2-5m")

bio2mr_mid <- raster("mrmidbi2.tif")
bio10mr_mid <- raster("mrmidbi10.tif")
bio11mr_mid <- raster("mrmidbi11.tif")
bio15mr_id <- raster("mrmidbi15.tif")
bio17mr_mid <- raster("mrmidbi17.tif")

stacked_mid_mr <- stack(bio2mr_mid, bio10mr_mid, bio11mr_mid, bio15mr_mid, bio17mr_mid)
crs(stacked_mid_mr) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_all_mr <- crop(stacked_mid_mr, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1]) 
unstacked <- cropped_all_mr

setwd("~/path/to/dir/WorldClim_data/cropped_midMIROC_climate_XlgExtent/")
writeRaster(unstacked[[1]], "~/path/to/dir/WorldClim_data/cropped_midMIROC_climate_XlgExtent/Bio2.asc")
writeRaster(unstacked[[2]], "~/path/to/dir/WorldClim_data/cropped_midMIROC_climate_XlgExtent/Bio10.asc")
writeRaster(unstacked[[3]], "~/path/to/dir/WorldClim_data/cropped_midMIROC_climate_XlgExtent/Bio11.asc")
writeRaster(unstacked[[4]], "~/path/to/dir/WorldClim_data/cropped_midMIROC_climate_XlgExtent/Bio15.asc")
writeRaster(unstacked[[5]], "~/path/to/dir/WorldClim_data/cropped_midMIROC_climate_XlgExtent/Bio17.asc")

bio2m <- raster("Bio2.asc")
bio10m <- raster("Bio10.asc")
bio11m <- raster("Bio11.asc")
bio15m <- raster("Bio15.asc")
bio17m <- raster("Bio17.asc")


env_stack_midMIROC <- stack(bio2m, bio10m, bio11m, bio15m, bio17m)
crs(env_stack_midMIROC) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"



P_raw_mid4_miroc <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_midMIROC, args="outputformat=raw")
plot(P_raw_mid4_miroc)
writeRaster(P_raw_mid4_miroc, "~/path/to/dir/pungens_raw_midHolocene_miroc.asc")
P_cloglog_mid4_miroc <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_midMIROC, args="outputformat=cloglog")
plot(P_cloglog_mid4_miroc)
writeRaster(P_cloglog_mid4_miroc, "~/path/to/dir/pungens_cloglog_midHolocene_miroc.asc")

## MPI-ESM-P
setwd("~/path/to/dir/WorldClim_data/memidbi_2-5m/")
bio2mp_mid <- raster("memidbi2.tif")
bio10mp_mid <- raster("memidbi10.tif")
bio11mp_mid <- raster("memidbi11.tif")
bio15mp_mid <- raster("memidbi15.tif")
bio17mp_mid <- raster("memidbi17.tif")

stacked_mid_mpi <- stack(bio2mp_mid, bio10mp_mid, bio11mp_mid, bio15mp_mid, bio17mp_mid)
crs(stacked_mid_mpi) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_mpi_mid <- crop(stacked_mid_mpi, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1]) 
unstacked <- unstack(cropped_mpi_mid)

setwd("~/path/to/dir/WorldClim_data/cropped_midMPI_climate_XlgExtent/")
writeRaster(unstacked[[1]], "~/Documents/TGG_SDM/WorldClim_data/cropped_midMPI_climate_XlgExtent/Bio2.asc")
writeRaster(unstacked[[2]], "~/Documents/TGG_SDM/WorldClim_data/cropped_midMPI_climate_XlgExtent/Bio10.asc")
writeRaster(unstacked[[3]], "~/Documents/TGG_SDM/WorldClim_data/cropped_midMPI_climate_XlgExtent/Bio11.asc")
writeRaster(unstacked[[4]], "~/Documents/TGG_SDM/WorldClim_data/cropped_midMPI_climate_XlgExtent/Bio15.asc")
writeRaster(unstacked[[5]], "~/Documents/TGG_SDM/WorldClim_data/cropped_midMPI_climate_XlgExtent/Bio16.asc")

setwd("~/path/to/dir/WorldClim_data/memidbi_2-5m/")
bio2mp <- raster("Bio2.asc")
bio10mp <- raster("Bio10.asc")
bio11mp <- raster("Bio11.asc")
bio15mp <- raster("Bio15.asc")
bio17mp <- raster("Bio17.asc")

env_stack_midMPI <- stack(bio2mp, bio10mp, bio11mp, bio15mp, bio17mp)
crs(env_stack_midMPI) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
env_stack_midMPI <- crop(env_stack_midMPI, ext)

P_raw_mid4_mpi <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_midMPI, args="outputformat=raw")
plot(P_raw_mid4_mpi)
writeRaster(P_raw_mid4_mpi, "~/path/to/dir/pungens_raw_midHolocene_mpi.asc")


#### Mid Holocene ensemble ####
avg_env_mid <- (env_stack_midCCSM + env_stack_midMIROC + env_stack_midMPI)/3

P_raw_mid_avg <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], avg_env_mid, args="outputformat=raw")
plot(P_raw_mid_avg)
writeRaster(P_raw_mid_avg, "~/path/to/dir/pungens_raw_midHolocene_ensemble.asc", overwrite=TRUE)


#### LGM ####
## CCSM4
setwd("~/path/to/dir/WorldClim_data/cclgmbi_2-5m/")
bio2c_lgm <- raster("cclgmbi2.tif")
bio10c_lgm <- raster("cclgmbi10.tif")
bio11c_lgm <- raster("cclgmbi11.tif")
bio15c_lgm <- raster("cclgmbi15.tif")
bio17c_lgm <- raster("cclgmbi17.tif")

stacked_lgm <- stack(bio2c_lgm, bio10c_lgm, bio11c_lgm, bio15c_lgm, bio17c_lgm)
crs(stacked_lgm) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_all <- crop(stacked_lgm, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1]) 
unstacked <- unstack(cropped_all)

setwd("~/path/to/dir/WorldClim_data/cropped_lgmCCSM_climate_XlgExtent/")
writeRaster(unstacked[[1]], "~/path/to/dir/WorldClim_data/cclgmbi_2-5m/Bio2.asc")
writeRaster(unstacked[[2]], "~/path/to/dir/WorldClim_data/cclgmbi_2-5m/Bio10.asc")
writeRaster(unstacked[[3]], "~/path/to/dir/WorldClim_data/cclgmbi_2-5m/Bio11.asc")
writeRaster(unstacked[[4]], "~/path/to/dir/WorldClim_data/cclgmbi_2-5m/Bio15.asc")
writeRaster(unstacked[[5]], "~/path/to/dir/WorldClim_data/cclgmbi_2-5m/Bio17.asc")


setwd("~/path/to/dir/WorldClim_data/cclgmbi_2-5m/")
bio2cc <- raster("Bio2.asc")
bio10cc <- raster("Bio10.asc")
bio11cc <- raster("Bio11.asc")
bio15cc <- raster("Bio15.asc")
bio17cc <- raster("Bio17.asc")

env_stack_lgmCCSM <- stack(bio2cc, bio10cc, bio11cc, bio15cc, bio17cc)

crs(env_stack_lgmCCSM) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
env_stack_lgmCCSM <- crop(env_stack_lgmCCSM, ext)


P_raw_lgm4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lgmCCSM, args="outputformat=raw")
plot(P_raw_lgm4)
writeRaster(P_raw_lgm4, "~/path/to/dir/pungens_raw_lgm_ccsm.asc", overwrite=TRUE)
P_cloglog_lgm4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lgmCCSM, args="outputformat=cloglog")
plot(P_cloglog_lgm4)
writeRaster(P_cloglog_lgm4, "~/path/to/dir/pungens_cloglog_lgm_ccsm.asc", overwrite=TRUE)

## MIROC-ESM
setwd("~/path/to/dir/WorldClim_data/mrlgmbi_2-5m")

bio2m_lgm <- raster("mrlgmbi2.tif")
bio10m_lgm <- raster("mrlgmbi10.tif")
bio11m_lgm <- raster("mrlgmbi11.tif")
bio15m_lgm <- raster("mrlgmbi15.tif")
bio17m_lgm <- raster("mrlgmbi17.tif")

stacked_lgm_m <- stack(bio2m_lgm, bio10m_lgm, bio11m_lgm, bio15m_lgm, bio17m_lgm)
crs(stacked_lgm_m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_all <- crop(stacked_lgm_m, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1]) 
unstacked <- unstack(cropped_all)

writeRaster(unstacked[[1]], "~/path/to/dir/WorldClim_data/cropped_lgmMIROC_climate_XlgExtent/Bio2.asc")
writeRaster(unstacked[[2]], "~/path/to/dir/WorldClim_data/cropped_lgmMIROC_climate_XlgExtent/Bio10.asc")
writeRaster(unstacked[[3]], "~/path/to/dir/WorldClim_data/cropped_lgmMIROC_climate_XlgExtent/Bio11.asc")
writeRaster(unstacked[[4]], "~/path/to/dir/WorldClim_data/cropped_lgmMIROC_climate_XlgExtent/Bio15.asc")
writeRaster(unstacked[[5]], "~/path/to/dir/WorldClim_data/cropped_lgmMIROC_climate_XlgExtent/Bio17.asc")


setwd("~/path/to/dir/WorldClim_data/cropped_lgmMIROC_climate_XlgExtent/")
bio2m <- raster("Bio2.asc")
bio10m <- raster("Bio10.asc")
bio11m <- raster("Bio11.asc")
bio15m <- raster("Bio15.asc")
bio17m <- raster("Bio17.asc")

env_stack_lgmMIROC <- stack(bio2m, bio10m, bio11m, bio15m, bio17m)
crs(env_stack_lgmMIROC) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"




P_raw_lgm4_miroc <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lgmMIROC, args="outputformat=raw")
plot(P_raw_lgm4_miroc)
writeRaster(P_raw_lgm4_miroc, "~/path/to/dir/pungens_raw_lgm_miroc.asc")
P_cloglog_lgm4_miroc <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lgmMIROC, args="outputformat=cloglog")
plot(P_cloglog_lgm4_miroc)
writeRaster(P_cloglog_lgm4_miroc, "~/path/to/dir/pungens_cloglog_lgm_miroc.asc")


## MPI-ESM-P
setwd("~/path/to/dir/WorldClim_data/melgmbi_2-5m/")
bio2mp_lgm <- raster("melgmbi2.tif")
bio10mp_lgm <- raster("melgmbi10.tif")
bio11mp_lgm <- raster("melgmbi11.tif")
bio15mp_lgm <- raster("melgmbi15.tif")
bio17mp_lgm <- raster("melgmbi17.tif")

stacked_lgm_mpi <- stack(bio2mp_lgm, bio10mp_lgm, bio11mp_lgm, bio15mp_lgm, bio17mp_lgm)
crs(stacked_lgm_mpi) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_all <- crop(stacked_lgm_mpi, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1]) 
unstacked <- unstack(cropped_all)

writeRaster(unstacked[[1]], "~/path/to/dir/WorldClim_data/melgmbi_2-5m/Bio2.asc")
writeRaster(unstacked[[2]], "~/path/to/dir/WorldClim_data/melgmbi_2-5m/Bio10.asc")
writeRaster(unstacked[[3]], "~/path/to/dir/WorldClim_data/melgmbi_2-5m/Bio11.asc")
writeRaster(unstacked[[4]], "~/path/to/dir/WorldClim_data/melgmbi_2-5m/Bio15.asc")
writeRaster(unstacked[[5]], "~/path/to/dir/WorldClim_data/melgmbi_2-5m/Bio17.asc")

bio2mp <- raster("Bio2.asc")
bio10mp <- raster("Bio10.asc")
bio11mp <- raster("Bio11.asc")
bio15mp <- raster("Bio15.asc")
bio17mp <- raster("Bio17.asc")

env_stack_lgmMPI <- stack(bio2mp, bio10mp, bio11mp, bio15mp, bio17mp)
crs(env_stack_lgmMPI) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

P_raw_lgm4_mpi <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lgmMPI, args="outputformat=raw")
plot(P_raw_lgm4_mpi)
writeRaster(P_raw_lgm4_mpi, "~/path/to/pungens_raw_lgm_mpi.asc")


#### LGM ensemble - raster mean ####
avg_env_lgm <- (env_stack_lgmCCSM + env_stack_lgmMIROC + env_stack_lgmMPI)/3
P_raw_lgm_avg <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], avg_env_lgm, args="outputformat=raw")
plot(P_raw_lgm_avg)
writeRaster(P_raw_lgm_avg, "~/path/to/pungens_raw_lgm_ensemble.asc", overwrite=TRUE)

plot(env_stack_lgmCCSM$Bio15)
plot(env_stack_lgmMIROC$Bio15)
plot(env_stack_lgmMPI$Bio15)
plot(avg_env_lgm$Bio15)

## Hindcast to Last Interglacial 
setwd("~/path/to/WorldClim_data/lig_30s_bio/")

bio2l <- raster("lig_30s_bio_2.bil")
lig2_crop <- crop(bio2l, ext)
bio10l <- raster("lig_30s_bio_10.bil")
bio11l <- raster("lig_30s_bio_11.bil")
bio15l <- raster("lig_30s_bio_15.bil")
bio17l <- raster("lig_30s_bio_17.bil")
plot(lig2_crop)

stacked <- stack(bio2l, bio10l, bio11l, bio15l, bio17l)
crs(stacked) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ext <- extent(-95, -60, 25, 55)
cropped_all <- crop(stacked, ext)
plot(cropped_all[[1]], main=names(cropped_all)[1])
setwd("~/path/to/WorldClim_data/cropped_lig_climate_XlgExtent/")

unstacked <- unstack(cropped_all)

writeRaster(unstacked[[1]], "~/Documents/TGG_SDM/WorldClim_data/cropped_lig_climate_XlgExtent/Bio2.asc")
writeRaster(unstacked[[2]], "~/Documents/TGG_SDM/WorldClim_data/cropped_lig_climate_XlgExtent/Bio10.asc")
writeRaster(unstacked[[3]], "~/Documents/TGG_SDM/WorldClim_data/cropped_lig_climate_XlgExtent/Bio11.asc")
writeRaster(unstacked[[4]], "~/Documents/TGG_SDM/WorldClim_data/cropped_lig_climate_XlgExtent/Bio15.asc")
writeRaster(unstacked[[5]], "~/Documents/TGG_SDM/WorldClim_data/cropped_lig_climate_XlgExtent/Bio17.asc")

bio2l <- raster("Bio2.asc")
plot(bio2l)
bio10l <- raster("Bio10.asc")
bio11l <- raster("Bio11.asc")
bio15l <- raster("Bio15.asc")
bio17l <- raster("Bio17.asc")

bio2a <- aggregate(bio2l, fact=5)
bio10a <- aggregate(bio10l, fact=5)
bio11a <- aggregate(bio11l, fact=5)
bio15a <- aggregate(bio15l, fact=5)
bio17a <- aggregate(bio17l, fact=5)


env_stack_lig <- stack(bio2a, bio10a, bio11a, bio15a, bio17a)
crs(env_stack_lig) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

P_raw_lig4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lig, args="outputformat=raw")
plot(P_raw_lig4)
writeRaster(P_raw_lig4, "~/path/to/pungens_raw_lig.asc")
P_cloglog_lig4 <- predict(eval.models(P_train)[["fc.LQ_rm.1"]], env_stack_lig, args="outputformat=cloglog")
plot(P_cloglog_lig4)
writeRaster(P_cloglog_lig4, "~/path/to/dir/pungens_cloglog_lig.asc")






#### P rigida hindcasts ####
setwd("~/path/to/dir/location_data")

R_xy <- read.csv("thinned_Rdata_thin1.csv")
R_loc <- R_xy[,2:3] 
View(R_loc)
# delete (-88.16528, 38.77611) (-84.00000, 32.00000) (-76.00981,35.13000)
R_loc <- R_xy[,2:3]  
index <- which(R_loc$Longitude <(-88))
index 
R_loc <- R_loc[-index,]
View(R_loc)
index2 <-which(R_loc$Latitude == (32.00000))
index2
R_loc <- R_loc[-index2,]
index3 <-which(R_loc$Latitude == (35.13000))
index3
R_loc <- R_loc[-index3,]

Longitude <- as.numeric(R_loc[,1])
Longitude
Latitude<- as.numeric(R_loc[,2])
yz <- cbind(Longitude, Latitude)
yz


bb.r<- bbox(yz)
summary(bb.r)
bb.buf.r <- extent(bb.r[1]-5,bb.r[3]+5, bb.r[2]-5,bb.r[4]+5)  
summary(bb.buf.r)
envs.backg.r <- crop(env_stack_now,bb.buf.r)  
plot(envs.backg.r[[1]], main=names(envs.backg.r)[1])
points(yz)
envs.backg.r

bg.r <- randomPoints(envs.backg.r[[1]], n=10000)
bg.r <- as.data.frame(bg.r)
plot(envs.backg.r[[1]], legend=FALSE)
points(bg.r, col='red')
head(bg.r)
colnames(bg.r) <- c("Longitude", "Latitude")


block.r <- get.block(yz, bg.r, orientation = "lat_lon")
table(block.r$occs.grp)

evalplot.grps(pts = yz, pts.grp = block.r$occs.grp, envs = envs.backg.r) + 
  ggplot2::ggtitle("Spatial block partitions: occurrences")

evalplot.grps(pts = bg.r, pts.grp = block.r$bg.grp, envs = envs.backg.r) + 
  ggplot2::ggtitle("Spatial block partitions: background")


occs.zR <- cbind(yz, raster::extract(env_stack_now, yz))
bg.zR <- cbind(bg.r, raster::extract(env_stack_now, bg.r))

tune.args <- list(fc = c("L","LQ","LQH","H"), rm = 1:5)
ps <- list(orientation = "lat_lon")

R_train5 <- ENMevaluate(occ=yz, env=env_stack_now, bg=bg.r, tune.args=tune.args, partitions = "block", other.settings= os, partition.settings = ps, doClamp=TRUE, parallel = FALSE, algorithm = 'maxent.jar',overlap = TRUE)
pred <- R_train5@predictions[[which(R_train5@results$delta.AICc==0)]] #LQH_3
plot(pred)
pred
m5.mxR <- eval.models(R_train5)[["fc.LQH_rm.3"]]
m5.mxR #save html as Maxent_model_rigida_Delta0
plot(m5.mxR, type = "raw")
# saved as Rplot_R_train5_var_contribution_delta0.pdf

eval.predictions(R_train5)

res_R <- eval.results(R_train5)
opt.aicc_R <- res_R$delta.AICc==0
opt.aicc_R



setwd("~/path/to/dir/")
write.csv(bg.r, "Prigida_ENM_background_pts.csv")
write.csv(R_train5@results,file="Prigida_ENM_eval_results.csv")

#### null model ####
mod.null.r <- ENMnulls(R_train5, mod.settings = list(fc = "LQH", rm = 3), no.iter = 100)
null.emp.results(mod.null.r)
null_emp_stats_r <- null.emp.results(mod.null.r)
write.csv(null_emp_stats_r, "null_versus_empirical_models_R_train5.csv")
evalplot.nulls(mod.null.r, stats = c("or.10p", "auc.val"), plot.type = "histogram")
# saved as Rplot_null_versus_empirical_models_rigida_R_train5.pdf

#### metadata for Rigida model ####
rmm_R <- eval.rmm(R_train5)
rangeModelMetadata::rmmToCSV(rmm_R, "rmm_metadata_Rigida_R_train5.csv")


## write occs.gp files and bg.grp files for ridida and pungens


dismo::response(eval.models(R_train5)[["fc.LQH_rm.3"]])
# saved as Rplot_R_train5_response_curves_delta0
dismo::response(m5.mxR, range="p")
dismo::response(m5.mxR, range="pa")

write.csv(block$occs.grp, "~/path/to/dir/Pungens_occs_grp_SDMs.csv")
write.csv(block$bg.grp, "~/path/to/dir/Pungens_bg_grp_SDMs.csv")
write.csv(block.r$occs.grp, "~/path/to/dir/rigida_occs_grp_SDMs.csv")
write.csv(block.r$bg.grp, "~/path/to/dir/rigida_bg_grp_SDMs.csv")

block.r$bg.grp

R_raw_now <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_now, args="outputformat=raw")
plot(R_raw_now)
writeRaster(R_raw_now, "~/path/to/dir/rigida_raw_current.asc")
R_cloglog_now <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_now, args="outputformat=cloglog")
plot(R_cloglog_now)
writeRaster(R_cloglog_now, "~/path/to/dir/rigida_cloglog_current.asc")


## mid-Holocene
## CCSM
R_raw_mid <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_midCCSM, args="outputformat=raw")
plot(R_raw_mid)
writeRaster(R_raw_mid, "~/path/to/dir/rigida_raw_midHolocene.asc")
plot(P_raw_mid4)
R_cloglog_mid <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_midCCSM, args="outputformat=cloglog")
plot(R_cloglog_mid)
writeRaster(R_cloglog_mid, "~/path/to/dir/rigida_cloglog_midHolocene.asc")

## MIROC
R_raw_midMiroc <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_midMIROC, args="outputformat=raw")
plot(R_raw_midMiroc)
plot(P_raw_mid4_miroc)
writeRaster(R_raw_midMiroc, "~/Documents/TGG_SDM/Corrected_Model/Updated_script_models/rigida_raw_midHolocene_miroc.asc")
R_cloglog_midMiroc <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_midMIROC, args="outputformat=cloglog")
plot(R_cloglog_midMiroc)
writeRaster(R_cloglog_midMiroc, "~/path/to/dir/rigida_cloglog_midHolocene_miroc.asc")

## MPI
R_raw_midMPI <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_midMPI, args="outputformat=raw")
plot(R_raw_midMPI)
writeRaster(R_raw_midMPI, "~/path/to/dir/rigida_raw_midHolocene_mpi.asc")

plot(P_raw_mid4_mpi)

#### Rigida midH ensemble ####
R_raw_mid_avg <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], avg_env_mid, args="outputformat=raw")
plot(R_raw_mid_avg)
plot(P_raw_mid_avg)
writeRaster(R_raw_mid_avg, "~/path/to/dir/rigida_raw_midHolocene_ensemble.asc", overwrite=TRUE)



## LGM
## CCSM
R_raw_lgm <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lgmCCSM, args="outputformat=raw")
plot(R_raw_lgm)
writeRaster(R_raw_lgm, "~/path/to/dir/rigida_raw_lgm.asc", overwrite=TRUE)
R_cloglog_lgm <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lgmCCSM, args="outputformat=cloglog")
plot(R_cloglog_lgm)
writeRaster(R_cloglog_lgm, "~/path/to/dir/rigida_cloglog_lgm.asc", overwrite=TRUE)
plot(P_cloglog_lgm4)
## MIROC
R_raw_lgmMiroc <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lgmMIROC, args="outputformat=raw")
plot(R_raw_lgmMiroc)
writeRaster(R_raw_lgmMiroc, "~/path/to/dir/rigida_raw_lgm_miroc.asc")
plot(P_raw_lgm4_miroc)
R_cloglog_lgmMiroc <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lgmMIROC, args="outputformat=cloglog")
plot(R_cloglog_lgmMiroc)
writeRaster(R_cloglog_lgmMiroc, "~/path/to/dir/rigida_cloglog_lgm_miroc.asc")
## MPI
R_raw_lgmMPI <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lgmMPI, args="outputformat=raw")
plot(R_raw_lgmMPI)
writeRaster(R_raw_lgmMPI, "~/path/to/dir/rigida_raw_lgm_mpi.asc")


#### Rigida lgm ensemble ####
R_raw_lgm_avg <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], avg_env_lgm, args="outputformat=raw")
plot(R_raw_lgm_avg)
plot(P_raw_lgm_avg)
writeRaster(R_raw_lgm_avg, "~/path/to/dir/rigida_raw_lgm_ensemble.asc", overwrite=TRUE)



## LIG

R_raw_lig <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lig, args="outputformat=raw")
plot(R_raw_lig)
writeRaster(R_raw_lig, "~/path/to/dir/rigida_raw_lig.asc")
R_cloglog_lig <- predict(eval.models(R_train5)[["fc.LQH_rm.3"]], env_stack_lig, args="outputformat=cloglog")
plot(R_cloglog_lig)
writeRaster(R_cloglog_lig, "~/path/to/dir/rigida_cloglog_lig.asc")
plot(P_raw_lig4)

#### standardize rasters in ENMTools ####
library(ENMTools)

## current SDMs
P_now_standard <- raster.standardize(P_raw_now4)
plot(P_now_standard)
plot(P_raw_now4)

R_now_standard <- raster.standardize(R_raw_now)
plot(R_now_standard)
plot(R_raw_now)

## mid-Holocene
P_mid_standard_miroc <- raster.standardize(P_raw_mid4_miroc)
plot(P_mid_standard_miroc)
P_mid_standard_ccsm <- raster.standardize(P_raw_mid4)
plot(P_mid_standard_ccsm)
P_mid_standard_mpi <- raster.standardize(P_raw_mid4_mpi)
plot(P_mid_standard_mpi)
P_mid_ens <- raster.standardize(P_raw_mid_avg)
writeRaster(P_mid_ens, "~/path/to/dir/standardized_SDM_predictions/Pungens_midH_standardized_ensemble.asc", overwrite=TRUE)


R_mid_standard_miroc <- raster.standardize(R_raw_midMiroc)
plot(R_mid_standard_miroc)
R_mid_standard_ccsm <- raster.standardize(R_raw_mid)
plot(R_mid_standard_ccsm)
R_mid_standard_mpi <- raster.standardize(R_raw_midMPI)
plot(P_mid_standard_mpi)
R_mid_ens <- raster.standardize(R_raw_mid_avg)
plot(R_mid_ens)
writeRaster(R_mid_ens, "~/path/to/dir/standardized_SDM_predictions/Rigida_midH_standardized_ensemble.asc", overwrite=TRUE)


P_lgm_standard_miroc <- raster.standardize(P_raw_lgm4_miroc)
plot(P_lgm_standard_miroc)
P_lgm_standard_ccsm <- raster.standardize(P_raw_lgm4)
plot(P_lgm_standard_ccsm)
P_lgm_standard_mpi <- raster.standardize(P_raw_lgm4_mpi)
plot(P_lgm_standard_mpi)

P_lgm_standard_ens <- raster.standardize(P_raw_lgm_avg)
plot(P_lgm_standard_ens)
writeRaster(P_lgm_standard_ens, "~/path/to/dir/standardized_SDM_predictions/Pungens_lgm_standardized_ensemble.asc", overwrite=TRUE)


R_lgm_standard_miroc <- raster.standardize(R_raw_lgmMiroc)
plot(R_lgm_standard_miroc)
R_lgm_standard_ccsm <- raster.standardize(R_raw_lgm)
plot(R_lgm_standard_ccsm)
R_lgm_standard_mpi <- raster.standardize(R_raw_lgmMPI)
plot(R_lgm_standard_mpi)
R_lgm_standard_ens <- raster.standardize(R_raw_lgm_avg)
plot(R_lgm_standard_ens)
writeRaster(R_lgm_standard_ens, "~/path/to/dir/standardized_SDM_predictions/Rigida_lgm_standardized_ensemble.asc", overwrite=TRUE)



P_lig_standard <- raster.standardize(P_raw_lig4)
plot(P_lig_standard)

R_lig_standard <- raster.standardize(R_raw_lig)
plot(R_lig_standard)

setwd("~/path/to/dir/standardized_SDM_predictions/")
writeRaster(P_now_standard, "Pungens_nowSDM_standardized_ENMTools.asc" , overwrite=TRUE)
writeRaster(P_mid_standard_ccsm, "Pungens_midHolocene_ccsmSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(P_mid_standard_miroc, "Pungens_midHolocene_mirocSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(P_mid_standard_mpi, "Pungens_midHolocene_mpiSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(P_lgm_standard_ccsm, "Pungens_lgm_ccsmSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(P_lgm_standard_miroc, "Pungens_lgm_mirocSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(P_lgm_standard_mpi, "Pungens_lgm_mpiSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(P_lig_standard, "Pungens_lig_SDM_standardized_ENMTools.asc", overwrite=TRUE)



writeRaster(R_now_standard, "Rigida_nowSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_mid_standard_ccsm, "Rigida_midHolocene_ccsmSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_mid_standard_miroc, "Rigida_midHolocene_mirocSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_mid_standard_mpi, "Rigida_midHolocene_mpiSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_lgm_standard_ccsm, "Rigida_lgm_ccsmSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_lgm_standard_miroc, "Rigida_lgm_mirocSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_lgm_standard_mpi, "Rigida_lgm_mpiSDM_standardized_ENMTools.asc", overwrite=TRUE)
writeRaster(R_lig_standard, "Rigida_lig_SDM_standardized_ENMTools.asc", overwrite=TRUE)


#### cumulative transformation of each SDM prediction ####
library(scales)
cumulative = function(r){
  r = round(1000000000 * r) # do raster *10^9 and round it to integer
  v = values(r) # get raster values as vector
  t = table(v) # get the number of pixels with each value
  t = data.frame(t) # make a data frame (pixel value, frequency)
  t$v = as.numeric(as.character(t$v)) # change the pixel values from factor to numeric
  t$product = t$v * t$Freq # multiply pixel values by their frequencies
  t$cumsum = cumsum(t$product) # cumsum: sum all pixel values <= than itself
  t$rescaled = rescale(t$cumsum, to = c(0,1)) # rescale the sum to 0-100
  r = subs(r, t, by="v", which="rescaled") # reclassify raster
  return(r)
}

P_now_cum <- cumulative(P_now_standard)
plot(P_now_cum)
plot(P_raw_now4)

R_now_cum <- cumulative(R_now_standard)
plot(R_now_cum)

P_mid_cum1 <- cumulative(P_mid_standard_ccsm)
plot(P_mid_cum1)
P_mid_cum2 <- cumulative(P_mid_standard_miroc)
plot(P_mid_cum2)
P_mid_cum3 <- cumulative(P_mid_standard_mpi)
plot(P_mid_cum3)
P_mid_ensemble <- cumulative(P_mid_ens)
plot(P_mid_ensemble)

R_mid_cum1 <- cumulative(R_mid_standard_ccsm)
plot(R_mid_cum1)
R_mid_cum2 <- cumulative(R_mid_standard_miroc)
plot(R_mid_cum2)
R_mid_cum3 <- cumulative(R_mid_standard_mpi)
plot(R_mid_cum3)
R_mid_ensemble <- cumulative(R_mid_ens)
plot(R_mid_ensemble)

P_lgm_cum1 <- cumulative(P_lgm_standard_ccsm)
plot(P_lgm_cum1)
P_lgm_cum2 <- cumulative(P_lgm_standard_miroc)
plot(P_lgm_cum2)
P_lgm_cum3 <- cumulative(P_lgm_standard_mpi)
plot(P_lgm_cum3)
P_lgm_ens_cum <- cumulative(P_lgm_standard_ens)
plot(P_lgm_ens_cum)

R_lgm_cum1 <- cumulative(R_lgm_standard_ccsm)
plot(R_lgm_cum1)
R_lgm_cum2 <- cumulative(R_lgm_standard_miroc)
plot(R_lgm_cum2)
R_lgm_cum3 <- cumulative(R_lgm_standard_mpi)
plot(R_lgm_cum3)
R_lgm_ens_cum <- cumulative(R_lgm_standard_ens)
plot(R_lgm_ens_cum)
plot(P_lgm_ens_cum)

P_lig_cum <- cumulative(P_lig_standard)
plot(P_lig_cum)

R_lig_cum <- cumulative(R_lig_standard)
plot(R_lig_cum)

setwd("~/path/to/dir/cumulative_SDM_predictions/")
writeRaster(P_now_cum, "Pungens_nowSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_mid_cum1, "Pungens_midHolocene_ccsmSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_mid_cum2, "Pungens_midHolocene_mirocSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_mid_cum3, "Pungens_midHolocene_mpiSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_lgm_cum1, "Pungens_lgm_ccsmSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_lgm_cum2, "Pungens_lgm_mirocSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_lgm_cum3, "Pungens_lgm_mpiSDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_lig_cum, "Pungens_lig_SDM_cumulative.asc", overwrite=TRUE)
writeRaster(P_lgm_ens_cum, "Pungens_lgm_standard_ensemble_cumulative.asc", overwrite=TRUE)
writeRaster(P_mid_ensemble, "Pungens_midH_standard_ensemble_cumulative.asc", overwrite=TRUE)


writeRaster(R_now_cum, "Rigida_nowSDM_cumulative.asc")
writeRaster(R_mid_cum1, "Rigida_midHolocene_ccsmSDM_cumulative.asc", overwrite=TRUE)
writeRaster(R_mid_cum2, "Rigida_midHolocene_mirocSDM_cumulative.asc")
writeRaster(R_mid_cum3, "Rigida_midHolocene_mpiSDM_cumulative.asc", overwrite=TRUE)
writeRaster(R_lgm_cum1, "Rigida_lgm_ccsmSDM_cumulative.asc", overwrite=TRUE)
writeRaster(R_lgm_cum2, "Rigida_lgm_mirocSDM_cumulative.asc")
writeRaster(R_lgm_cum3, "Rigida_lgm_mpiSDM_cumulative.asc", overwrite=TRUE)
writeRaster(R_lig_cum, "Rigida_lig_SDM_cumulative.asc")
writeRaster(R_lgm_ens_cum, "Rigida_lgm_standard_ensemble_cumulative.asc", overwrite=TRUE)
writeRaster(R_mid_ensemble, "Rigida_midH_standard_ensemble_cumulative.asc", overwrite=TRUE)



#### plot cumulative SDM rasters ####

setwd("~/path/to/dir/location_data")

R_xy <- read.csv("thinned_Rdata_thin1.csv")
R_loc <- R_xy[,2:3] 
View(R_loc)
# delete (-88.16528, 38.77611) (-84.00000, 32.00000) (-76.00981,35.13000)
R_loc <- R_xy[,2:3]  
index <- which(R_loc$Longitude <(-88))
index 
R_loc <- R_loc[-index,]
View(R_loc)
index2 <-which(R_loc$Latitude == (32.00000))
index2
R_loc <- R_loc[-index2,]
index3 <-which(R_loc$Latitude == (35.13000))
index3
R_loc <- R_loc[-index3,]

Longitude <- as.numeric(R_loc[,1])
Latitude<- as.numeric(R_loc[,2])
yz <- cbind(Longitude, Latitude)
yz

## pungens occurrence data
P_xy <- read.csv("thinned_Pdata_thin1.csv")
P_loc <- P_xy[,2:3]
Longitude <- as.numeric(P_loc[,1])
Latitude<- as.numeric(P_loc[,2])
xy <- cbind(Longitude, Latitude)


library(rgdal)
library(RColorBrewer)
library(raster)
setwd("~/path/to/dir/ice_extent_LGM_shapefile/")
## shapefiles from Dyke, 2003 .. see email with github link to A. Wickert
ice <- readOGR("ice018000.shp")
crs(ice) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

setwd("~/path/to/dir/cumulative_SDM_predictions/")
ext_plot <- extent(-95, -60, 28, 55)
P_lig_cum <- raster("Pungens_lig_SDM_cumulative.asc")
crs(P_lig_cum) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_lig_crop <- crop(P_lig_cum, ext_plot)
plot(P_lig_crop)
points(xy, pch=20, col="gray33", cex=0.5)


R_lig_cum <- raster("Rigida_lig_SDM_cumulative.asc")
crs(R_lig_cum) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_lig_crop <- crop(R_lig_cum, ext_plot)
plot(R_lig_crop)
points(yz, pch=20, col="gray33", cex=0.5)



## LGM plots

P_lgm_cum1c <- raster("Pungens_lgm_ccsmSDM_cumulative.asc")
crs(P_lgm_cum1c) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_lgm_crop1c <- crop(P_lgm_cum1c, ext_plot)
plot(P_lgm_crop1c, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)

P_lgm_cum1m <- raster("Pungens_lgm_mirocSDM_cumulative.asc")
crs(P_lgm_cum1m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_lgm_crop1m <- crop(P_lgm_cum1m, ext_plot)
plot(P_lgm_crop1m, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)

P_lgm_cum1mp <- raster("Pungens_lgm_mpiSDM_cumulative.asc")
crs(P_lgm_cum1mp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_lgm_crop1mp <- crop(P_lgm_cum1mp, ext_plot)
plot(P_lgm_crop1mp, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)


P_lgm_ens_cum <- raster("Pungens_lgm_standard_ensemble_cumulative.asc")
crs(P_lgm_ens_cum) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_lgm_crop3 <- crop(P_lgm_ens_cum, ext_plot)
plot(P_lgm_crop3, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)

R_lgm_cum1c <- raster("Rigida_lgm_ccsmSDM_cumulative.asc")
crs(R_lgm_cum1c) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_lgm_crop1c <- crop(R_lgm_cum1c, ext_plot)
plot(R_lgm_crop1c)
points(yz, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)

R_lgm_cum1m <- raster("Rigida_lgm_mirocSDM_cumulative.asc")
plot(R_lgm_cum1m)
crs(R_lgm_cum1m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_lgm_crop1m <- crop(R_lgm_cum1m, ext_plot)
plot(R_lgm_crop1m)
points(yz, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)

R_lgm_cum1mp <- raster("Rigida_lgm_mpiSDM_cumulative.asc")
crs(R_lgm_cum1mp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_lgm_crop1mp <- crop(R_lgm_cum1mp, ext_plot)
plot(R_lgm_crop1mp)
points(yz, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)

R_lgm_ens_cum <- raster("Rigida_lgm_standard_ensemble_cumulative.asc")
crs(R_lgm_ens_cum) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_lgm_crop3 <- crop(R_lgm_ens_cum, ext_plot)
plot(R_lgm_crop3)
points(yz, pch=20, col="gray33", cex=0.5)
plot(ice, add=TRUE, border="lightblue", lwd=1.5)


## Mid Holocene plots
P_mid_cum1c <- raster("Pungens_midHolocene_ccsmSDM_cumulative.asc")
crs(P_mid_cum1c) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_mid_crop1c <- crop(P_mid_cum1c, ext_plot)
plot(P_mid_crop1c, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)


P_mid_cum1m <- raster("Pungens_midHolocene_mirocSDM_cumulative.asc")
crs(P_mid_cum1m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_mid_crop1m <- crop(P_mid_cum1m, ext_plot)
plot(P_mid_crop1m, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)

P_mid_cum1mp <- raster("Pungens_midHolocene_mpiSDM_cumulative.asc")
crs(P_mid_cum1mp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_mid_crop1mp <- crop(P_mid_cum1mp, ext_plot)
plot(P_mid_crop1mp, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)

P_mid_ens <- raster("Pungens_midH_standard_ensemble_cumulative.asc")
crs(P_mid_ens) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_mid_crop3 <- crop(P_mid_ens, ext_plot)
plot(P_mid_crop3, ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)


R_mid_cum1c <- raster("Rigida_midHolocene_ccsmSDM_cumulative.asc")
crs(R_mid_cum1c) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_mid_crop1c <- crop(R_mid_cum1c, ext_plot)
plot(R_mid_crop1c)
points(yz, pch=20, col="gray33", cex=0.5)


R_mid_cum1m <- raster("Rigida_midHolocene_mirocSDM_cumulative.asc")
crs(R_mid_cum1m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_mid_crop1m <- crop(R_mid_cum1m, ext_plot)
plot(R_mid_crop1m)
points(yz, pch=20, col="gray33", cex=0.5)

R_mid_cum1mp <- raster("Rigida_midHolocene_mpiSDM_cumulative.asc")
crs(R_mid_cum1mp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_mid_crop1mp <- crop(R_mid_cum1mp, ext_plot)
plot(R_mid_crop1mp)
points(yz, pch=20, col="gray33", cex=0.5)


R_mid_ens <- raster("Rigida_midH_standard_ensemble_cumulative.asc")
crs(R_mid_ens) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_mid_crop3 <- crop(R_mid_ens, ext_plot)
plot(R_mid_crop3)
points(yz, pch=20, col="gray33", cex=0.5)

plot(P_mid_crop3)

## current plots

P_now_cum <- raster("Pungens_nowSDM_cumulative.asc")
crs(P_now_cum) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
P_now_crop <- crop(P_now_cum, ext_plot)
plot(P_now_crop, xlab="Longitude", ylab="Latitude")
points(xy, pch=20, col="gray33", cex=0.5)


R_now_cum <- raster("Rigida_nowSDM_cumulative.asc")
crs(R_now_cum) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
R_now_crop <- crop(R_now_cum, ext_plot)
plot(R_now_crop, xlab="Longitude")
points(yz, pch=20, col="gray33", cex=0.5)




#### Form a hypothesis for pop size change  ####
P_now_pts <- rasterToPoints(P_now_cum, fun=function(P_now_cum){P_now_cum>0.5}, spatial = TRUE)
P_mid_pts <- rasterToPoints(P_mid_ensemble, fun=function(P_mid_ensemble){P_mid_ensemble>0.5}, spatial = TRUE)
P_lgm_pts <- rasterToPoints(P_lgm_ens_cum, fun=function(P_lgm_ens_cum){P_lgm_ens_cum>0.5}, spatial = TRUE)
P_lig_pts <- rasterToPoints(P_lig_cum, fun=function(P_lig_cum){P_lig_cum>0.5}, spatial = TRUE)


R_now_pts <- rasterToPoints(R_now_cum, fun=function(R_now_cum){R_now_cum>0.5}, spatial = TRUE)
R_mid_pts <- rasterToPoints(R_mid_ensemble, fun=function(R_mid_ensemble){R_mid_ensemble>0.5}, spatial = TRUE)
R_lgm_pts <- rasterToPoints(R_lgm_ens_cum, fun=function(R_lgm_ens_cum){R_lgm_ens_cum>0.5}, spatial = TRUE)
R_lig_pts <- rasterToPoints(R_lig_cum, fun=function(R_lig_cum){R_lig_cum>0.5}, spatial = TRUE)


dim(P_now_pts) #6632
dim(P_mid_pts) #4924
dim(P_lgm_pts) #2947
dim(P_lig_pts) # 3577

dim(R_now_pts) #11128 
dim(R_mid_pts) #3725
dim(R_lgm_pts) #4400
dim(R_lig_pts) # 1983





#### Form a hypothesis for gene flow ####
## make venn diagram of overlap ##
install.packages("VennDiagram")   # Install & load VennDiagram package
library("VennDiagram")
library(dplyr)



#LIG SDM overlap
R_lig50_df <- as.data.frame(R_lig_pts)
dim(R_lig50_df) #1983
P_lig50_df <- as.data.frame(P_lig_pts)
dim(P_lig50_df) #3577
common_coords_lig <- inner_join(R_lig50_df[,2:3], P_lig50_df[,2:3])
dim(common_coords_lig) #640

grid.newpage()                    # Create new plotting page
lig_pts <- draw.pairwise.venn(area1 = 3577,    # Draw pairwise venn diagram
                              area2 = 1983,
                              cross.area = 640, col=c("blue", "orange"), fill=c("blue", "orange"), alpha=c(0.25, 0.25), euler.d = TRUE, scaled=TRUE)




# lgm SDM overlap... from ensembled predictions
R_lgm50_df <- as.data.frame(R_lgm_pts)
dim(R_lgm50_df) #4400
P_lgm50_df <- as.data.frame(P_lgm_pts)
dim(P_lgm50_df) #2947
common_coords_lgm <- inner_join(R_lgm50_df[,2:3], P_lgm50_df[,2:3])
dim(common_coords_lgm) #606

grid.newpage()
lgm_ens_pts <- draw.pairwise.venn(area1 = 2947,    # Draw pairwise venn diagram
                                  area2 = 4400 ,
                                  cross.area = 606, col=c("blue", "orange"), fill=c("blue", "orange"), alpha=c(0.25, 0.25), euler.d = TRUE, scaled=TRUE)



#mid-Holocene SDM overlap ... ensemble from raw predictions
R_mid50_df <- as.data.frame(R_mid_pts)
dim(R_mid50_df) #3725
P_mid50_df <- as.data.frame(P_mid_pts)
dim(P_mid50_df) #4924
common_coords_mid <- inner_join(R_mid50_df[,2:3], P_mid50_df[,2:3])
dim(common_coords_mid) #1412

grid.newpage()
midH_ens_pts <- draw.pairwise.venn(area1 = 3725,    # Draw pairwise venn diagram
                                   area2 = 4924,
                                   cross.area = 1412, col=c("blue", "orange"), fill=c("blue", "orange"), alpha=c(0.25, 0.25), euler.d = TRUE, scaled=TRUE)


#current SDM overlap
R_now50_df <- as.data.frame(R_now_pts)
dim(R_now50_df) #11128
P_now50_df <- as.data.frame(P_now_pts)
dim(P_now50_df) #6632
common_coords_now <- inner_join(R_now50_df[,2:3], P_now50_df[,2:3])
dim(common_coords_now) #2498

grid.newpage()
current_pts <- draw.pairwise.venn(area1 = 6632,    # Draw pairwise venn diagram
                                  area2 = 11128,
                                  cross.area = 2498, col=c("blue", "orange"), fill=c("blue", "orange"), alpha=c(0.25, 0.25), euler.d = TRUE, scaled=TRUE)





#### raster overlap in ENMTools ####

now_overlap <- raster.overlap(P_now_cum, R_now_cum)
now_overlap
##$D
#[1] 0.6124016

##$I
#[1] 0.8763198

##$rank.cor
#[1] 0.9092813

mid_overlap <- raster.overlap(P_mid_ensemble, R_mid_ensemble)
mid_overlap
#$D
#[1] 0.4249661

#$I
#[1] 0.6882606

#$rank.cor
#[1] 0.5410139



lgm_overlap_ens <- raster.overlap(P_lgm_ens_cum, R_lgm_ens_cum)
lgm_overlap_ens
#$D
#[1] 0.1700541

#$I
#[1] 0.4175774

#$rank.cor
#[1] 0.712554


lig_overlap <- raster.overlap(P_lig_cum, R_lig_cum)
lig_overlap
#$D
#[1] 0.2879556

#$I
#[1] 0.5277512

#$rank.cor
#[1] 0.5580765








#### Raster overlap of GCMs for LGM and Mid-Holocene ####
lgm_models_ccsm_miroc <- raster.overlap(P_lgm_cum1, P_lgm_cum2)
lgm_models_ccsm_miroc
#$D
#[1] 0.05240278

#$I
#[1] 0.1261928

#$rank.cor
#[1] 0.3636107


lgm_models_ccsm_mpi <- raster.overlap(P_lgm_cum1, P_lgm_cum3)
lgm_models_ccsm_mpi
#$D
#[1] 0.1237215

#$I
#[1] 0.2782551

#$rank.cor
#[1] 0.3975498

lgm_models_miroc_mpi <- raster.overlap(P_lgm_cum2, P_lgm_cum3)
lgm_models_miroc_mpi
#$D
#[1] 0.02486987

#$I
#[1] 0.09501053

#$rank.cor
#[1] 0.2307604

midH_models_ccsm_miroc <- raster.overlap(P_mid_cum1, P_mid_cum2)
midH_models_ccsm_miroc
#$D
#[1] 0.3302237

#$I
#[1] 0.5703702

#$rank.cor
#[1] 0.9025507

midH_models_ccsm_mpi <- raster.overlap(P_mid_cum1, P_mid_cum3)
midH_models_ccsm_mpi
#$D
#[1] 0.5045245

#$I
#[1] 0.7212972

#$rank.cor
#[1] 0.8678716

midH_models_miroc_mpi <- raster.overlap(P_mid_cum2, P_mid_cum3)
midH_models_miroc_mpi
#$D
#[1] 0.4009482

#$I
#[1] 0.6181731

#$rank.cor
#[1] 0.8597167



