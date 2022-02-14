#####
#
## create fuels map
#
#####

rm(list=ls())

# wd inherited from project

# load libraries
library(raster) # version 3.0-2
library(sp) # version 1.3-1
library(rgeos) # version 0.5-1
library(rgdal) # version 1.4-4
library(sf) # version 0.7-7
library(dplyr) # version 0.8.3
library(ggplot2) # version 3.2.1

####
# 1a.  load and mosaic lidar and naip rasters
####

# ### individual tiles not included in data deposit, so this
# # section commented out. final mosaic of predictor variables are included
# 
# ### all lidar metrics
# 
# path="processed_data/GRTE_rasters/lidar_tiles"
# 
# r.list = list.files(path, pattern=".tif", full.names = TRUE)
# 
# # read in first raster
# r.first = raster(r.list[1])
# r.out = stack(r.first, extent(r.first)-60) # crop it
# 
# for(i in c(2:length(r.list))) {
# 
#   r.in = stack(r.list[i])
#   
#   # skip if raster is less than 3 by 3
#   if(nrow(r.in)<3 | ncol(r.in)<3) {
#     next
#   }
#   
#   # remove outermost edge to omit cells that may have been half-in, half-out in underlying data
#   # crop by 60 m to remove cells on each side
#   r.crop = crop(r.in, extent(r.in)-60)
#   
#   # merge with master raster
#   r.out = merge(r.out, r.in)
#   
#   print(paste("done with ",i))
# 
# }
# 
# # names for future reference
# # names(r.out) = c("zmax","zmean","zcv","zp10","zp25","zp50","zp75","zp90","zp99","zfcc","zfsc")  
# 
# writeRaster(r.out, "processed_data/GRTE_rasters/grte_30m_lidar_full.tif", overwrite=TRUE)
# 
# 
# path="processed_data/GRTE_rasters/naip_tiles"
# 
# # functionalize this
# 
# mosaic_naip = function(var.in) {
#   
#   r.list = list.files(path, pattern=var.in, full.names = TRUE)
#   
#   # read in first raster
#   r.first = raster(r.list[1])
#   r.out = stack(r.first, extent(r.first)-60) # crop it
#   
#   for(i in c(2:length(r.list))) {
#     
#     r.in = stack(r.list[i])
#     
#     # skip if raster is less than 3 by 3
#     if(nrow(r.in)<3 | ncol(r.in)<3) {
#       next
#     }
#     
#     # remove outermost edge to omit cells that may have been half-in, half-out in underlying data
#     # crop by 60 min to remove cells on each side
#     r.crop = crop(r.in, extent(r.in)-60)
#     
#     # merge with master raster
#     r.out = merge(r.out, r.in)
#     
#     print(paste("done with ",i))
#     
#   }
#   
#   writeRaster(r.out, paste0("processed_data/GRTE_rasters/grte_30m_naip_",var.in,".tif"))
#   
#   return(r.out)
# 
# }
# 
# 
# # naip: max, min, mean, sd, hom
# 
# mosaic_naip("max")
# mosaic_naip("min")
# mosaic_naip("mean")
# mosaic_naip("sd")
# mosaic_naip("hom")
# 
# 
# ### naip: cv
# 
# # read in mean and sd
# naip.sd = stack("processed_data/GRTE_rasters/grte_30m_naip_sd.tif")
# naip.mean = stack("processed_data/GRTE_rasters/grte_30m_naip_mean.tif")
# 
# naip.cv = (naip.sd/naip.mean) * 100
# 
# writeRaster(naip.cv, "processed_data/GRTE_rasters/grte_30m_naip_cv.tif")

####
# 1b.  create final mask using study area boundary
####

### read in fuels map mask, created in las_statistics_extraction
# masked to conifer, deciduous, or shrubland
# removed disturbances, management after 2014

grte.maskin = raster("processed_data/GRTE_rasters/grte_30m_final_mask_nad2011.tif")

### read in extent shapefile

full.shp = st_read("processed_data/GRTE_shps/fuels_map_outline.shp")

# create final mask using shapefile

grte.mask = mask(crop(grte.maskin, as(full.shp, "Spatial"), snap="out"),as(full.shp,"Spatial"))

####
# 2. check lidar predictors for fuels map creation
####

### read in plot data
fuels.mets = read.csv("processed_data/fuels_map_variables/lidar_metrics.csv")

### read in and mask lidar

lidar.full = stack("processed_data/GRTE_rasters/grte_30m_lidar_full.tif")

# names
names(lidar.full) = c("zmax","zmean","zcv","zp10","zp25","zp50","zp75","zp90","zp99","zfcc","zfsc")

# crop and mask
lidar.mask = mask(crop(lidar.full,grte.mask),grte.mask)

# cleaning
# look at histograms and compare to plot data
# omit cells that should not be included in analysis
# note any strange values
# set problematic values to NA

# zmax
plot(lidar.mask[[1]])
par(mfrow=c(1,2)) 
hist(lidar.mask[[1]])
hist(fuels.mets$zmax) 
summary(lidar.mask[[1]])
# looks good except for outliers > 60 m in height
# remove these cells from analysis
length(lidar.mask[[1]][lidar.mask[[1]] > 60])  # 34 cells
lidar.out = lidar.mask[[1]]
lidar.out[lidar.out>60] = NA

# update lidar mask
lidar.mask2 = mask(lidar.mask, lidar.out)

# zmax check
par(mfrow=c(1,1)) 
plot(lidar.mask2[[1]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[1]])
hist(fuels.mets$zmax) 
summary(lidar.mask2[[1]])
# looks good 

# zmean
par(mfrow=c(1,1)) 
plot(lidar.mask2[[2]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[2]])
hist(fuels.mets$zmean) 
summary(lidar.mask2[[2]])
# looks good

# zcv
par(mfrow=c(1,1)) 
plot(lidar.mask2[[3]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[3]])
hist(fuels.mets$zcv) 
summary(lidar.mask2[[3]])
length(lidar.mask2[[3]][lidar.mask2[[3]] > 250])  # 5406 cells, keep an eye on anything that needs cv
# some additional NAs, remove these from analysis
lidar.mask1 = lidar.mask2
lidar.mask2 = mask(lidar.mask1,lidar.mask1[[3]])

# zp10
par(mfrow=c(1,1)) 
plot(lidar.mask2[[4]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[4]])
hist(fuels.mets$zp10) 
summary(lidar.mask2[[4]])
length(lidar.mask[[4]][lidar.mask[[4]] < 0])  # 30,088 cells.
# no top models use zp10

# zp25
par(mfrow=c(1,1)) 
plot(lidar.mask2[[5]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[5]])
hist(fuels.mets$zp25) 
summary(lidar.mask2[[5]])
length(lidar.mask[[5]][lidar.mask[[5]] < 0])  # 193 cells, very few
# no top models use zp25

# zp50
par(mfrow=c(1,1)) 
plot(lidar.mask2[[6]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[6]])
hist(fuels.mets$zp50) 
summary(lidar.mask2[[6]])
length(lidar.mask[[6]][lidar.mask[[6]] < 0])  # 1 cell
# some equations use zp50

# zp75
par(mfrow=c(1,1)) 
plot(lidar.mask2[[7]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[7]])
hist(fuels.mets$zp75) 
summary(lidar.mask2[[7]])
# looks good

# zp90
par(mfrow=c(1,1)) 
plot(lidar.mask2[[8]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[8]])
hist(fuels.mets$zp90) 
summary(lidar.mask2[[8]])
# looks good

# zp99
par(mfrow=c(1,1)) 
plot(lidar.mask2[[9]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[9]])
hist(fuels.mets$zp99) 
summary(lidar.mask2[[9]])
# looks good

# zfcc
par(mfrow=c(1,1)) 
plot(lidar.mask2[[10]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[10]])
hist(fuels.mets$zfcc) 
summary(lidar.mask2[[10]])
# definitely some issues with areas where length(Z) returns
# a very large number instead of 0 when no canopy first returns
length(lidar.mask2[[10]][lidar.mask2[[10]] > 100]) # 242328
length(lidar.mask2[[10]][lidar.mask2[[10]] == 100]) # 0
# set problematic cells to 0
lidar.mask2[[10]][lidar.mask2[[10]] > 100] = 0
# reexamine
par(mfrow=c(1,1)) 
plot(lidar.mask2[[10]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[10]])
hist(fuels.mets$zfcc) 
summary(lidar.mask2[[10]])
# looks good

# zfsc
par(mfrow=c(1,1)) 
plot(lidar.mask2[[11]])
par(mfrow=c(1,2)) 
hist(lidar.mask2[[11]])
hist(fuels.mets$zfsc) 
summary(lidar.mask2[[11]])
# looks pretty good, take closer look at 100%
length(lidar.mask2[[11]][lidar.mask2[[11]] > 100]) # 0
length(lidar.mask2[[11]][lidar.mask2[[11]] == 100]) # 2, leaving as is


### read in and mask naip

naip.max = stack("processed_data/GRTE_rasters/grte_30m_naip_max.tif")
naip.min = stack("processed_data/GRTE_rasters/grte_30m_naip_min.tif")
naip.mean = stack("processed_data/GRTE_rasters/grte_30m_naip_mean.tif")
naip.sd = stack("processed_data/GRTE_rasters/grte_30m_naip_sd.tif")
naip.cv = stack("processed_data/GRTE_rasters/grte_30m_naip_cv.tif")
naip.hom = stack("processed_data/GRTE_rasters/grte_30m_naip_hom.tif")

# names
names(naip.max) = c("Rmax","Gmax","Bmax","NIRmax","NDVImax")
names(naip.min) = c("Rmin","Gmin","Bmin","NIRmin","NDVImin")
names(naip.mean) = c("Rmean","Gmean","Bmean","NIRmean","NDVImean")
names(naip.sd) = c("Rsd","Gsd","Bsd","NIRsd","NDVIsd")
names(naip.cv) = c("Rcv","Gcv","Bcv","NIRcv","NDVIcv")
names(naip.hom) = c("Rhom","Ghom","Bhom","NIRhom","NDVIhom")

# combine
naip.full = stack(crop(naip.max,naip.hom),
                  crop(naip.min,naip.hom),
                  crop(naip.mean,naip.hom),
                  crop(naip.sd,naip.hom),
                  crop(naip.cv,naip.hom),
                  naip.hom)


# crop and mask
naip.mask1 = mask(crop(naip.full,grte.mask),grte.mask)

# also mask based on updated lidar
naip.mask = mask(naip.mask1,lidar.mask2[[1]])

length(lidar.mask2[[1]][!is.na(lidar.mask2[[3]])])
length(naip.mask[[1]][!is.na(naip.mask[[1]])])

# cleaning
# omit cells that should not be included in analysis
# note any strange values
# set problematic values to NA

fuels.naip=read.csv("processed_data/fuels_map_variables/naip_metrics.csv")

# max
plot(naip.mask[[c(1:5)]])
par(mfrow=c(2,2)) 
hist(naip.mask[[c(1:2)]])
hist(fuels.naip$Rmax)
hist(fuels.naip$Gmax)
hist(naip.mask[[c(3:4)]])
hist(fuels.naip$Bmax)
hist(fuels.naip$NIRmax)
hist(naip.mask[[5]])
hist(fuels.naip$NDVImax)
summary(naip.mask[[c(1:5)]])
# looks good except for NDVI outliers
# remove these cells from analysis
length(naip.mask[[5]][naip.mask[[5]] > 1])  # 47 cells, remove these from analysis
# looks like these are high-elevation, issues with snow
length(naip.mask[[5]][naip.mask[[5]] < 0])  # 111 cells, leaving these in for now
naip.out = naip.mask[[5]]
naip.out[naip.out>1] = NA
# mask to this
naip.mask2 = mask(naip.mask,naip.out)
# check ndvi now
hist(naip.mask2[[5]])
hist(fuels.naip$NDVImax) # looks good

# min
plot(naip.mask2[[c(6:10)]])
par(mfrow=c(2,2)) 
hist(naip.mask2[[c(6:7)]])
hist(fuels.naip$Rmin)
hist(fuels.naip$Gmin)
hist(naip.mask2[[c(8:9)]])
hist(fuels.naip$Bmin)
hist(fuels.naip$NIRmin)
hist(naip.mask2[[10]])
hist(fuels.naip$NDVImin)
summary(naip.mask2[[c(6:10)]])
# looks good 

# mean
plot(naip.mask2[[c(11:15)]])
par(mfrow=c(2,2)) 
hist(naip.mask2[[c(11:12)]])
hist(fuels.naip$Rmean)
hist(fuels.naip$Gmean)
hist(naip.mask2[[c(13:14)]])
hist(fuels.naip$Bmean)
hist(fuels.naip$NIRmean)
hist(naip.mask2[[15]])
hist(fuels.naip$NDVImean)
summary(naip.mask2[[c(11:15)]])
# looks good 

# sd
plot(naip.mask2[[c(16:20)]])
par(mfrow=c(2,2)) 
hist(naip.mask2[[c(16:17)]])
hist(fuels.naip$Rsd)
hist(fuels.naip$Gsd)
hist(naip.mask2[[c(18:19)]])
hist(fuels.naip$Bsd)
hist(fuels.naip$NIRsd)
hist(naip.mask2[[20]])
hist(fuels.naip$NDVIsd)
summary(naip.mask2[[c(16:20)]])
# looks good 

# cv
plot(naip.mask2[[c(21:25)]])
par(mfrow=c(2,2)) 
hist(naip.mask2[[c(21:22)]])
hist(fuels.naip$Rcv)
hist(fuels.naip$Gcv)
hist(naip.mask2[[c(23:24)]])
hist(fuels.naip$Bcv)
hist(fuels.naip$NIRcv)
hist(naip.mask2[[25]])
hist(fuels.naip$NDVIcv)
# both look a little crazy, no top models use NDVIcv
length(naip.mask2[[25]][naip.mask[[25]] > 1000])  # 15,334 cells
length(naip.mask2[[25]][naip.mask[[25]] < -500])  # 28,203 cells
summary(naip.mask2[[c(21:25)]])
# looks good 

# hom
plot(naip.mask2[[c(26:30)]])
par(mfrow=c(2,2)) 
hist(naip.mask2[[c(26:27)]])
hist(fuels.naip$Rhom)
hist(fuels.naip$Ghom)
hist(naip.mask2[[c(28:29)]])
hist(fuels.naip$Bhom)
hist(fuels.naip$NIRhom)
hist(naip.mask2[[30]])
hist(fuels.naip$NDVIhom)
summary(naip.mask2[[c(26:30)]])
# looks good 

# finally re-mask lidar with new naip.mask2
lidar.mask3 = mask(lidar.mask2,naip.mask2[[1]])

length(lidar.mask3[[1]][!is.na(lidar.mask3[[1]])]) # 846099
length(naip.mask2[[1]][!is.na(naip.mask2[[1]])]) # 846099

# export final versions
writeRaster(lidar.mask3, "processed_data/GRTE_rasters/grte_30m_lidar_masked.tif", overwrite=TRUE)
writeRaster(naip.mask2, "processed_data/GRTE_rasters/grte_30m_naip_masked.tif", overwrite=TRUE)

#### prep vegetation maps

# read in
grte.con = raster("processed_data/GRTE_rasters/grte_30m_veg_con_nad2011.tif")
grte.dec = raster("processed_data/GRTE_rasters/grte_30m_veg_dec_nad2011.tif")
grte.shr = raster("processed_data/GRTE_rasters/grte_30m_veg_shr_nad2011.tif")

# stack
grte.veg = stack(grte.con,grte.dec,grte.shr)
names(grte.veg) = c("conifer","deciduous","shrubland")

# mask
veg.mask = mask(crop(grte.veg,lidar.mask3[[1]]),lidar.mask3[[1]])

# need to binarize
# if majority forest, set to 1
# first create forested relative to shrubland
veg.mask[[4]] = (veg.mask[[1]] + veg.mask[[2]] + 0.01)/(veg.mask[[3]] + 0.01)
veg.mask[[4]][veg.mask[[4]] >= 1] = 1
veg.mask[[4]][veg.mask[[4]] < 1] = 0
# use this to subset conifer and deciduous
veg.mask[[1]] = veg.mask[[1]] * veg.mask[[4]]
veg.mask[[2]] = veg.mask[[2]] * veg.mask[[4]]
# use inverse to subset shrubland
veg.mask[[3]] = veg.mask[[3]] * (1-veg.mask[[4]])

# look at output
plot(veg.mask)
summary(veg.mask)

# set to 0 or 1
# anything still classified as shrubland gets assigned to shrubland
length(veg.mask[[3]][veg.mask[[3]]<0.5 & veg.mask[[3]]>0]) # 532
veg.mask[[3]][veg.mask[[3]] > 0] = 1
# conifer or deciduous is based on which is greater, ties go to conifer
# divide, if >= 1, assign conifer
veg.mask[[5]] = (veg.mask[[1]]+0.01) / (veg.mask[[2]] + 0.01)
veg.mask[[5]][veg.mask[[5]]>=1] = 1
veg.mask[[5]][veg.mask[[5]]<1] = 0
# multiply by conifer and deciuous
veg.mask[[1]] = veg.mask[[1]] * veg.mask[[5]]
veg.mask[[2]] = veg.mask[[2]] * (1-veg.mask[[5]])
# set to 0 or 1
veg.mask[[1]][veg.mask[[1]]>0] = 1
veg.mask[[2]][veg.mask[[2]]>0] = 1


# subset back to conifer, deciduous and shrubland
veg.final = veg.mask[[c(1:3)]]

names(veg.final) = c("conifer","deciduous","shrubland")
plot(veg.final)

length(veg.final[[1]][veg.final[[1]] == 1]) # 505240
length(veg.final[[2]][veg.final[[2]] == 1]) # 28812
length(veg.final[[3]][veg.final[[3]] == 1]) # 312075

505240+28812+312075
# match

# save this too
writeRaster(veg.final,"processed_data/GRTE_rasters/grte_30m_veg_masked.tif")

####
# 3. fuels map
####

### read in final predictors if necessary

lidar.map = stack("processed_data/GRTE_rasters/grte_30m_lidar_masked.tif")
naip.map = stack("processed_data/GRTE_rasters/grte_30m_naip_masked.tif")
veg.map = stack("processed_data/GRTE_rasters/grte_30m_veg_masked.tif")

names(lidar.map) = c("zmax","zmean","zcv","zp10","zp25","zp50","zp75","zp90","zp99","zfcc","zfsc")
names(naip.map) = c("Rmax","Gmax","Bmax","NIRmax","NDVImax",
                    "Rmin","Gmin","Bmin","NIRmin","NDVImin",
                    "Rmean","Gmean","Bmean","NIRmean","NDVImean",
                    "Rsd","Gsd","Bsd","NIRsd","NDVIsd",
                    "Rcv","Gcv","Bcv","NIRcv","NDVIcv",
                    "Rhom","Ghom","Bhom","NIRhom","NDVIhom")
full.map = stack(lidar.map,naip.map)
names(veg.map) = c("conifer","deciduous","shrubland")

### to create unmasked version of fuels map for landfire comparison, load full predictor set
# commented out because this is for comparison purposes only
# full.map = stack(lidar.full,mask(crop(naip.full,lidar.full[[1]]),lidar.full[[1]]))

### read final selected equations
model.fin = read.csv("analysis/fuels_prediction_map/model_final_selection.csv")
model.pred = read.csv("analysis/fuels_prediction_map/model_final_predictors.csv")

# make veg map into a mask
veg.map[veg.map==0] = NA
# add forest mask
veg.map[["forest"]] = max(veg.map[["conifer"]],veg.map[["deciduous"]],na.rm=TRUE)
plot(veg.map)


# go through all models

model.in = model.fin

# create separate raster for each vegetation subset
conifer.out = raster()
forest.out = raster()
shrub.out = raster()

### uncomment to create unmasked map of CC, CH, CBD, CBH for comparison with field data
# full.map = stack(lidar.full,mask(crop(naip.full,lidar.full[[1]]),lidar.full[[1]]))

# doing this as for loop to track progress
for(i in 1:dim(model.in)[1]) {
  print(i)
  
  model.form = model.pred %>%
    filter(fuels_metric == model.in[i,]$fuels_metric,
           veg_types == model.in[i,]$veg_types,
           model_number == model.in[i,]$model_number) %>%
    dplyr::select(c(pred,coef))
  
  # set up equation for each possible number of predictors
  
  if (model.in[i,]$n_pred == 1) {
    
    print("1 predictor")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,")"))
    
    fuelsvar = coefi + (full.map[[pred1]] * coef1) 
    
  } else if (model.in[i,]$n_pred == 2) {
    
    print("2 predictors")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    pred2 = as.character(model.form[3,]$pred)
    coef2 = round(as.numeric(model.form[3,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,") + (",coef2,pred2,")"))
    
    fuelsvar = coefi + (full.map[[pred1]] * coef1) + (full.map[[pred2]] * coef2) 
    
  } else if (model.in[i,]$n_pred == 3) {
    
    print("3 predictors")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    pred2 = as.character(model.form[3,]$pred)
    coef2 = round(as.numeric(model.form[3,]$coef),4)
    pred3 = as.character(model.form[4,]$pred)
    coef3 = round(as.numeric(model.form[4,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,") + (",coef2,pred2,") + (",coef3,pred3,")"))
    
    fuelsvar = coefi + (full.map[[pred1]] * coef1) + (full.map[[pred2]] * coef2) + (full.map[[pred3]] * coef3)
    
  } else if (model.in[i,]$n_pred == 4) {
    
    print("4 predictors")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    pred2 = as.character(model.form[3,]$pred)
    coef2 = round(as.numeric(model.form[3,]$coef),4)
    pred3 = as.character(model.form[4,]$pred)
    coef3 = round(as.numeric(model.form[4,]$coef),4)
    pred4 = as.character(model.form[5,]$pred)
    coef4 = round(as.numeric(model.form[5,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,") + (",coef2,pred2,") + (",coef3,pred3,") + (",coef4,pred4,")"))
    
    fuelsvar = coefi + (full.map[[pred1]] * coef1) + (full.map[[pred2]] * coef2) + (full.map[[pred3]] * coef3) + (full.map[[pred4]] * coef4)
    
  } else if (model.in[i,]$n_pred == 5) {
    
    print("5 predictors")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    pred2 = as.character(model.form[3,]$pred)
    coef2 = round(as.numeric(model.form[3,]$coef),4)
    pred3 = as.character(model.form[4,]$pred)
    coef3 = round(as.numeric(model.form[4,]$coef),4)
    pred4 = as.character(model.form[5,]$pred)
    coef4 = round(as.numeric(model.form[5,]$coef),4)
    pred5 = as.character(model.form[6,]$pred)
    coef5 = round(as.numeric(model.form[6,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,") + (",coef2,pred2,") + (",coef3,pred3,") + (",coef4,pred4,") + (",coef5,pred5,")"))
    
    fuelsvar = coefi + (full.map[[pred1]] * coef1) + (full.map[[pred2]] * coef2) + (full.map[[pred3]] * coef3) + (full.map[[pred4]] * coef4) + (full.map[[pred5]]*coef5)
    
  }
  
  names(fuelsvar) = model.in[i,]$fuels_metric
  
  print(model.in[i,]$fuels_metric)
  print(model.in[i,]$veg_types)
  
  # stack with correct raster
  if (model.in[i,]$veg_types == "Conifer") {
    conifer.out = stack(conifer.out,fuelsvar)
  } else if (model.in[i,]$veg_types == "Conifer_deciduous") {
    forest.out = stack(forest.out,fuelsvar)
  } else if (model.in[i,]$veg_types == "Shrubland") {
    shrub.out = stack(shrub.out,fuelsvar)
  }
  
}


# check forest vs. conifer equations
length(forest.out[[3]][!is.na(forest.out[[1]])])
length(conifer.out[[1]][!is.na(conifer.out[[1]])])

# mask
conifer.mask = mask(conifer.out,veg.map[["conifer"]])
forest.mask = mask(forest.out, veg.map[["forest"]])
shrub.mask = mask(shrub.out, veg.map[["shrubland"]])

# check forest vs. conifer equations
length(forest.mask[[3]][!is.na(forest.mask[[1]])])
length(conifer.mask[[1]][!is.na(conifer.mask[[1]])])


# transform
conifer.mask[["CBD_kg_m3"]] = conifer.mask[["CBD_kg_m3"]]^3
conifer.mask[["CBH_m"]] = conifer.mask[["CBH_m"]]^2
conifer.mask[["CFL_kg_m2"]] = conifer.mask[["CFL_kg_m2"]]^2
conifer.mask[["CWD_Mg_ha"]] = exp(conifer.mask[["CWD_Mg_ha"]])

forest.mask[["BA_m2_ha"]] = forest.mask[["BA_m2_ha"]]^2
forest.mask[["TPH_largetrees"]] = forest.mask[["TPH_largetrees"]]^2
forest.mask[["TPH_trees"]] = exp(forest.mask[["TPH_trees"]])

shrub.mask[["Shrub_pctDead"]] = exp(shrub.mask[["Shrub_pctDead"]])

# examine - note outlier values
# read in fuels for comparison
fuels.comb = read.csv("processed_data/fuels_map_variables/field_plot_2019_fuels.csv") 

par(mfrow=c(1,2))
hist(conifer.mask[["CBD_kg_m3"]])
hist(fuels.comb$CBD_kg_m3)
hist(conifer.mask[["CBH_m"]])
hist(fuels.comb$CBH_m)
hist(conifer.mask[["CFL_kg_m2"]])
hist(fuels.comb$CFL_kg_m2)
hist(conifer.mask[["CWD_cover_pct"]])
hist(fuels.comb$CWD_cover_pct)
hist(conifer.mask[["CWD_Mg_ha"]])
hist(fuels.comb$CWD_Mg_ha)
hist(conifer.mask[["Shrub_cover_pct"]])
hist(fuels.comb$Shrub_cover_pct)
hist(conifer.mask[["Shrub_ht_m"]])
hist(fuels.comb$Shrub_ht_m)

# outliers 
# set min or max to natural value OR approx. twice the value of max observed in the field
conifer.mod = conifer.mask
conifer.mod[["CBD_kg_m3"]] = reclassify(conifer.mod[["CBD_kg_m3"]],c(-Inf,0,0,0.4,Inf,0.4))
conifer.mod[["CBH_m"]] = reclassify(conifer.mod[["CBH_m"]],c(-Inf,0,0,5,Inf,5))
conifer.mod[["CFL_kg_m2"]] = reclassify(conifer.mod[["CFL_kg_m2"]],c(-Inf,0,0,4,Inf,4))
conifer.mod[["CWD_cover_pct"]] = reclassify(conifer.mod[["CWD_cover_pct"]],c(-Inf,0,0,50,Inf,50))
conifer.mod[["CWD_Mg_ha"]] = reclassify(conifer.mod[["CWD_Mg_ha"]],c(-Inf,0,0,300,Inf,300))
conifer.mod[["Shrub_cover_pct"]] = reclassify(conifer.mod[["Shrub_cover_pct"]],c(-Inf,0,0,100,Inf,100))
conifer.mod[["Shrub_ht_m"]] = reclassify(conifer.mod[["Shrub_ht_m"]],c(-Inf,0,0,2,Inf,2))


# check

hist(conifer.mod[["CBD_kg_m3"]], main = "fuels map CBD", xlab="CBD", xlim=c(0,0.4), breaks=20)
hist(fuels.comb[fuels.comb$Vegetation_type=="Conifer",]$CBD_kg_m3, main = "field observations CBD", xlab="CBD", xlim = c(0,0.4), breaks=10)
hist(conifer.mod[["CBH_m"]])
hist(fuels.comb$CBH_m)
hist(conifer.mod[["CFL_kg_m2"]])
hist(fuels.comb$CFL_kg_m2)
hist(conifer.mod[["CWD_cover_pct"]])
hist(fuels.comb$CWD_cover_pct)
hist(conifer.mod[["CWD_Mg_ha"]])
hist(fuels.comb$CWD_Mg_ha)
hist(conifer.mod[["Shrub_cover_pct"]])
hist(fuels.comb$Shrub_cover_pct)
hist(conifer.mod[["Shrub_ht_m"]])
hist(fuels.comb$Shrub_ht_m)


# forest

# examine
par(mfrow=c(1,2))
hist(forest.mask[["BA_m2_ha"]])
hist(fuels.comb$BA_m2_ha)
hist(forest.mask[["CC_pct"]])
hist(fuels.comb$CC_pct)
hist(forest.mask[["SH_m"]])
hist(fuels.comb$SH_m)
hist(forest.mask[["TPH_trees"]])
hist(fuels.comb$TPH_trees)


# set min or max to natural value OR approx. twice the value of max observed in the field
forest.mod = forest.mask
forest.mod[["CC_pct"]] = reclassify(forest.mod[["CC_pct"]],c(-Inf,0,0,100,Inf,100))
forest.mod[["SH_m"]] = reclassify(forest.mod[["SH_m"]],c(-Inf,0,0))
forest.mod[["TPH_trees"]] = reclassify(forest.mod[["TPH_trees"]],c(80000,Inf,80000))

# check
par(mfrow=c(1,2))
hist(forest.mod[["BA_m2_ha"]])
hist(fuels.comb$BA_m2_ha)
hist(forest.mod[["CC_pct"]])
hist(fuels.comb$CC_pct)
hist(forest.mod[["SH_m"]])
hist(fuels.comb$SH_m)
hist(forest.mod[["TPH_trees"]])
hist(fuels.comb$TPH_trees)

# shrubland

# examine
par(mfrow=c(1,2))
hist(shrub.mask[["prop_sagebrush"]])
hist(fuels.comb$Sagebrush_cover_pct/fuels.comb$Shrub_cover_pct)
hist(shrub.mask[["Shrub_cover_pct"]])
hist(fuels.comb$Shrub_cover_pct)
hist(shrub.mask[["Shrub_ht_m"]])
hist(fuels.comb$Shrub_ht_m)
hist(shrub.mask[["Shrub_pctDead"]])
hist(fuels.comb$Shrub_pctDead)

# set min or max to natural value OR approx. twice the value of max observed in the field
shrub.mod = shrub.mask
shrub.mod[["prop_sagebrush"]] = reclassify(shrub.mod[["prop_sagebrush"]],c(-Inf,0,0,1,Inf,1))
shrub.mod[["Shrub_cover_pct"]] = reclassify(shrub.mod[["Shrub_cover_pct"]],c(-Inf,0,0,100,Inf,100))
shrub.mod[["Shrub_ht_m"]] = reclassify(shrub.mod[["Shrub_ht_m"]],c(-Inf,0,0))
shrub.mod[["Shrub_pctDead"]] = reclassify(shrub.mod[["Shrub_pctDead"]],c(-Inf,0,0,100,Inf,100))

# check
par(mfrow=c(1,2))
hist(shrub.mod[["prop_sagebrush"]])
hist(fuels.comb$Sagebrush_cover_pct/fuels.comb$Shrub_cover_pct)
hist(shrub.mod[["Shrub_cover_pct"]])
hist(fuels.comb$Shrub_cover_pct)
hist(shrub.mod[["Shrub_ht_m"]])
hist(fuels.comb$Shrub_ht_m)
hist(shrub.mod[["Shrub_pctDead"]])
hist(fuels.comb$Shrub_pctDead)

### quantify outliers
# how many values had to be reclassified?
out.pct = data.frame(var = c("CBD_kg_m3","CBH_m","CFL_kg_m2","CWD_cover_pct","CWD_Mg_ha","Shrub_cover_pct","Shrub_ht_m","BA_m2_ha","CC_pct","SH_m","TPH_trees","prop_sagebrush","Shrub_cover_pct","Shrub_ht_m","Shrub_pctDead"),
                     veg_class = c(rep("conifer",7),rep("forest",4),rep("shrubland",4)),
                     below_0 = c(length(conifer.mask[["CBD_kg_m3"]][conifer.mask[["CBD_kg_m3"]]<0]),
                                 length(conifer.mask[["CBH_m"]][conifer.mask[["CBH_m"]]<0]),
                                 length(conifer.mask[["CFL_kg_m2"]][conifer.mask[["CFL_kg_m2"]]<0]),
                                 length(conifer.mask[["CWD_cover_pct"]][conifer.mask[["CWD_cover_pct"]]<0]),
                                 length(conifer.mask[["CWD_Mg_ha"]][conifer.mask[["CWD_Mg_ha"]]<0]),
                                 length(conifer.mask[["Shrub_cover_pct"]][conifer.mask[["Shrub_cover_pct"]]<0]),
                                 length(conifer.mask[["Shrub_ht_m"]][conifer.mask[["Shrub_ht_m"]]<0]),
                                 length(forest.mask[["BA_m2_ha"]][forest.mask[["BA_m2_ha"]]<0]),
                                 length(forest.mask[["CC_pct"]][forest.mask[["CC_pct"]]<0]),
                                 length(forest.mask[["SH_m"]][forest.mask[["SH_m"]]<0]),
                                 length(forest.mask[["TPH_trees"]][forest.mask[["TPH_trees"]]<0]),
                                 length(shrub.mask[["prop_sagebrush"]][shrub.mask[["prop_sagebrush"]]<0]),
                                 length(shrub.mask[["Shrub_cover_pct"]][shrub.mask[["Shrub_cover_pct"]]<0]),
                                 length(shrub.mask[["Shrub_ht_m"]][shrub.mask[["Shrub_ht_m"]]<0]),
                                 length(shrub.mask[["Shrub_pctDead"]][shrub.mask[["Shrub_pctDead"]]<0])),
                     above_max = c(length(conifer.mask[["CBD_kg_m3"]][conifer.mask[["CBD_kg_m3"]]>0.4]),
                                   length(conifer.mask[["CBH_m"]][conifer.mask[["CBH_m"]]>5]),
                                   length(conifer.mask[["CFL_kg_m2"]][conifer.mask[["CFL_kg_m2"]]>4]),
                                   length(conifer.mask[["CWD_cover_pct"]][conifer.mask[["CWD_cover_pct"]]>50]),
                                   length(conifer.mask[["CWD_Mg_ha"]][conifer.mask[["CWD_Mg_ha"]]>300]),
                                   length(conifer.mask[["Shrub_cover_pct"]][conifer.mask[["Shrub_cover_pct"]]>100]),
                                   length(conifer.mask[["Shrub_ht_m"]][conifer.mask[["Shrub_ht_m"]]>2]),
                                   length(forest.mask[["BA_m2_ha"]][forest.mask[["BA_m2_ha"]]>75]),
                                   length(forest.mask[["CC_pct"]][forest.mask[["CC_pct"]]>100]),
                                   length(forest.mask[["SH_m"]][forest.mask[["SH_m"]]>70]),
                                   length(forest.mask[["TPH_trees"]][forest.mask[["TPH_trees"]]>80000]),
                                   length(shrub.mask[["prop_sagebrush"]][shrub.mask[["prop_sagebrush"]]>1]),
                                   length(shrub.mask[["Shrub_cover_pct"]][shrub.mask[["Shrub_cover_pct"]]>100]),
                                   length(shrub.mask[["Shrub_ht_m"]][shrub.mask[["Shrub_ht_m"]]>3.2]),
                                   length(shrub.mask[["Shrub_pctDead"]][shrub.mask[["Shrub_pctDead"]]>100])),
                     tot_cells = c(rep(length(conifer.mask[[1]][!is.na(conifer.mask[[1]])]),7),
                                   rep(length(forest.mask[[1]][!is.na(forest.mask[[1]])]),4),
                                   rep(length(shrub.mask[[1]][!is.na(shrub.mask[[1]])]),4))) %>%
  mutate(below_0 = 100*below_0/tot_cells,
         above_max = 100*above_max/tot_cells,
         tot_out = below_0+above_max)

out.pct %>%
  filter(veg_class %in% c("forest","conifer")) %>%
  summary()

out.pct %>%
  filter(veg_class %in% c("shrubland")) %>%
  summary()

out.pct %>%
  summary()

### write final rasters

# conifer
writeRaster(conifer.mod[["CBD_kg_m3"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_CBD-kg-m3_nad2011.tif", overwrite=TRUE)
writeRaster(conifer.mod[["CBH_m"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_CBH-m_nad2011.tif", overwrite=TRUE)
writeRaster(conifer.mod[["CFL_kg_m2"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_CFL-kg-m2_nad2011.tif", overwrite=TRUE)
writeRaster(conifer.mod[["CWD_cover_pct"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_CWD-cover-pct_nad2011.tif", overwrite=TRUE)
writeRaster(conifer.mod[["CWD_Mg_ha"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_CWD-Mg-ha_nad2011.tif", overwrite=TRUE)
writeRaster(conifer.mod[["Shrub_cover_pct"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_Shrub-cover-pct_nad2011.tif", overwrite=TRUE)
writeRaster(conifer.mod[["Shrub_ht_m"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_conifer_Shrub-ht-m_nad2011.tif", overwrite=TRUE)

# forest
writeRaster(forest.mod[["BA_m2_ha"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_forest_BA-m2-ha_nad2011.tif", overwrite=TRUE)
writeRaster(forest.mod[["CC_pct"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_forest_CC-pct_nad2011.tif", overwrite=TRUE)
writeRaster(forest.mod[["SH_m"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_forest_CH-m_nad2011.tif", overwrite=TRUE)
writeRaster(forest.mod[["TPH_trees"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_forest_TPH-all-trees_nad2011.tif", overwrite=TRUE)

# shrubland
writeRaster(shrub.mod[["prop_sagebrush"]]*shrub.mod[["Shrub_cover_pct"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_shrubland_Sagebrush-cover-pct_nad2011.tif", overwrite=TRUE)
writeRaster(shrub.mod[["Shrub_cover_pct"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_shrubland_Shrub-cover-pct_nad2011.tif", overwrite=TRUE)
writeRaster(shrub.mod[["Shrub_ht_m"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_shrubland_Shrub-ht-m_nad2011.tif", overwrite=TRUE)
writeRaster(shrub.mod[["Shrub_pctDead"]],"processed_data/GRTE_rasters/final_fuels_map/grte_30m_shrubland_Shrub-dead-pct_nad2011.tif", overwrite=TRUE)

# # unmasked rasters for landfire comparison
# # note need to uncomment different version of full.map above
# writeRaster(conifer.mod[["CBD_kg_m3"]],"processed_data/GRTE_rasters/fuels_map_comparison/grte_30m_conifer_CBD-kg-m3_nad2011_unmasked.tif", overwrite=TRUE)
# writeRaster(conifer.mod[["CBH_m"]],"processed_data/GRTE_rasters/fuels_map_comparison/grte_30m_conifer_CBH-m_nad2011_unmasked.tif", overwrite=TRUE)
# writeRaster(forest.mod[["CC_pct"]],"processed_data/GRTE_rasters/fuels_map_comparison/grte_30m_forest_CC-pct_nad2011_unmasked.tif", overwrite=TRUE)
# writeRaster(forest.mod[["SH_m"]],"processed_data/GRTE_rasters/fuels_map_comparison/grte_30m_forest_CH-m_nad2011_unmasked.tif", overwrite=TRUE)