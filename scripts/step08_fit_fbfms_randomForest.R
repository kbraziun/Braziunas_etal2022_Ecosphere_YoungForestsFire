#####
#
## fit las metrics to fbfm characteristics with random forest
#
#####

rm(list=ls())

# load libraries
library(dplyr) # version 0.8.3
library(openxlsx) # version 4.1.3
library(sf) # version 0.7-7
library(tidyr) # version 1.0.2
library(caret) # version 6.0-84
library(ranger) # version 0.12.1
library(raster) # version 3.0-2

# plotting
library(ggplot2) # version 3.2.1
library(corrplot) # version 0.84
library(leaflet) # version 2.0.3
library(viridis) # version 0.5.1
library(RColorBrewer) # version 1.1-2
library(sjPlot) # version 2.8.6

####
# 1. load data
####

### 2019 field data

# fuels and plot locations
fuels.in = read.csv("processed_data/fuels_map_variables/field_plot_2019_fuels.csv") 
plots.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                     sheet = "General_plot_measurements",
                     colNames=TRUE)

fuels = plots.in %>%
  dplyr::select(c(Plot_code,Postprocessed_easting,Postprocessed_northing)) %>%
  left_join(fuels.in, by="Plot_code")


head(fuels)
summary(fuels)
summary(fuels$FBFM)


### reclassified 2002-3 fbfms

gfuels.in = read.xlsx("data/GRTE_Vegmap_Plots/FBFM_Reclassification/GRTE_FBFM_VegMap_2002-3COMPLETE.xlsx",
                     sheet = "Surface_FBFMs",
                     colNames=TRUE) %>%
  filter(!is.na(FBFM), FBFM != "TU")  # remove one plot without fully classified FBFM

gplots.in1 = st_read("data/GRTE_Vegmap_Plots/2002_Plots/Spatial_data/2002_plots_shpf/plots_2002.shp") %>%
  # remove prefix from plot code
  mutate(PlotNumber = gsub("GT-","",PLOT_CODE)) %>%
  dplyr::select(c(PlotNumber,X_NAD83,Y_NAD83))

gplots.in2 = st_read("data/GRTE_Vegmap_Plots/2003_Plots/Spatial_data/2003_plots_shpf/Plots_2003.shp") %>%
  # remove prefix from plot code
  mutate(PlotNumber = gsub("GT-","",PLOT_CODE)) %>%
  dplyr::select(c(PlotNumber,X_NAD83,Y_NAD83))

gplots.in = rbind(gplots.in1,gplots.in2)


gfuels = gplots.in %>%
  inner_join(gfuels.in, by="PlotNumber") %>%
  rename(Plot_code = PlotNumber) %>%
  dplyr::select(-c(Primary_surface_carriers,Comments))


####
# 2. convert coordinates and map plot locations
####

gfuels.sf = st_transform(gfuels, crs="+init=epsg:4326") %>%
  dplyr::select(-c(X_NAD83,Y_NAD83))

fuels.2019 = st_as_sf(x = fuels, coords = c("Postprocessed_easting", "Postprocessed_northing"),
                      crs = "+init=epsg:26912") %>%
  st_transform(crs = "+init=epsg:4326") %>%
  dplyr::select(c(Plot_code,FBFM,geometry))

fuels.sf = fuels.2019 %>%
  rbind(gfuels.sf)

# just 2019 plots
pal = colorNumeric(palette="Greens",domain=c(1:10))

leaflet(fuels.2019) %>%
  addTiles()  %>%
  addCircleMarkers(
    color = ~pal(as.numeric(FBFM)),
    # popup = mapLab,
    stroke = FALSE, fillOpacity = 0.8
  )


# all plots

leaflet(fuels.sf) %>%
  addTiles()  %>%
  addCircleMarkers(
    color = ~pal(as.numeric(FBFM)),
    # popup = mapLab,
    stroke = FALSE, fillOpacity = 0.8
  ) 


####
# 3. extract predictors from 30m rasters
####

### lidar + naip + veg type

grtel = stack("processed_data/GRTE_rasters/grte_30m_lidar_full.tif") 
# naip
grten1 = stack("processed_data/GRTE_rasters/grte_30m_naip_cv.tif")
grten2 = stack("processed_data/GRTE_rasters/grte_30m_naip_hom.tif")
grten3 = stack("processed_data/GRTE_rasters/grte_30m_naip_max.tif")
grten4 = stack("processed_data/GRTE_rasters/grte_30m_naip_mean.tif")
grten5 = stack("processed_data/GRTE_rasters/grte_30m_naip_min.tif")
grten6 = stack("processed_data/GRTE_rasters/grte_30m_naip_sd.tif")
# veg type, load unmasked version
grtev1 = raster("processed_data/GRTE_rasters/grte_30m_veg_con_nad2011.tif")
grtev2 = raster("processed_data/GRTE_rasters/grte_30m_veg_dec_nad2011.tif")
grtev3 = raster("processed_data/GRTE_rasters/grte_30m_veg_shr_nad2011.tif")

# create master veg type raster combining all three
grtev1.bin = reclassify(grtev1, c(-Inf,0.499,NA,0.499,1,1))
grtev2.bin = reclassify(grtev2, c(-Inf,0.499,NA,0.499,1,2))
grtev3.bin = reclassify(grtev3, c(-Inf,0.499,NA,0.499,1,3))
grtev = max(stack(grtev1.bin,grtev2.bin,grtev3.bin),na.rm=TRUE)


fuels.nad11 = fuels.sf %>%
  st_transform(crs = st_crs(grtel))

# extract coords for closest raster cell
fuels.pred = fuels.nad11 %>%
  cbind(extract(grtel, fuels.nad11, method="simple")) %>%
  cbind(extract(grten1, fuels.nad11, method="simple")) %>%
  cbind(extract(grten2, fuels.nad11, method="simple"))%>%
  cbind(extract(grten3, fuels.nad11, method="simple"))%>%
  cbind(extract(grten4, fuels.nad11, method="simple"))%>%
  cbind(extract(grten5, fuels.nad11, method="simple"))%>%
  cbind(extract(grten6, fuels.nad11, method="simple"))%>%
  cbind(extract(grtev, fuels.nad11, method="simple"))

head(fuels.pred)

### aggregate by veg type

fuels.veg = fuels.pred %>%
  rename(Veg_type = extract.grtev..fuels.nad11..method....simple..)

fuels.veg %>%
  filter(is.na(Veg_type)) # need to relabel Con_BerryGlade to conifer

fuels.fin = fuels.veg %>%
  mutate(Veg_type = ifelse(Plot_code == "Con_BerryGlade_Unburned3",1,Veg_type)) %>%
  filter(!is.na(Veg_type), !is.na(grte_30m_lidar_full.1))  # take out points with NA lidar values

names(fuels.fin) = c("Plot_code","FBFM",
                     "zmax","zmean","zcv","zp10","zp25","zp50","zp75","zp90","zp99","zfcc","zfsc",
                     "Rcv","Gcv","Bcv","NIRcv","NDVIcv",
                     "Rhom","Ghom","Bhom","NIRhom","NDVIhom",
                     "Rmax","Gmax","Bmax","NIRmax","NDVImax",
                     "Rmean","Gmean","Bmean","NIRmean","NDVImean",
                     "Rmin","Gmin","Bmin","NIRmin","NDVImin",
                     "Rsd","Gsd","Bsd","NIRsd","NDVIsd",
                     "Veg_type","geometry")

head(fuels.fin)


####
# 4. fbfm exploration
####

summary(as.factor(fuels.fin$FBFM))

fuels.fin %>%
  group_by(FBFM) %>%
  tally() %>%
  ggplot(aes(y=n*100/165,x=FBFM)) +
  geom_bar(stat="identity") +
  ylab("Frequency (%)") +
  theme_bw()

# some poorly represented FBFMs, will likely need to lump together
# based on similarities between fire behavior and structure in these FBFMs
# want to have a minimum # of representative plots, even if still unbalanced

# look at relationship between FBFM and predictor variables
fuels.fin %>%
  as.data.frame() %>%
  dplyr::select(-c(geometry,Veg_type)) %>%
  pivot_longer(cols=c("zmean":"NDVIsd"),names_to="predictors") %>%
  ggplot(aes(x=value, y=FBFM)) +
  facet_wrap(~predictors, scales="free_x") +
  geom_point()

fuels.fin %>%
  ggplot(aes(x=Veg_type, y=FBFM)) +
  geom_point()
# veg type likely to be important

# main issues with zp10 and NDVI cv (clear outliers) and with zfcc (unrealistic values > 100 for GS and other FBFMs)
# fix: remove zp10 and NDVIcv from analysis
# fix: set zfcc values > 1 or NA to 0

# finally, include TU4 with TU5 and SH1 with SH2, due to low # of observations
  
fuels.sub = fuels.fin %>%
  mutate(zfcc = ifelse(is.na(zfcc),0,
                       ifelse(zfcc>100,0,zfcc)),
         FBFM = ifelse(FBFM %in% c("SH1", "SH2"), "GS1",
                       ifelse(FBFM=="TU4","TU5",
                              ifelse(FBFM=="TL4","TL3",as.character(FBFM))))) %>%
  dplyr::select(-c(zp10,NDVIcv)) %>%
  as.data.frame() %>%
  # add index for later resampling
  mutate(index = seq.int(nrow(.)))

summary(fuels.sub)
summary(as.factor(fuels.fin$FBFM))
summary(as.factor(fuels.sub$FBFM))

fuels.sub %>%
  dplyr::select(-geometry) %>%
  write.csv("processed_data/fuels_map_variables/fbfms_6categories.csv",row.names=FALSE)

####
# 5. relationships among predictors
####

### correlation matrix, examine lidar and imagery separately for ease of viewing
start.l = which(colnames(fuels.sub)=="zmean")
end.l = which(colnames(fuels.sub)=="zfsc")
start.i = which(colnames(fuels.sub)=="Rcv")
end.i = which(colnames(fuels.sub)=="NDVIsd")

# Creating pairwise correlation matrix among predictors
corrMat.all = cor(fuels.sub[start.l:end.i])

# Color ramp
pal = brewer.pal(10, name = "RdYlBu")

## Making correlation plot
corrplot(corrMat.all, method = "circle", diag = F, cl.lim = c(-1, 1), 
         tl.col = "black", col = pal, bg = "grey90")
corrMat.all
# some correlations among lidar and naip predictors, but worst culprits within lidar and within naip

corrMat.lidar = cor(fuels.sub[start.l:end.l])
corrplot(corrMat.lidar, method = "circle", diag = F, cl.lim = c(-1, 1), 
         tl.col = "black", col = pal, bg = "grey90")
corrMat.lidar

## lots of strongly correlated predictors
# zmean, zp50, zp75, zp90, zp99,zfcc all strongly positively correlated
# zp25 strongly correlated with zp50 and zmean

corrMat.naip = cor(fuels.sub[start.i:end.i])
corrplot(corrMat.naip, method = "circle", diag = F, cl.lim = c(-1, 1), 
         tl.col = "black", col = pal, bg = "grey90")
corrMat.naip

## lots of correlations
# for visible bands, mins, means, maxes, etc. tend to be strongly correlated
# min, mean, cv, Rhom, NIRsd seem to tend to have lots of multicollinearity
# max, sd, Ghom, Bhom, Bcv, NIRhom less so

####
# 6. fit RF for 6 categories
####

# get column positions across all variables
start.c = which(colnames(fuels.sub)=="zmean")
end.c = which(colnames(fuels.sub)=="Veg_type")

# change FBFM and veg_type to factor
fuels.sub$FBFM = as.factor(fuels.sub$FBFM)
fuels.sub$Veg_type = as.factor(fuels.sub$Veg_type)

# use this to get random numbers
# set.seed ensures same result if code rerun, choose initial value based on random number between 1 and 1000
sample(1:1000,1)

### step one, create 10 data subsets, using upsampling to generate equal representation 
# in each set
# first generate 10 data partitions, each with 0.8 observations
set.seed(792)
part.samples = fuels.sub$FBFM %>%
  createDataPartition(p = 0.8, times = 10, list = TRUE)

# second use upsampling to add in additional observations from poorly represented classes
part.data=list()

samp.seeds = c(498, 951, 327,  85, 428, 345, 808, 477, 903, 480)

for(i in 1:length(part.samples)){
  
  print(i)
  part = fuels.sub[part.samples[[i]],]
  set.seed(samp.seeds[i])
  part.up = upSample(part, part$FBFM)
  
  # check for equal sampling
  print(summary(part.up$FBFM))
  
  part.data[[i]] = part.up$index
  
}

part.data

# subset dataset for only FBFM and predictors
fuels.mod = fuels.sub[c(2,start.c:end.c)]
ncol(fuels.mod) # 39 potential predictors

### rfe
# rfe controls
controls = rfeControl(method = "cv", 
                      number = 10,
                      index = part.data,
                      functions = rfFuncs)

set.seed(347)
# recursive feature elimination
caretRFE = rfe(y = fuels.mod$FBFM,  # response
               x = within(fuels.mod, rm(FBFM)),  # predictors
               rfeControl = controls,  # CV method and other parameters
               sizes = c(1:39))  # create models for all potential number of predictors

caretRFE # best model has 19 predictors (full set), 
plot(caretRFE) #  but starts leveling out earlier

str(caretRFE)

predictors(caretRFE)
# pred.new = predictors(caretRFE) # select top predictors
# or to minimize potential overfitting, using ~n/10
pred.new = predictors(caretRFE)[1:16] # ~n/10

### tune model
# pivot back to using original dataset
fuels.train = fuels.sub %>%
  dplyr::select(c(FBFM,c(pred.new)))

# calculate case weights from data, based on inverse of class proportion
cprobs = summary(fuels.train$FBFM)/length(fuels.train$FBFM) # proportion of observations by class
inv.weight = 1/cprobs  # inverse of proportion
rpart.weights = left_join(fuels.train, data.frame(FBFM = names(inv.weight), weight = inv.weight), by="FBFM") %>%
  dplyr::select(weight) %>%
  unlist() # list of inverse weights of same length and corresponding with observations
rpart.weights

ncol(fuels.train) 

# set control as cv, no upsampling because using case weights
control = trainControl(method="cv",number=10) 

# default mtry is either p/3 (breiman) or sqrt(p)
# recommend try default, default/2, default * 2, adding in multi possibilities here
tgrid = expand.grid(
  mtry=c(1:5,10),
  # mtry=c(1:6,12),
  splitrule = "gini",  
  min.node.size = c(1:5)
)

# model tuning
set.seed(491)
model_caret = train(FBFM~.,data=fuels.train,
                    method="ranger",
                    weights=rpart.weights,
                    trControl=control,
                    tuneGrid = tgrid,
                    importance="permutation",
                    num.trees=1000)

model_caret

## Getting optimal hyperparameters from cross-validation above
optMtry = as.numeric(model_caret$bestTune[1,1])
optSplit = as.character(model_caret$bestTune[1,2])
optMinNode = as.numeric(model_caret$bestTune[1,3])

### First assess predictor-response relationships with partial plots
model = ranger(formula =FBFM~., data = fuels.train,
               case.weights = rpart.weights,
               importance = "permutation",num.trees = 1000, 
               min.node.size = optMinNode, splitrule = optSplit,
               mtry = optMtry, probability = TRUE, seed = 460)

# partial values for each FBFM, also for groups of GS, TL, TU
# create function to apply to each FBFM or group
partial_plot = function(class.in) {
  partialVals = lapply(pred.new, FUN = function(x){
    autoplot(pdp::partial(model, pred.var = as.character(x),
                          type = "classification", which.class = class.in, prob = T))})
  
  sjPlot::plot_grid(partialVals, margin = c(0.1,0.1,0.1,0.1))
}

# each of these takes a while to run, use to assess whether relationships
# between predictors and probability of fuel model make sense
# assess each FBFM individually but also within FBFM groups (GS, TL, TU)
partial_plot("GS1")
partial_plot("GS2")
partial_plot(c("GS1","GS2"))

partial_plot("TL1")
partial_plot("TL3") 
partial_plot(c("TL1","TL3")) 

partial_plot("TU1")
partial_plot("TU5")

## Developing final RF model of disturbance severity with selected parameters
# and full dataset, case weights for sampling unbalanced data
model = ranger(formula =FBFM~., data = fuels.train,
               case.weights = rpart.weights,
               importance = "permutation",num.trees = 1000, 
               min.node.size = optMinNode, splitrule = optSplit,
               mtry = optMtry, seed = 460)

model # oob prediction error 


### assess model performance with repeated cross-validation
# repeated 5-fold CV, 10 repeats
# here have to create cvs by class to ensure 5 distinct cross-validation datasets

nfolds=5
nreps=10
ntot=50

# create seeds for each FBFM and repeat
cv.seedlist = list(c(675, 197, 974,  72, 588, 630),
                c(347, 173, 495, 253, 804, 306),
                c(360, 627, 629,  13, 728, 273),
                c(81, 876, 877, 964, 240, 558),
                c(514, 460, 679, 240, 476, 865),
                c(524, 875, 360, 576, 884, 772),
                c(224, 819, 585, 537, 669, 104),
                c(481, 144, 820, 405, 411,  48),
                c(508, 995, 980, 380, 546, 128),
                c(555,  80, 494, 811, 786, 302))

fuels.cvout = vector(mode="list",length=0)

for(j in 1:nreps) {
  
  # create empty list of nfolds x nreps
  fuels.cvrep = vector(mode="list",length=nfolds)
  
  cv.seeds = cv.seedlist[[j]]
  
  for(i in 1:length(levels(fuels.sub$FBFM))) {
    
    # break each FBFM class into 5 equal proportions
    levels(fuels.sub$FBFM)[i]
    fuels.class = fuels.sub[fuels.sub$FBFM == levels(fuels.sub$FBFM)[i],]
    
    set.seed(cv.seeds[i])
    fuels.cv = fuels.class$FBFM %>%
      createFolds(k=5,list=TRUE) 
    
    # then match these all up with their index from the full dataset
    for(k in 1:nfolds){
      
      fuels.cvrep[[k]] = c(fuels.cvrep[[k]],fuels.class[fuels.cv[[k]],]$index)
      
    }
  }

  fuels.cvout = c(fuels.cvout,fuels.cvrep)
  print(j)
}


fuels.cvout

# use model fit to rest of data to predict each cross-validation data set
cv.seeds2 = c(916, 308, 301, 255, 126,
              782, 193, 359, 747, 312,
              258, 716,  19, 629,  42,
              625, 892, 867, 731, 836,
              925, 557, 689, 404, 753,
              962, 846, 289, 138, 185,
              417, 369, 695, 844, 148,
              738, 422, 142, 141, 727,
              168, 719, 851, 965, 726,
              24, 732, 615, 250, 611)
obsVals = c()
predVals = c()

for(k in 1:ntot){
  ## Subsetting data
  testData = fuels.train[fuels.cvout[[k]],] ## Subsetting to test set
  trainData = fuels.train[-fuels.cvout[[k]],] ## Subsetting to training set
  
  print("test data summary:")
  print(summary(testData$FBFM))
  print(paste0("npred = ",ncol(fuels.train)-1))
  
  ## calculate case weights from data, based on inverse of class proportion
  trainProbs = summary(trainData$FBFM)/length(trainData$FBFM) # proportion of observations by class
  inv.tweight = 1/trainProbs  # inverse of proportion
  train.weights = left_join(trainData, data.frame(FBFM = names(inv.tweight), weight = inv.tweight), by="FBFM") %>%
    dplyr::select(weight) %>%
    unlist() # list of inverse weights of same length and corresponding with observations
  train.weights
  
  ## Fitting model to training data and predicting to test set
  mod = ranger(formula =FBFM~., data = trainData,
               case.weights = train.weights,
               importance = "permutation",num.trees = 1000, 
               min.node.size = optMinNode, splitrule = optSplit,
               mtry = optMtry, seed = cv.seeds2[k])
  
  ## And getting those values
  obsVals = c(obsVals, as.character(testData$FBFM))
  predVals = c(predVals, as.character(predict(mod, data = testData)$predictions))
}

# 10x repeated 5-fold CV accuracy
confusionMatrix(data= as.factor(predVals), reference = as.factor(obsVals))

# also assess balanced accuracy
cmat = table(obsVals, predVals)
cmat
sum(diag(cmat))/sum(cmat)  # accuracy, 67%
mean(diag(cmat/10)/summary(fuels.train$FBFM))  # balanced accuracy, 63%

### variable importance
vimp = stack(model$variable.importance)
vimp$ind = factor(vimp$ind, levels = vimp$ind[order(vimp$values, decreasing = F)])
vimp$values = vimp$values/sum(vimp$values)*100

# Then plotting variable importance
# note that correlated variables may "share" importance
ggplot(vimp, aes(x = ind, y = values, fill = ind)) + 
  geom_bar(stat = "identity", color = "black") + 
  coord_flip() + theme_bw() + xlab("Variable") + ylab("Relative Importance (%)") +
  scale_fill_manual(values = viridis(n = nrow(vimp), direction = -1)) +
  theme(legend.position = "none") 

### is variable importance stable?
# fit 10 models, compare
vimp.out = data.frame()
vimp.seeds = c(460, 386, 210,  58, 697, 865, 813, 347, 336, 149)

for(i in 1:length(vimp.seeds)) {
  
  model.v = ranger(formula =FBFM~., data = fuels.train,
                   case.weights = rpart.weights,
                   importance = "permutation",num.trees = 1000, 
                   min.node.size = optMinNode, splitrule = optSplit,
                   mtry = optMtry, seed = vimp.seeds[i])
  
  vimp.v = stack(model.v$variable.importance)
  vimp.v$ind = factor(vimp.v$ind, levels = vimp.v$ind[order(vimp.v$values, decreasing = F)])
  vimp.v$values = vimp.v$values/sum(vimp.v$values)*100
  
  vimp.out = rbind(vimp.out, cbind(data.frame(mod.n = i), vimp.v))
  
}

ggplot(vimp.out, aes(x = ind, y = values, fill = ind)) + 
  facet_wrap(~mod.n) +
  geom_bar(stat = "identity", color = "black") + 
  coord_flip() + theme_bw() + xlab("Variable") + ylab("Relative Importance (%)") +
  scale_fill_manual(values = viridis(n = nrow(vimp), direction = -1)) +
  theme(legend.position = "none") # general consistency among which variables are closer to top, closer to bottom

# average variable importance across these 10 models
vimp.avg = vimp.out %>%
  group_by(ind) %>%
  summarise(values=mean(values)) %>%
  mutate(ind = factor(as.character(ind), levels = ind[order(values, decreasing = F)]))

vimp.avg %>% 
  ggplot(aes(x = ind, y = values, fill = ind)) + 
  geom_bar(stat = "identity", color = "black") + 
  coord_flip() + theme_bw() + xlab("Variable") + ylab("Relative Importance (%)") +
  scale_fill_manual(values = viridis(n = nrow(vimp), direction = -1)) +
  theme(legend.position = "none")

# compare importance of LiDAR v. NAIP variables
vimp.group = vimp.avg %>%
  mutate(ind = as.character(ind)) %>%
  mutate(grouped_vars = ifelse(substr(ind,1,1)=="z","lidar",
                               ifelse(ind=="Veg_type","Veg_type","naip"))) %>%
  group_by(grouped_vars) %>%
  summarise(values=sum(values)) %>%
  mutate(grouped_vars = factor(grouped_vars, levels = grouped_vars[order(values, decreasing = F)])) 

vimp.group %>%
  ggplot(aes(x = grouped_vars, y = values, fill = grouped_vars)) + 
  geom_bar(stat = "identity", color = "black") + 
  coord_flip() + theme_bw() + xlab("Variable") + ylab("Relative Importance (%)") +
  scale_fill_manual(values = viridis(n = nrow(vimp.group), direction = -1)) +
  theme(legend.position = "none")


####
# 7. create FBFM fuels map
####

### create master fuels map prediction rasters

# read in masked and corrected lidar and naip predictors
lidar.map = stack("processed_data/GRTE_rasters/grte_30m_lidar_masked.tif")
naip.map = stack("processed_data/GRTE_rasters/grte_30m_naip_masked.tif")

names(lidar.map) = c("zmax","zmean","zcv","zp10","zp25","zp50","zp75","zp90","zp99","zfcc","zfsc")
names(naip.map) = c("Rmax","Gmax","Bmax","NIRmax","NDVImax",
                    "Rmin","Gmin","Bmin","NIRmin","NDVImin",
                    "Rmean","Gmean","Bmean","NIRmean","NDVImean",
                    "Rsd","Gsd","Bsd","NIRsd","NDVIsd",
                    "Rcv","Gcv","Bcv","NIRcv","NDVIcv",
                    "Rhom","Ghom","Bhom","NIRhom","NDVIhom")

# mask veg raster
grtev.mask = mask(crop(grtev,lidar.map[[1]]),lidar.map[[1]])
names(grtev.mask) = "Veg_type"

# full predictor set
full.map = stack(lidar.map,naip.map,grtev.mask)

summary(full.map[[42]])

### predict FBFMs

grte.fbfm = raster::predict(object=full.map,model=model, type="response",
                        fun = function(model, ...) predict(model, ...)$predictions,
                        na.rm=TRUE)

### reclassify based on landfire lookup table
grte.fbfm@data@attributes

# GS1: 121
# GS2: 122
# TL1: 181
# TL3: 183
# TU1: 161
# TU5: 165

grte.out = reclassify(grte.fbfm, c(0.5,1.5,121,
                                   1.5,2.5,122,
                                   2.5,3.5,181,
                                   3.5,4.5,183,
                                   4.5,5.5,161,
                                   5.5,6.5,165))

### write out
writeRaster(grte.out,"processed_data/GRTE_rasters/final_fuels_map/grte_30m_fbfms_nad2011.tif",overwrite=TRUE)

### compare distribution of FBFMs by vegetation type

# lookup table for veg types
veg.lookup = data.frame(Veg_type = factor(c(1:3)), Veg_name = c("Conifer","Deciduous","Shrubland"))
fbfm.lookup = data.frame(layer = c(121,122,161,165,181,183),
                         FBFM = c("GS1","GS2","TU1","TU5","TL1","TL3"))

# field data
# tot obs by veg type
field.tot = fuels.train %>%
  group_by(Veg_type) %>%
  summarise(n_veg = n())

# pct fbfm by veg type
field.veg = fuels.train %>%
  group_by(Veg_type,FBFM) %>%
  summarise(n = n()) %>%
  left_join(veg.lookup, by="Veg_type") %>%
  left_join(field.tot, by="Veg_type") %>%
  mutate(pct_veg = (n*100)/n_veg)

# mapped data
fbfm.stack = stack(grte.out, grtev.mask)

map.vals = values(fbfm.stack) %>%
  as.data.frame() %>%
  filter(!is.na(layer)) %>%
  mutate(Veg_type = factor(Veg_type))

# total for entire map and by veg type
tot = dim(map.vals)[1]
tot_veg = map.vals %>%
  group_by(Veg_type) %>%
  summarise(n_veg = n())

map.veg = map.vals %>%
  group_by(Veg_type,layer) %>%
  summarise(n = n(), pct_tot=(n*100)/tot) %>%
  left_join(veg.lookup, by="Veg_type") %>%
  left_join(fbfm.lookup, by="layer") %>%
  left_join(tot_veg, by="Veg_type") %>%
  mutate(pct_veg = (n*100)/n_veg)

# compare
field.veg
map.veg

# Color ramp
palf = brewer.pal(10, name = "Paired")[c(7:8,9:10,3:4)]

field.veg %>%
  ggplot(aes(x=Veg_name, y=pct_veg, fill=FBFM)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palf)

map.veg %>%
  ggplot(aes(x=Veg_name, y=pct_veg, fill=FBFM)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palf)
