#####
#
## Berry Fire fuels and fire severity analyses (Q1 & Q2)
#
#####

rm(list=ls())

# load libraries
library(dplyr) # version 0.8.3
library(tidyr) # version 1.0.2
library(ggplot2) # version 3.2.1
library(corrplot) # version 0.84
library(RColorBrewer) # version 1.1-2
library(car) # version 3.0-7
library(lme4) # version 1.1-21
library(lmerTest) # version 3.1-1
library(ape) # version 5.4-1
library(MuMIn) # version 1.43.15
library(MatchIt) # version 4.1.0
library(optimx) # version 2020-4.2


####
# 1. load data
####

# Q1: young v. mature forest
age.in = read.csv("processed_data/Berry_Fire_sample_polys/berry_fire_age_lidar.csv") %>%
  # remove these columns because repeated in naip predictors
  dplyr::select(-c(elev_m,slope_deg,aspect_deg)) %>%
  # join to naip predictors
  left_join(read.csv("processed_data/Berry_Fire_sample_polys/berry_fire_age_naip.csv"), 
            by=c("sample_id","age","rdnbr")) %>%
  dplyr::select(c(sample_id:rdnbr,elev_m:aspect_deg,zmax:zfsc,Rmax:NDVIhom))

# Q2: moderate v. extreme fire weather
wx.in = read.csv("processed_data/Berry_Fire_sample_polys/berry_fire_wx_lidar.csv") %>%
  # remove these columns because repeated in naip predictors
  dplyr::select(-c(elev_m,slope_deg,aspect_deg)) %>%
  # join to naip predictors
  left_join(read.csv("processed_data/Berry_Fire_sample_polys/berry_fire_wx_naip.csv"), 
            by=c("sample_id","weather","date","rdnbr")) %>%
  dplyr::select(c(sample_id:rdnbr,elev_m:aspect_deg,zmax:zfsc,Rmax:NDVIhom))

####
# 2. apply equations to predict fuel loads
####

# read final selected equations
model.fin = read.csv("analysis/fuels_prediction_map/model_final_selection.csv") %>%
  filter(veg_types %in% c("Conifer","Conifer_deciduous"))
model.pred = read.csv("analysis/fuels_prediction_map/model_final_predictors.csv")%>%
  filter(veg_types %in% c("Conifer","Conifer_deciduous"))

# load same function used to create fuels map
model_predict = function(metric.in, pred.in) {
  
  model.form = model.pred %>%
    filter(fuels_metric == metric.in) %>%
    dplyr::select(c(pred,coef))
  
  model.in = model.fin %>%
    filter(fuels_metric == metric.in)
  
  # set up equation for each possible number of predictors
  
  if (model.in$n_pred == 1) {
    
    print("1 predictor")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,")"))
    
    fuelsvar = coefi + (pred.in[[pred1]] * coef1) 
    
  } else if (model.in$n_pred == 2) {
    
    print("2 predictors")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    pred2 = as.character(model.form[3,]$pred)
    coef2 = round(as.numeric(model.form[3,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,") + (",coef2,pred2,")"))
    
    fuelsvar = coefi + (pred.in[[pred1]] * coef1) + (pred.in[[pred2]] * coef2) 
    
  } else if (model.in$n_pred == 3) {
    
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
    
    fuelsvar = coefi + (pred.in[[pred1]] * coef1) + (pred.in[[pred2]] * coef2) + (pred.in[[pred3]] * coef3)
    
  } else if (model.in$n_pred == 4) {
    
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
    
    fuelsvar = coefi + (pred.in[[pred1]] * coef1) + (pred.in[[pred2]] * coef2) + (pred.in[[pred3]] * coef3) + (pred.in[[pred4]] * coef4)
    
  } else if (model.in$n_pred == 5) {
    
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
    
    fuelsvar = coefi + (pred.in[[pred1]] * coef1) + (pred.in[[pred2]] * coef2) + (pred.in[[pred3]] * coef3) + (pred.in[[pred4]] * coef4) + (pred.in[[pred5]]*coef5)
    
  }

  return(fuelsvar)

  
}


# calculate fuel loads for age plots
age.fuels = age.in
# back-transform ln() for tree density
age.fuels$TPH_trees = exp(model_predict("TPH_trees",age.in))
# back-transform sqrt() for BA
age.fuels$BA_m2_ha = model_predict("BA_m2_ha",age.in)^2
age.fuels$CC_pct = model_predict("CC_pct",age.in)
age.fuels$CH_m = model_predict("SH_m",age.in)
# back-transform sqrt() for CFL
age.fuels$CFL_kg_m2 = model_predict("CFL_kg_m2",age.in)^2
# back-transform cube root for CBD
age.fuels$CBD_kg_m3 = model_predict("CBD_kg_m3",age.in)^3
# back-transform sqrt() for CBH
age.fuels$CBH_m = model_predict("CBH_m",age.in)^2
age.fuels$CWD_cover_pct = model_predict("CWD_cover_pct",age.in)
# back-transform ln() for CWD biomass
age.fuels$CWD_Mg_ha = exp(model_predict("CWD_Mg_ha",age.in))
age.fuels$Shrub_cover_pct = model_predict("Shrub_cover_pct",age.in)
age.fuels$Shrub_ht_m = model_predict("Shrub_ht_m",age.in)

summary(age.fuels)

# calculate fuel loads for wx plots
wx.fuels = wx.in
# back-transform ln() for tree density
wx.fuels$TPH_trees = exp(model_predict("TPH_trees",wx.in))
# back-transform sqrt() for BA
wx.fuels$BA_m2_ha = model_predict("BA_m2_ha",wx.in)^2
wx.fuels$CC_pct = model_predict("CC_pct",wx.in)
wx.fuels$CH_m = model_predict("SH_m",wx.in)
# back-transform sqrt() for CFL
wx.fuels$CFL_kg_m2 = model_predict("CFL_kg_m2",wx.in)^2
# back-transform cube root for CBD
wx.fuels$CBD_kg_m3 = model_predict("CBD_kg_m3",wx.in)^3
# back-transform sqrt() for CBH
wx.fuels$CBH_m = model_predict("CBH_m",wx.in)^2
wx.fuels$CWD_cover_pct = model_predict("CWD_cover_pct",wx.in)
# back-transform ln() for CWD biomass
wx.fuels$CWD_Mg_ha = exp(model_predict("CWD_Mg_ha",wx.in))
wx.fuels$Shrub_cover_pct = model_predict("Shrub_cover_pct",wx.in)
wx.fuels$Shrub_ht_m = model_predict("Shrub_ht_m",wx.in)

summary(wx.fuels)

### apply the same cut-offs as used in the fuels map when applicable
summary(age.fuels)
# TPH btwn 0-80000 good, BA 0+ good, CC 0-100 good, CH 0+ good, CFL 0-4 good,
# CBD needs fixed, <0 value, upper limit 0.4 good
# CBH 0-5 good, CWD cover 0-50 good
# CWD biomass needs fixed, lower limit 0 good, upper limit > 300 present
# Shrub_cover needs fixed, <0 value, upper limit 100 good
# Shrub_ht needs fixed, <0 value, upper limit 2 good

age.sub = age.fuels %>%
  mutate(CBD_kg_m3=ifelse(CBD_kg_m3<0,0,CBD_kg_m3),
         CWD_cover_pct=ifelse(CWD_cover_pct<0,0,CWD_cover_pct),
         CWD_Mg_ha = ifelse(CWD_Mg_ha>300,300,CWD_Mg_ha),
         Shrub_cover_pct=ifelse(Shrub_cover_pct<0,0,Shrub_cover_pct),
         Shrub_ht_m=ifelse(Shrub_ht_m<0,0,Shrub_ht_m)) %>%
  # remove lidar and naip predictors
  dplyr::select(-c(zmax:NDVIhom)) %>%
  # order factors logically
  mutate(age=factor(age,levels=c("young","old")))

summary(age.sub) # looks good

# wx
summary(wx.fuels)
# TPH between 0-80000 good, BA 0+ good
# CC needs fixed, <0 value, upper limit 100 good
# CH needs fixed, <0 value
# CFL 0-4 good
# CBD needs fixed, <0 value, upper limit 0.4 good
# CBH needs fixed, lower limit 0 good, upper limit >5 value
# CWD cover needs fixed, <0 value and >50 value
# CWD biomass needs fixed, lower limit good, upper limit > 300 value
# Shrub cover needs fixed, <0 value, upper limit 100 good
# Shrub ht needs fixed, <0 value, upper limit 2 good

wx.sub = wx.fuels %>%
  mutate(CC_pct=ifelse(CC_pct<0,0,CC_pct),
         CH_m=ifelse(CH_m<0,0,CH_m),
         CBD_kg_m3=ifelse(CBD_kg_m3<0,0,CBD_kg_m3),
         CBH_m = ifelse(CBH_m>5,5,CBH_m),
         CWD_cover_pct=ifelse(CWD_cover_pct<0,0,CWD_cover_pct),
         CWD_cover_pct = ifelse(CWD_cover_pct>50,50,CWD_cover_pct),
         CWD_Mg_ha = ifelse(CWD_Mg_ha>300,300,CWD_Mg_ha),
         Shrub_cover_pct=ifelse(Shrub_cover_pct<0,0,Shrub_cover_pct),
         Shrub_ht_m=ifelse(Shrub_ht_m<0,0,Shrub_ht_m)) %>%
  # remove lidar and naip predictors
  dplyr::select(-c(zmax:NDVIhom)) %>%
  # make burn date a factor
  mutate(date = factor(date))

summary(wx.sub)

# will match up age comparison in next step
# write out wx data here
write.csv(wx.sub, "analysis/q2_fuels_fire_severity_models/weather_fuels_severity.csv",row.names=FALSE)

####
# 3. q1: fuels and fire severity in young v. mature forest
####

### use propensity scores based on slope and elevation to match up age class samples
age.sub$age_class = as.numeric(age.sub$age)-1
match = matchit(age_class~elev_m+slope_deg, method="nearest", data=age.sub, distance="glm", link="probit",discard="both")
matched_data =match.data(match)

# initial plot comparison
age.plots = matched_data %>%
  dplyr::select(-c(x_ctr,y_ctr)) %>%
  pivot_longer(cols=c(TPH_trees:Shrub_ht_m,rdnbr:aspect_deg))

age.plots %>%
  mutate(age = factor(age, levels=c("young","old"))) %>%
  ggplot(aes(y=value,x=age)) +
  facet_wrap(~name,scales="free") +
  geom_boxplot() +
  theme_bw() 

### wilcoxon rank sum test for differences in distribution
# essentially whether median values differ
young.sub = matched_data %>% filter(age=="young")
old.sub = matched_data %>% filter(age=="old")

# default values are two-sided
wilcox.test(young.sub$TPH_trees, old.sub$TPH_trees) # p < 0.001, p = 0.000000014
wilcox.test(young.sub$BA_m2_ha, old.sub$BA_m2_ha) # ns, p = 0.22
wilcox.test(young.sub$CC_pct, old.sub$CC_pct) # p < 0.01, p = 0.0020
wilcox.test(young.sub$CH_m, old.sub$CH_m) # p < 0.001, p = 0.0000000000031
wilcox.test(young.sub$CFL_kg_m2, old.sub$CFL_kg_m2) # ns, p = 0.64
wilcox.test(young.sub$CBD_kg_m3, old.sub$CBD_kg_m3) # p < 0.001, p = 0.0000066
wilcox.test(young.sub$CBH_m, old.sub$CBH_m) # ns, p = 0.60
wilcox.test(young.sub$CWD_cover_pct, old.sub$CWD_cover_pct) # p < 0.001, p = 0.000072
wilcox.test(young.sub$CWD_Mg_ha, old.sub$CWD_Mg_ha) # ns, p = 0.72
wilcox.test(young.sub$Shrub_cover_pct, old.sub$Shrub_cover_pct) # p < 0.05, p = 0.0028
wilcox.test(young.sub$Shrub_ht_m, old.sub$Shrub_ht_m) # ns, p = 0.38
wilcox.test(young.sub$rdnbr, old.sub$rdnbr) # ns, p = 0.34
wilcox.test(young.sub$elev_m, old.sub$elev_m) # ns, p = 0.55
wilcox.test(young.sub$slope_deg, old.sub$slope_deg) # ns, p = 0.10

# examine median rdnbr
# stand replacing severity is RdNBR = 675 per Harvey et al. 2016
summary(young.sub$rdnbr) # median > 675, ~709
summary(old.sub$rdnbr) # median > 675, ~758

summary(young.sub)
summary(old.sub)

# write out data
matched_data %>%
  # remove unneeded columns
  dplyr::select(-c(age_class:subclass)) %>%
  write.csv("analysis/q1_young_old_forest_comparison/young_old_fuels_severity.csv",row.names=FALSE)

# write out median values
matched_data %>%
  dplyr::select(-subclass) %>%
  group_by(age) %>%
  summarise_all(median) %>%
  write.csv("analysis/q1_young_old_forest_comparison/young_old_forest_medians.csv",row.names=FALSE)


####
# 4. q2: predicting fire severity under moderate v. extreme fire weather
####

### first compare RdNBR
wxHi = wx.sub %>%
  filter(weather=="extreme") %>%
  mutate(date=as.factor(as.character(date)))
summary(wxHi)
# mean rdnbr is 705, median 717 so majority stand-replacing
# 43 cells from mature forest, 7 reburns

wxMod = wx.sub %>%
  filter(weather=="moderate") %>%
  mutate(date=as.factor(as.character(date)))
summary(wxMod)
# mean rdnbr is 493, median 440, majority not stand-replacing
# 42 cells from mature forest, 8 reburns

hist(wx.sub$rdnbr)
hist(wxMod$rdnbr)
hist(wxHi$rdnbr)
wilcox.test(wxHi$rdnbr,wxMod$rdnbr) # p < 0.05, p = 0.013

### use transformations as applicable for predictor variables
# if transformations applied during fuels regression fits improve normality,
# carry them over here

hist(wx.sub$TPH_trees)
hist(log(wx.sub$TPH_trees)) # ok
hist(wx.sub$BA_m2_ha)
hist(sqrt(wx.sub$BA_m2_ha)) # ok
hist(wx.sub$CFL_kg_m2)
hist(sqrt(wx.sub$CFL_kg_m2)) # ok
hist(wx.sub$CBD_kg_m3)
hist((wx.sub$CBD_kg_m3)^(1/3)) # ok
hist(wx.sub$CBH_m)
hist(sqrt(wx.sub$CBH_m)) # ok
hist(wx.sub$CWD_Mg_ha)
hist(log(wx.sub$CWD_Mg_ha)) # ok
hist(wx.sub$CWD_cover_pct)

### start with correlations, remove highly correlated (r > 0.7) variables
# subset and apply transformations
wx.corr = wx.sub %>%
  mutate(TPH_trees = log(TPH_trees),
         BA_m2_ha = sqrt(BA_m2_ha),
         CFL_kg_m2 = sqrt(CFL_kg_m2),
         CBD_kg_m3 = CBD_kg_m3^(1/3),
         CBH_m = sqrt(CBH_m),
         CWD_Mg_ha = log(CWD_Mg_ha),
         CWD_cover_pct = sqrt(CWD_cover_pct)) %>%
  dplyr::select(c(rdnbr,TPH_trees:Shrub_ht_m))

# Creating pairwise correlation matrix among predictors
corrMat = wx.corr %>%
  cor()

# Color ramp
pal = brewer.pal(10, name = "RdYlBu")

## Making correlation plot
corrplot(corrMat, method = "circle", diag = F, cl.lim = c(-1, 1), 
         tl.col = "black", col = pal, bg = "grey90")
abs(corrMat)>0.7

## highly correlated variables
# BA x CC, CH, CBH
# CC x CBD
# CBD x CBH
# CC x CBH, Shrub ht in moderate severity
# CH x CBH, CWD cover in moderate severity

# of these, CBH expected to influence fire severity over BA
# CBD, CBH good rationale for expecting them to influence fire severity; retain CBH based on relationship with rdnbr
# CBH expected to influence fire severity over CC, CH
# look at predictors x response
wx.corr %>%
  pivot_longer(cols=c(TPH_trees:Shrub_ht_m),names_to="fuel") %>%
  ggplot(aes(x=(value),y=(rdnbr)))+
  facet_wrap(~fuel, scales="free") +
  geom_point() +
  theme_bw()
# looks like a lot of mess!

# how closely correlated with rdnbr
abs(corrMat)[1,] # CBH is higher
# drop BA, CBD, CC, CH

corrMat = wx.corr %>%
  dplyr::select(-c(BA_m2_ha,CBD_kg_m3,CC_pct, CH_m)) %>%
  cor()
abs(corrMat)>0.7  # good

# now examine extreme and moderate weather separately to ensure no
# additional correlated predictors
wxHi.corr = wxHi %>%
  mutate(TPH_trees = log(TPH_trees),
         BA_m2_ha = sqrt(BA_m2_ha),
         CFL_kg_m2 = sqrt(CFL_kg_m2),
         CBD_kg_m3 = CBD_kg_m3^(1/3),
         CBH_m = sqrt(CBH_m),
         CWD_Mg_ha = log(CWD_Mg_ha)) %>%
  dplyr::select(c(rdnbr,TPH_trees:Shrub_ht_m))

corrMat = wxHi.corr %>%
  dplyr::select(-c(BA_m2_ha,CBD_kg_m3,CC_pct, CH_m)) %>%
  cor()
abs(corrMat)>0.7  # good

# moderate wx
wxMod.corr = wxMod %>%
  mutate(TPH_trees = log(TPH_trees),
         BA_m2_ha = sqrt(BA_m2_ha),
         CFL_kg_m2 = sqrt(CFL_kg_m2),
         CBD_kg_m3 = CBD_kg_m3^(1/3),
         CBH_m = sqrt(CBH_m),
         CWD_Mg_ha = log(CWD_Mg_ha)) %>%
  dplyr::select(c(rdnbr,TPH_trees:Shrub_ht_m))

corrMat = wxMod.corr %>%
  dplyr::select(-c(BA_m2_ha,CBD_kg_m3,CC_pct, CH_m)) %>%
  cor()
abs(corrMat)>0.7  # good

### fit models
### extreme weather days

# mixed model
model.x=lmer(rdnbr~log(TPH_trees)+sqrt(CFL_kg_m2)+sqrt(CBH_m)+CWD_cover_pct+log(CWD_Mg_ha)+
               Shrub_cover_pct+Shrub_ht_m+(1|date), data=wxHi, REML=FALSE,
             # specify optimizer for model convergence
             control = lmerControl(
               optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

# assess assumptions
plot(model.x) # seems ok, a couple extreme values, some unequal variance, no clear non-linear trend
qqPlot(residuals(model.x)) # looks a little funky at upper tail, within normal limits
summary(model.x)

# spatial autocorrelation
plot(residuals(model.x)~wxHi$x_ctr) # look good
plot(residuals(model.x)~wxHi$y_ctr) # look good
sp.mat = as.matrix(dist(cbind(wxHi$x_ctr,wxHi$y_ctr)))
sp.inv = 1/sp.mat
diag(sp.inv) = 0
# moran's i test for each variable
Moran.I(residuals(model.x),sp.inv) # nope

# model selection
options(na.action = "na.fail")
x.sel  = dredge(model.x, beta="none", rank="AICc") 

head(x.sel)
dim(x.sel[x.sel$delta<2])
x.sel[x.sel$delta<2] # 2 top models

# ID top models within AICc <2 and all predictors significant at p<0.05
summary(get.models(x.sel,1)[[1]]) # yes
summary(get.models(x.sel,2)[[1]]) # no

topmod.x = lmer(rdnbr~log(TPH_trees)+Shrub_cover_pct+(1|date), data=wxHi, REML=FALSE,
                control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

# check assumptions
plot(topmod.x, which=1) # variance a little unequal, but ok
qqPlot(residuals(topmod.x)) # a little skew, but ok

# test significance of random effect
redmod.x = lm((rdnbr)~log(TPH_trees)+Shrub_cover_pct,data=wxHi)
anova(topmod.x,redmod.x)

# model results
summary(topmod.x)
r.squaredGLMM(topmod.x)  # r2 fixed 0.27, total 0.27

summary(redmod.x) # r2 = 0.27

### moderate weather days

# mixed model
model.m=lmer(rdnbr~log(TPH_trees)+sqrt(CFL_kg_m2)+sqrt(CBH_m)+CWD_cover_pct+log(CWD_Mg_ha)+
               Shrub_cover_pct+Shrub_ht_m+(1|date), data=wxMod, REML=FALSE,
             control = lmerControl(
               optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

# assess assumptions
plot(model.m) # seems ok, some unequal variance, no clear non-linear trend
qqPlot(residuals(model.m)) # ok
summary(model.m)

# spatial autocorrelation
plot(residuals(model.m)~wxMod$x_ctr) # look good
plot(residuals(model.m)~wxMod$y_ctr) # look good
sp.mat = as.matrix(dist(cbind(wxMod$x_ctr,wxMod$y_ctr)))
sp.inv = 1/sp.mat
diag(sp.inv) = 0
# moran's i test for each variable
Moran.I(residuals(model.m),sp.inv) # nope

# model selection
options(na.action = "na.fail")
m.sel  = dredge(model.m, beta="none", rank="AICc") 

head(m.sel)
m.sel[m.sel$delta<2] # 
dim(m.sel[m.sel$delta<2]) # 5 top models

# ID top models within AICc <2 and all predictors significant at p<0.05
summary(get.models(m.sel,1)[[1]]) # no
summary(get.models(m.sel,2)[[1]]) # yes
summary(get.models(m.sel,3)[[1]]) # no
summary(get.models(m.sel,4)[[1]]) # yes
summary(get.models(m.sel,5)[[1]]) # no

# check assumptions and fit of top models
# model 1
topmod.m = lmer(rdnbr~log(CWD_Mg_ha)+log(TPH_trees)+Shrub_cover_pct+(1|date), REML=FALSE, data=wxMod,
                control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

plot(topmod.m, which=1) # a little unequal variance
qqPlot(residuals(topmod.m)) # ok

# test significance of random effect
redmod.m = lm(rdnbr~log(CWD_Mg_ha)+log(TPH_trees)+Shrub_cover_pct,data=wxMod)
anova(topmod.m,redmod.m)

# model results
summary(topmod.m)
r.squaredGLMM(topmod.m)  # r2 fixed 0.15, r2 full 0.51

summary(redmod.m)  # r2 0.16

# model 2
topmod.m2 = lmer(rdnbr~log(TPH_trees)+(1|date), REML=FALSE, data=wxMod)
plot(topmod.m2, which=1) # a little unequal variance
qqPlot(residuals(topmod.m2)) #

# test significance of random effect
redmod.m2 = lm(rdnbr~log(TPH_trees),data=wxMod)

anova(topmod.m2,redmod.m2)

# model results
summary(topmod.m2)
r.squaredGLMM(topmod.m2)  # r2 fixed 0.06, r2 full 0.46

summary(redmod.m2)  # r2 0.05
