#####
#
## fit las metrics to fuels characteristics
#
#####

rm(list=ls())

# load libraries
library(tidyr) # version 1.0.2
library(dplyr) # version 0.8.3
library(tibble) # version 2.1.3
library(leaps) # version 3.0
library(qpcR) # version  1.4-1
library(caret) # version 6.0-84
library(car) # version 3.0-7
library(asbio) # version 1.6-3

####
# 1. load data
####

fuels.in = read.csv("processed_data/fuels_map_variables/field_plot_2019_fuels.csv") 

ldr.in = read.csv("processed_data/fuels_map_variables/lidar_metrics.csv")
naip.in = read.csv("processed_data/fuels_map_variables/naip_metrics.csv")

# join dfs
fuels = fuels.in %>%
  left_join(ldr.in, by="Plot_code") %>%
  left_join(naip.in, by="Plot_code")

fuels.f = fuels %>%
  filter(Vegetation_type %in% c("Conifer","Deciduous"))

fuels.c = fuels %>%
  filter(Vegetation_type=="Conifer")

fuels.s = fuels %>%
  filter(Vegetation_type=="Shrub")

head(fuels.c)
summary(fuels.c)

####
# 2. functions
####

### write formulas for initial exploration, full model fit
# ind vars
# have to break into different groups to avoid overfitting
names(fuels)
names(fuels)[22:31] # lidar metrics, omitting zmax
names(fuels)[c(32:37,50:61)]  # R, NIR, NDVI
names(fuels)[38:49]  # B & G

# function to write lm formula
fuels.lm1 = function(dep.var) {
  as.formula(paste0(dep.var, "~", paste(names(fuels.fit)[22:31], collapse="+")))
}

# function to write lm formula
fuels.lm2 = function(dep.var) {
  as.formula(paste0(dep.var, "~", paste(names(fuels.fit)[c(32:37,50:61)], collapse="+")))
}

# function to write lm formula
fuels.lm3 = function(dep.var) {
  as.formula(paste0(dep.var, "~", paste(names(fuels.fit)[38:49], collapse="+")))
}

# function to write full lm formula
fuels.lm4 = function(dep.var) {
  as.formula(paste0(dep.var, "~", paste(names(fuels.fit)[c(22:61)], collapse="+")))
}

### helper functions for VIF, LOOCV, partial R2

# function to extract max VIF value and correlation
get_max_vif = function(dep.var, subset.name) {
  
  vif.out = data.frame()
  
  for(i in 1:length(summary(subset.name)$bic)) {
    
    # skip models with fewer than 2 predictors in addition to intercept
    if(length(coef(subset.name,i))<3) {
      vif.sub = data.frame(x = i,vif = 0, corr=0)
      vif.out = rbind(vif.out,vif.sub)
    } else {
      # fit each model
      red.form = as.formula(paste0(dep.var,"~",paste(names(coef(subset.name,i))[-1], collapse="+")))
      red.model = (lm(red.form, data=fuels.fit))
      corr.mat = cor((fuels.fit)[c(paste(names(coef(subset.name,i))[-1]))])
      # get maximum vif and maximum correlation
      vif.sub = data.frame(x = i, vif = max(vif(red.model)),
                           corr = max(abs(corr.mat[corr.mat<1])))
      vif.out = rbind(vif.out,vif.sub)
    }
  }
  return(vif.out)
}

# training for leave one out cross-validation
train.control = trainControl(method = "LOOCV")

# calc partial R2
partial_R2 = function(dep.variable,coeff.list) {
  
  # fit model
  mod.form = as.formula(paste0(dep.variable,"~",paste(coeff.list, collapse="+")))
  
  full.model = lm(mod.form, data=fuels.fit)
  
  partial.out = data.frame()
  
  for(i in 1:length(coeff.list)) {
    coeff=coeff.list[i]
    red.list = coeff.list[!coeff.list %in% coeff.list[i]]
    # fit reduced model without coeff of interest
    rmod.form = as.formula(paste0(dep.variable,"~",paste(red.list, collapse="+")))
    red.model = lm(rmod.form, data=fuels.fit)
    
    # calc partial R2 based on comparison of reduced and full model
    partr = partial.R2(red.model,full.model)
    
    partial.temp = data.frame(coeff=coeff,partial_r2=partr)
    partial.out = rbind(partial.out,partial.temp)
  }
  
  return(partial.out)
  
}

### functions for summarizing model fit and assessment stats
# create model summary data
# requires model.metric, model.veg, model.n, model.transform, var.in, model.df, model.subsets, fuels.fit to exist
# set up pipeline for rmse transform based on log or power transformation
# power.in allows specifying transformation of predictions to calculate rmse_transform
# allow specification of VIF, max corr cutoffs
summarize_models = function(power.in=1,
                            vif.cut = 5, corr.cut = 0.7) {
  
  # print input information as double check
  print(model.transform)
  print(var.in)
  print(power.in)

  # first create model dataframe with info on BIC, adjusted R2, VIF, max correlation
  # subset this by VIF and max correlation to get model short list
  model.df = data.frame(x = seq(1:length(summary(model.subsets)$bic))) %>%
    mutate(bic = summary(model.subsets)$bic,
           adj_rsq = summary(model.subsets)$adjr2) %>%
    # also extract max VIF and correlation
    left_join(get_max_vif(var.in,model.subsets), by="x") %>%
    filter(vif < vif.cut, corr < corr.cut)
  
  # create output dataframe
  model.out = model.df %>%
    add_column(fuels_metric = model.metric,
               veg_types = model.veg, 
               n = model.n, 
               transform=model.transform,
               n_pred = as.integer(NA),
               .before=1) %>%
    rename(model_number=x,
           max_correlation=corr) %>%
    add_column(rsq = as.numeric(NA),
               rmse = as.numeric(NA),
               rmse_transform = as.numeric(NA),
               loocv_rmse = as.numeric(NA),
               loocv_rsq = as.numeric(NA),
               rmse_diff = as.numeric(NA))
  
  # next iterate through each model in shortlist to determine parsimony, RMSE, LOOCV stats
  for(i in 1:dim(model.df)[1]) {
    
    # select 1 model at a time
    model.nbr = model.df$x[i] 
    
    # fit this model using var.in for any transformations
    model.formula = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.nbr))[-1], collapse="+")))
    model.reduced = lm(model.formula, data=fuels.fit)
    
    # number of predictors
    model.out$n_pred[i] = length(coef(model.subsets,model.nbr)[-1])
    
    # rsq
    model.out$rsq[i] = summary(model.reduced)$r.squared
    
    # RMSE
    model.out$rmse[i] = qpcR::RMSE(model.reduced)  
    
    # manually calculate for transformed
    # options for power or for log
    if(model.transform=="log") {
      model.out$rmse_transform[i] = sqrt((sum((exp(predict(model.reduced)) - fuels.fit$fuels_metric)^2))/model.n)
    } else if(model.transform=="log_plus1") {
      model.out$rmse_transform[i] = sqrt((sum(((exp(predict(model.reduced))-1) - fuels.fit$fuels_metric)^2))/model.n)
    } else {
      model.out$rmse_transform[i] = sqrt((sum(((predict(model.reduced)^(1/power.in)) - fuels.fit$fuels_metric)^2))/model.n)
      
    }
    
    # leave one out cross validation
    model.loocv = train(model.formula, data=fuels.fit, method="lm",
                        trControl=train.control)  
    
    model.out$loocv_rmse[i] = model.loocv$results$RMSE
    model.out$loocv_rsq[i] = model.loocv$results$Rsquared

    # RMSE difference
    model.out$rmse_diff[i] = model.out$loocv_rmse[i] - model.out$rmse[i]
    
  }
  
  return(model.out)
  
}
    
# create predictor summary data
summarize_predictors = function() {
  
  # calculate partial r2 for each predictor
  model.partr2 = partial_R2(var.in,names(coef(model.subsets,model.sel))[-1]) %>%
    mutate(coeff = as.character(coeff))
  
  # create output dataframe
  pred.out = data.frame(fuels_metric = model.metric, veg_types = model.veg, n = model.n,
                        n_pred = length(coef(model.subsets,model.sel)[-1]),
                        transform=model.transform,
                        model_number=model.sel,
                        pred = names(coef(model.subsets,model.sel)),
                        coef = coef(model.subsets,model.sel)) %>%
    mutate(pred = as.character(pred)) %>%
    left_join(model.partr2, by = c("pred" = "coeff"))
  
  return(pred.out)
}


####
# 3. Final fits for conifer and deciduous forest combined
####

# master output dataframes
summary.out = data.frame()
predictors.out = data.frame()

# change fuels.fit to df of interest
fuels.fit = fuels.f
str(fuels.fit)
summary(fuels.fit)
head(fuels.fit)

# predictors
predictors = names(fuels.fit)[c(22:61)]

model.veg = "Conifer_deciduous"
model.n = dim(fuels.fit)[1]

### a. stand density

model.metric = "TPH_trees"
model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("log(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("log(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("log(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("log(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("log(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("log(fuels_metric)"), data=fuels.fit))  # 

# transformation chosen based on boxcox
model.transform = "log"
var.in = "log(fuels_metric)"

# create formula
model.full = as.formula(paste0(var.in, "~", paste(predictors, collapse="+"))) 

# exhaustive best subsets selection
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power=(1),
                              vif.cut = 5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions

model.sel=45
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
partial_R2("log(fuels_metric)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=log(fuels_metric))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=log(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name
fuels.fit = fuels.f


### b. stand basal area

model.metric = "BA_m2_ha"
model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("sqrt(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("sqrt(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("sqrt(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("sqrt(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("sqrt(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("sqrt(fuels_metric)"), data=fuels.fit))  # 

# boxcox suggests sqrt transform
model.transform = "sqrt"
var.in = "sqrt(fuels_metric)"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("NIRmean","zp99","zp90","zfsc")] # remove
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 20, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation, adj r2 cutoffs
model.summ = summarize_models(power.in=(1/2),
                              vif.cut=5, corr.cut=0.7)  

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=41
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("sqrt(fuels_metric)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=sqrt(fuels_metric))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=sqrt(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name
fuels.fit = fuels.f


### c. canopy cover

model.metric = "CC_pct"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("(fuels_metric)"), data=fuels.fit))  # 

model.transform = "none"
var.in = "fuels_metric"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zfsc")] 
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 30, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 
# first model did not meet linearity assumption

model.sel=122
model.notes = "bestLOOCVRmse_meetsassumptions"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name
fuels.fit = fuels.f


### d. stand height

model.metric = "SH_m"
model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("(fuels_metric)"), data=fuels.fit))  # 

model.transform = "none"
var.in = "fuels_metric"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("NDVIcv")] 
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+"))) 

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions

model.sel= 47
model.notes="bestLOOCVrmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
partial_R2("(fuels_metric)",names(coef(model.subsets,model.sel))[-1])


fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=(fuels_metric))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))


# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name
fuels.fit = fuels.f


####
# 4. Final fits for conifer forest only
####

# change fuels.fit to df of interest
fuels.fit = fuels.c
str(fuels.fit)
summary(fuels.fit)
head(fuels.fit)

predictors = names(fuels.fit)[c(22:61)]

model.veg = "Conifer"
model.n = dim(fuels.fit)[1]

### a. canopy fuel load

model.metric = "CFL_kg_m2"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric) 

head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("(fuels_metric)^(1/2)"), data=fuels.fit))  # 

model.transform = "sqrt"
var.in = "fuels_metric^(1/2)"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zp90","zmean","zp50","zp75")]
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting (consult with MGT)
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=(1/2),
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 
# first 3 models do not meet linearity assumption

model.sel=31
model.notes="bestLOOCVrmse_meetsassumptions"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric^(1/2)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=(fuels_metric)^(1/2))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric)^(1/2))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original 
fuels.fit = fuels.c


### b. canopy bulk density

model.metric = "CBD_kg_m3"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric) 

head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric)^(1/3)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("(fuels_metric)^(1/3)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric)^(1/3)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("(fuels_metric)^(1/3)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric)^(1/3)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("(fuels_metric)^(1/3)"), data=fuels.fit))  # 

model.transform = "cube_root"
var.in = "fuels_metric^(1/3)"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zp25")]
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 20, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=(1/3),
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 
# first model does not meet linearity assumption

model.sel=61
model.notes="bestLOOCVRmse_meetsassumptions"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
# cor(fuels.fit[c(paste(predictors))])
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric^(1/3)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=(fuels_metric)^(1/3))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric)^(1/3))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))


# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original 
fuels.fit = fuels.c


### c. canopy base height

model.metric = "CBH_m"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric) 

head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric)^(1/2)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("(fuels_metric)^(1/2)"), data=fuels.fit))  # 

model.transform = "sqrt"
var.in = "fuels_metric^(1/2)"

# create formula
model.full = fuels.lm4(var.in)

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=(1/2),
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=42
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
# cor(fuels.fit[c(paste(predictors))])
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("sqrt(fuels_metric)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=sqrt(fuels_metric))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=sqrt(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original 
fuels.fit = fuels.c


### d. cwd cover

model.metric = "CWD_cover_pct"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("fuels_metric"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  # 

model.transform = "none"
var.in = "fuels_metric"

# create formula
model.full = fuels.lm4(var.in)

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation, adj r2 cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=41
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.c


### e. cwd biomass

model.metric = "CWD_Mg_ha"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("log(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("log(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("log(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("log(fuels_metric)"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("log(fuels_metric)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("log(fuels_metric)"), data=fuels.fit))  # 

model.transform = "log"
var.in = "log(fuels_metric)"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zfsc")]
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 30, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=122
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
partial_R2("log(fuels_metric)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=log(fuels_metric))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=log(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.c


### f. shrub cover 

model.metric = "Shrub_cover_pct"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  
plot(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  
boxcox(lm(fuels.lm2("fuels_metric"), data=fuels.fit)) 
plot(lm(fuels.lm2("fuels_metric"), data=fuels.fit)) 
boxcox(lm(fuels.lm3("fuels_metric"), data=fuels.fit)) 
plot(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  

model.transform = "none"
var.in = "fuels_metric"

# create formula
model.full = fuels.lm4(var.in)  

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation, adj r2 cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 
# first 2 models not meet assumption of linearity

model.sel=45
model.notes="bestLOOCVRmse_meetsassumptions"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.c


### g. shrub height 

model.metric = "Shrub_ht_m"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  
plot(lm(fuels.lm1("fuels_metric"), data=fuels.fit)) 
boxcox(lm(fuels.lm2("fuels_metric"), data=fuels.fit))  
plot(lm(fuels.lm2("fuels_metric"), data=fuels.fit)) 
boxcox(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  
plot(lm(fuels.lm3("fuels_metric"), data=fuels.fit)) 

model.transform = "none"
var.in = "fuels_metric"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zp10","Bcv")] 
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 5 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 20, nvmax = 5, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=88
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.c


####
# 5. Final fits for shrubland
####

# restrict to 3 potential predictors

# change fuels.fit to df of interest
fuels.fit = fuels.s
str(fuels.fit)
summary(fuels.fit)
head(fuels.fit)

predictors = names(fuels.fit)[c(22:61)]

model.veg = "Shrubland"
model.n = dim(fuels.fit)[1]


### a. shrub cover 

model.metric = "Shrub_cover_pct"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  # 
plot(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("fuels_metric"), data=fuels.fit))  # 
plot(lm(fuels.lm2("fuels_metric"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  # 
plot(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  # 

model.transform = "none"
var.in = "fuels_metric"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zp10","Bcv","Rsd","zp99","Gmean")] 
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 3 predictors to minimize overfitting
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 40, nvmax = 3, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

model.summ = summarize_models(power.in=(1),
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=87
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.s


### b. shrub height 

model.metric = "Shrub_ht_m"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("fuels_metric"), data=fuels.fit)) 
plot(lm(fuels.lm1("fuels_metric"), data=fuels.fit)) 
boxcox(lm(fuels.lm2("fuels_metric"), data=fuels.fit))  
plot(lm(fuels.lm2("fuels_metric"), data=fuels.fit)) 
boxcox(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  
plot(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  

model.transform = "none"
var.in = "fuels_metric"

# create formula
# remove predictors causing issues with collinearity and model linearity
pred.new = predictors[!predictors %in% c("zcv")] 
model.full = as.formula(paste0(var.in, "~", paste(pred.new, collapse="+")))

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 3 predictors to minimize overfitting (consult with MGT)
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 20, nvmax = 3, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=56
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.s


### c. shrub percent dead 

model.metric = "Shrub_pctDead"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric+1)^(1/3)"), data=fuels.fit))  
plot(lm(fuels.lm1("(fuels_metric)^(1/3)"), data=fuels.fit)) 
boxcox(lm(fuels.lm2("(fuels_metric+1)^(1/3)"), data=fuels.fit))  
plot(lm(fuels.lm2("(fuels_metric)^(1/3)"), data=fuels.fit))  
boxcox(lm(fuels.lm3("(fuels_metric+1)^(1/3)"), data=fuels.fit))  
plot(lm(fuels.lm3("(fuels_metric)^(1/3)"), data=fuels.fit))  

model.transform = "cube_root"
var.in = "fuels_metric^(1/3)"

# create formula
model.full = fuels.lm4(var.in)

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 3 predictors to minimize overfitting 
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 10, nvmax = 3, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=(1/3),
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=22
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric^(1/3)",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric^(1/3))) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric)^(1/3))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.s


### d. Sagebrush cover from proportion sagebrush

fuels.fit$prop_sagebrush = fuels.fit$Sagebrush_cover_pct/fuels.fit$Shrub_cover_pct 
# set to max of 1 - 1 plot where overlap bumped up
fuels.fit$prop_sagebrush = ifelse(fuels.fit$prop_sagebrush<=1, fuels.fit$prop_sagebrush, 1)

model.metric = "prop_sagebrush"

model.final=data.frame()
predictors.final = data.frame()

# step 1 is to rename the column so can use the same script in every section
fuels.fit = rename(fuels.fit,fuels_metric = model.metric)
head(fuels.fit)

# look at mini models for potential transformations, assumptions, outliers
boxcox(lm(fuels.lm1("(fuels_metric+0.01)"), data=fuels.fit))  # 
plot(lm(fuels.lm1("fuels_metric"), data=fuels.fit))  # 
boxcox(lm(fuels.lm2("(fuels_metric+1)"), data=fuels.fit))  # 
plot(lm(fuels.lm2("fuels_metric"), data=fuels.fit))  # 
boxcox(lm(fuels.lm3("(fuels_metric+0.01)"), data=fuels.fit))  # 
plot(lm(fuels.lm3("fuels_metric"), data=fuels.fit))  # 

model.transform = "none"
var.in = "(fuels_metric)"

# create formula
model.full = fuels.lm4(var.in)

# exhaustive best subsets selection
# look at top models using bic and r2
# narrowing down to 3 predictors to minimize overfitting 
model.subsets = regsubsets(model.full, data = fuels.fit, nbest = 20, nvmax = 3, method="exhaustive")

plot(model.subsets, scale="bic")
plot(model.subsets, scale="adjr2")

# summarize models that pass VIF, max correlation cutoffs
model.summ = summarize_models(power.in=1,
                              vif.cut=5, corr.cut=0.7)

arrange(model.summ,desc(adj_rsq))
arrange(model.summ, loocv_rmse)

# select based on best LOOCVrmse that meets linear model assumptions 

model.sel=55
model.notes="bestLOOCVRmse"

# check coefficients and relationships
names(coef(model.subsets,model.sel))[-1]
coef(model.subsets,model.sel)
cor((fuels.fit)[c(paste(names(coef(model.subsets,model.sel))[-1]))])
partial_R2("fuels_metric",names(coef(model.subsets,model.sel))[-1])

fuels.fit %>%
  dplyr::select(c("fuels_metric",names(coef(model.subsets,model.sel))[-1])) %>%
  pivot_longer(-fuels_metric) %>%
  ggplot(aes(x=value,y=fuels_metric)) +
  facet_wrap(~name,scales="free") +
  geom_point() +
  theme_bw()

# check model assumptions
model.form = as.formula(paste0(var.in,"~",paste(names(coef(model.subsets,model.sel))[-1], collapse="+")))
model.red = lm(model.form, data=fuels.fit)

boxcox(model.red)
# assumptions: check relatively equal variance, normality, minimize leverage, influence (cook's distance) < 1
plot(model.red, which=1)  # ok
qqPlot(model.red) # ok
plot(model.red, which=5)  # ok

# look at correlations among predictors and transformed fuels metric
fuels.cor = fuels.fit %>%
  mutate(fuels_metric=(fuels_metric))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))]))
abs(cor((fuels.cor)[c("fuels_metric",paste(names(coef(model.subsets,model.sel))[-1]))], method="spearman"))

# add to list of final models
model.final = rbind(model.final, data.frame(model.summ[model.summ$model_number==model.sel,],
                                            model_notes=model.notes))

# add to list of final predictors
predictors.final = rbind(predictors.final, summarize_predictors())

# examine
model.final
predictors.final

# add to master output

summary.out = rbind(summary.out,model.final)
predictors.out = rbind(predictors.out, predictors.final)

# revert back to original column name, needs to be updated
fuels.fit = fuels.s

# write current iteration
write.csv(summary.out, "analysis/fuels_prediction_map/model_final_selection.csv", row.names=FALSE)
write.csv(predictors.out, "analysis/fuels_prediction_map/model_final_predictors.csv", row.names=FALSE)
