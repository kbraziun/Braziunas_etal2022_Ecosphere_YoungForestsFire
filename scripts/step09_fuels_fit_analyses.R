#####
#
## summary stats, predicted v. observed for Q1
#
#####

rm(list=ls())

# load libraries

library(dplyr) # version 0.8.3
library(openxlsx) # version 4.1.3
library(tidyr) # version 1.0.2

####
# 1. load fuels data and create summary tables
####

plots.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                     sheet = "General_plot_measurements",
                     colNames=TRUE)

plot.summ = plots.in %>%
  mutate(Elevation = Postprocessed_elevation_m - 1.2) %>%  # gps unit was 1.2 m off the ground
  dplyr::select(c(Vegetation_type,Aspect_deg,Slope_deg,Elevation)) %>%
  pivot_longer(cols=c(Aspect_deg,Slope_deg,Elevation)) %>%
  group_by(Vegetation_type,name) %>%
  summarise_at("value",c(mean=mean,sd=sd,min=min,max=max),na.rm=TRUE)

fuels.in = read.csv("processed_data/fuels_map_variables/field_plot_2019_fuels.csv") %>%
  mutate(Sagebrush_cover_pct = ifelse(Sagebrush_cover_pct>Shrub_cover_pct,Shrub_cover_pct,Sagebrush_cover_pct),
         prop_sagebrush = Sagebrush_cover_pct/Shrub_cover_pct)

fuels.summ = fuels.in %>%
  dplyr::select(-c(Plot_code,FBFM)) %>%
  pivot_longer(cols=c(Shrub_cover_pct:prop_sagebrush)) %>%
  group_by(Vegetation_type,name) %>%
  summarise_at("value",c(mean=mean,sd=sd,min=min,max=max),na.rm=TRUE) %>%
  # remove columns that do not apply to certain vegetation types
  filter(!(Vegetation_type=="Conifer" & name %in% c("prop_sagebrush","Sagebrush_cover_pct","Sagebrush_ht_m")),
         !(Vegetation_type=="Deciduous" & name %in% c("prop_sagebrush","Sagebrush_cover_pct","Sagebrush_ht_m",
                                                      "CBD_kg_m3","CBH_m","CFL_kg_m2")),
         !(Vegetation_type=="Shrub" & name %in% c("CBD_kg_m3","CBH_m","CFL_kg_m2",
                                                      "CWD_cover_pct","CWD_Mg_ha",
                                                      "TPH_trees","BA_m2_ha","TPH_largetrees","TPH_saps","Sap_ht_m",
                                                      "SH_m","CC_pct")))
  

write.csv(plot.summ, "analysis/fuels_prediction_map/plot_summary_stats.csv", row.names=FALSE)
write.csv(fuels.summ, "analysis/fuels_prediction_map/fuels_summary_stats.csv", row.names=FALSE)

####
# 2. predicted versus observed
####

# load final models
model.in = read.csv("analysis/fuels_prediction_map/model_final_selection.csv")
model.pred = read.csv("analysis/fuels_prediction_map/model_final_predictors.csv")

# load predictors
lidar.in = read.csv("processed_data/fuels_map_variables/lidar_metrics.csv")
naip.in = read.csv("processed_data/fuels_map_variables/naip_metrics.csv")

# match up with fuels
fuels.pred = fuels.in %>%
  dplyr::select(c(Plot_code,Vegetation_type)) %>%
  left_join(lidar.in, by="Plot_code") %>%
  left_join(naip.in, by="Plot_code")

fuels.for = fuels.pred

fuels.shr = fuels.pred 

# loop through all predictors and use final model to predict fuels
# doing this as for loop to track progress
for(i in 1:dim(model.in)[1]) {
  print(i)
  
  model.form = model.pred %>%
    filter(fuels_metric == model.in[i,]$fuels_metric,
           veg_types == model.in[i,]$veg_types,
           model_number == model.in[i,]$model_number) %>%
    dplyr::select(c(pred,coef))
  
  fuels.met = model.in[i,]$fuels_metric
  
  print(fuels.met)
  
  # set up equation for each possible number of predictors
  
  if (model.in[i,]$n_pred == 1) {
    
    print("1 predictor")
    print(model.form)
    
    predi = as.character(model.form[1,]$pred)
    coefi = round(as.numeric(model.form[1,]$coef),4)
    pred1 = as.character(model.form[2,]$pred)
    coef1 = round(as.numeric(model.form[2,]$coef),4)
    
    print(paste0(coefi," + (",coef1,pred1,")"))
    
    fuels.pred[,paste0(fuels.met)] = coefi + (fuels.pred[,pred1] * coef1) 
    
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
    
    fuels.pred[,paste0(fuels.met)] = coefi + (fuels.pred[,pred1] * coef1) + (fuels.pred[,pred2] * coef2)
    
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
    
    fuels.pred[,paste0(fuels.met)] = coefi + (fuels.pred[,pred1] * coef1) + (fuels.pred[,pred2] * coef2) + (fuels.pred[,pred3] * coef3)
    
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
    
    fuels.pred[,paste0(fuels.met)] = coefi + (fuels.pred[,pred1] * coef1) + (fuels.pred[,pred2] * coef2) + (fuels.pred[,pred3] * coef3) + (fuels.pred[,pred4] * coef4)
    
    
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
    
    fuels.pred[,paste0(fuels.met)] = coefi + (fuels.pred[,pred1] * coef1) + (fuels.pred[,pred2] * coef2) + (fuels.pred[,pred3] * coef3) + (fuels.pred[,pred4] * coef4) + (fuels.pred[,pred5] * coef5)
    
  }
  
  # add to output
  if (model.in[i,]$veg_types %in% c("Conifer","Conifer_deciduous")) {
    fuels.for[,paste0(fuels.met)] = fuels.pred[,paste0(fuels.met)]
  } else if (model.in[i,]$veg_types == "Shrubland") {
    fuels.shr[,paste0(fuels.met)] = fuels.pred[,paste0(fuels.met)]
  }
  
}


# observed values
fuels.field = fuels.in %>%
  dplyr::select(-c(Sagebrush_ht_m,FBFM,TPH_largetrees,TPH_saps,Sap_ht_m)) %>%
  pivot_longer(c(Shrub_cover_pct:prop_sagebrush),names_to="variable",values_to="observed") 

# predicted values for forest
fuels.outf = fuels.for %>%
  dplyr::select(-c(zmax:NDVIhom)) %>%
  mutate(BA_m2_ha = BA_m2_ha^2,
         CBD_kg_m3 = CBD_kg_m3^3,
         CBH_m = CBH_m^2,
         CFL_kg_m2 = CFL_kg_m2^2,
         CWD_Mg_ha = exp(CWD_Mg_ha),
         TPH_trees = exp(TPH_trees)) %>%
  pivot_longer(cols=c(TPH_trees:Shrub_ht_m),names_to="variable",values_to="predicted") 

# predicted values for shrubland
fuels.outs = fuels.shr %>%
  filter(Vegetation_type=="Shrub") %>%
  dplyr::select(-c(zmax:NDVIhom)) %>%
  mutate(Shrub_pctDead = Shrub_pctDead^3) %>%
  pivot_longer(cols=c(Shrub_cover_pct:prop_sagebrush),names_to="variable",values_to="predicted") 

# join for comparison, remove NA fuels
# forest fuels, remove shrubland, remove canopy fuels from deciduous
fuels.compf = left_join(fuels.outf,fuels.field,by=c("Plot_code","Vegetation_type","variable")) %>%
  filter(Vegetation_type != "Shrub",
         !(Vegetation_type=="Deciduous" & variable %in% c("CBD_kg_m3","CBH_m","CFL_kg_m2",
                                                          "CWD_cover_pct","CWD_Mg_ha",
                                                          "Shrub_cover_pct","Shrub_ht_m")))


# shrubland
fuels.comps = left_join(fuels.outs,fuels.field,by=c("Plot_code","Vegetation_type","variable"))

# write out
write.csv(fuels.compf, "analysis/fuels_prediction_map/forest_predicted_observed.csv", row.names=FALSE)
write.csv(fuels.comps, "analysis/fuels_prediction_map/shrubland_predicted_observed.csv", row.names=FALSE)
