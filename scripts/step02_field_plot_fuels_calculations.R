#####
#
## fuel characteristics calculations for 2019 fuels plots
#
#####

rm(list=ls())

# wd inherited from project

# load libraries
library(openxlsx) # version 4.1.3
library(dplyr) # version 0.8.3
library(tidyr) # version 1.0.2
library(ggplot2) # version 3.2.1
library(cowplot) # version 1.0.0

####
# 1. surface fuels
####

### load data
# plots
plots = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                  sheet="General_plot_measurements", colNames=TRUE)

# shrubs
shrub.sm = read.csv("data/Field_plots_2019/Cleaned_data/shrub_intercept_cover.csv")

quadrats.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                       sheet="Shrub_quadrats", colNames=TRUE)

shrub.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                     sheet="Shrub_intercepts", colNames=TRUE) %>%
  arrange(Plot_code, T, Intersect_begin)

# cwd
cwd.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                   sheet="CWD", colNames=TRUE) %>%
  arrange(Plot_code, T, Intersect_begin)

### shrub cover
# using shrub intercept
# reference: Sampling Vegetation Attributes (USDA FS rev. 1999)
shrub.cover = shrub.sm %>%
  # calc each intersect length
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code) %>%
  # sum by intersect and divide by total intersect length
  # which is 50 m
  summarise(Shrub_cover = (sum(Intersect_length)/50) * 100)

### shrub height
# weighted avg of all shrubs + trees < 2 m in height
# total overlapping cover
shrub.coverlap = shrub.in %>%
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code) %>%
  summarise(Overlap_cover = (sum(Intersect_length, na.rm=TRUE)/50)*100)
  
shrub.ht1 = shrub.in %>%
  # first calculate cover by species
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code,Species) %>%
  summarise(Spec_ht = mean(Height_m), Spec_cover = (sum(Intersect_length, na.rm=TRUE)/50)*100) %>%
  left_join(shrub.coverlap, by="Plot_code") %>%
  # species proportion of cover
  mutate(Prop_cover = Spec_cover/Overlap_cover) %>%
  # weighted average
  group_by(Plot_code) %>%
  summarise(Shrub_smtree_avg_ht = sum(Spec_ht * Prop_cover))

### also calculate cover for sagebrush specifically,
# in the interest of sage grouse habitat
# this will include any overlap among different sage species,
# but this overlap is usually minimal

sage.cover = shrub.in %>%
  filter(Species %in% c("ARCA","ARTRV","ARTRS2","ARTR4")) %>%
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code) %>%
  summarise(Sagebrush_cover_pct = (sum(Intersect_length, na.rm=TRUE)/50)*100)

sage.prop = shrub.in %>%
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code) %>%
  summarise(Total_cover_pct = (sum(Intersect_length, na.rm=TRUE)/50)*100) %>%
  left_join(sage.cover,by="Plot_code") %>%
  mutate(Prop_sagebrush = Sagebrush_cover_pct/Total_cover_pct) %>%
  dplyr::select(-c(Total_cover_pct:Sagebrush_cover_pct))

sage.ht = shrub.in %>%
  filter(Species %in% c("ARCA","ARTRV","ARTRS2","ARTR4")) %>%
  # first calculate cover by species
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code,Species) %>%
  summarise(Spec_ht = mean(Height_m), Spec_cover = (sum(Intersect_length, na.rm=TRUE)/50)*100) %>%
  left_join(shrub.coverlap, by="Plot_code") %>%
  # species proportion of cover
  mutate(Prop_cover = Spec_cover/Overlap_cover) %>%
  # weighted average
  group_by(Plot_code) %>%
  summarise(Sagebrush_ht_m = sum(Spec_ht * Prop_cover))


### shrub % dead
# using quadrats
# total of 8 .25-m2 quadrats for a total of 2-m2 per plot
# just averaging across them
quad.dead = quadrats.in %>%
  # add ht in proportion to cover for weighted avg
  mutate(Dead_sum = Pct_dead * Pct_cover) %>%
  group_by(Plot_code) %>%
  # Pct dead is weighted average of dead_sum divided by sum of total cover
  summarise(Sum_cover_q = sum(Pct_cover),
            Pct_dead = sum(Dead_sum, na.rm=TRUE)/Sum_cover_q)

### cwd cover
cwd.cover = cwd.in %>%
  # drop extra entries if multi intersect
  filter(!is.na(Intersect_end)) %>%
  # calc each intersect length
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code) %>%
  # sum by intersect and divide by total intersect length
  # which is 50 m
  summarise(Cwd_cover = (sum(Intersect_length)/50) * 100)

### cwd biomass (Mg/ha)
# source: Brown 1974

# 2 separate tallies, one for sound wood (classes 1-3) with specific gravity = 0.4,
# one for rotten (classes 4-5) with specific gravity = 0.3
# sound wood: w [tons/acre] = 4.656 * sum(dbh^2) [in^2] * c / N * l [ft]
# rotten wood: w = 3.492 * sum(dbh^2) * c / N * l
# where: w = weight of downed woody material
# dbh = dbh
# c = slope correction factor. c = sqrt(1 + (% slope/100)^2)
# N * l = length of sampling plane, with N as # of transects and l as length of each
# constant incorporates differences in specific gravity between sound and rotten

# I have converted these to metric units so that
# sound: w [Mg/ha] = 0.493 * sum(dbh^2) [cm^2] * c / N * l [m]
# rotten: w [Mg/ha] = 0.370 * sum(dbh^2) [cm^2] * c / N * l [m]
# this is constant * 0.155 [cm^2 to in^2] * 0.3048 [1/m to 1/ft] * 2.242 [Mg/ha to tons/acre]

cwd.biomass = plots %>%
  # first add slope from general plot measurements
  dplyr::select(Plot_code, Slope_deg) %>%
  right_join(cwd.in, by="Plot_code") %>%
  # remove unneeded columns, position not used but retained for now
  dplyr::select(-c(Intersect_begin,Intersect_end,Multi_intersect)) %>%
  # add sound or rotten
  mutate(Sound = ifelse(Class %in% c(1:3),TRUE,FALSE)) %>%
  # add dbh [in] and squared
  mutate(Diam_in = Diam_cm/2.54,
         Diam_cm2 = Diam_cm^2,
         Diam_in2 = Diam_in^2) %>%
  # sum dbh squared by plot
  group_by(Plot_code,Slope_deg,Sound) %>%
  summarise(Diam_cm2_sum = sum(Diam_cm2),
            Diam_in2_sum = sum(Diam_in2)) %>%
  # caculate slope correction factor
  mutate(Slope_rad = Slope_deg * pi/180,  # convert to radians
         Slope_pct = tan(Slope_rad) * 100,  # conver to percent
         c = sqrt(1 + (Slope_pct/100)^2)) %>%
  # add transect length, english conversions as double check
  mutate(Length_m = 50,
         Length_ft = 50/0.3048) %>%
  # add multiplier constant
  mutate(Sound_metric = ifelse(isTRUE(Sound),0.493,0.370),
         Sound_english = ifelse(isTRUE(Sound),4.656,3.492)) %>%
  # calculate CWD biomass in metric and english withough slope correction,
  # also metric with slope correction
  mutate(CWD_bm_met = (Sound_metric * Diam_cm2_sum)/Length_m,
         CWD_bm_met_slope = (Sound_metric * Diam_cm2_sum * c)/Length_m,
         CWD_bm_eng = (Sound_english * Diam_in2_sum)/Length_ft) %>%
  # add sound plus rotten
  group_by(Plot_code,Slope_deg,Slope_pct,c) %>%
  summarise_at(vars(CWD_bm_met:CWD_bm_eng), sum) %>%
  # add ck value for english to metric
  mutate(CWD_bm_ck = CWD_bm_eng * 2.242)
  
plot(CWD_bm_met~CWD_bm_ck, data=cwd.biomass)
abline(a=0,b=1, col="red")  # looks good

plot(CWD_bm_met_slope~CWD_bm_met, data=cwd.biomass)
abline(a=0,b=1, col="red")  # some minor differences

### FBFM

fbfm.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx", 
                    sheet="Surface_FBFMs", colNames = TRUE) %>%
  dplyr::select(Plot_code, Photo_FBFM) %>%
  rename(FBFM = Photo_FBFM)


####
# 2. stand structure and composition
####

# live trees, including df
trees.in = read.csv("data/Field_plots_2019/Cleaned_data/trees_cleaned_with_heights_grouped.csv") %>%
  filter(Status %in% c("L","DF"))

# plot size multiplier for trees
# multiply x 4 for Con_BerryGlade_Unburned3!!!
plot_size_ha = (pi * 13^2)/10000

# live saplings
saps.in = read.csv("data/Field_plots_2019/Cleaned_data/saplings_cleaned.csv")

# plot size multiplier for saps
# only measured on T1 and T3!!
sap_transect_size = ((2 * 26))/10000

# tree density, ba, dominant height (90th pctile)
trees.stand = trees.in %>%
  # add tph multiplier based on plot size
  mutate(tph_adj = ifelse(Plot_code != "Con_BerryGlade_Unburned3",(1/plot_size_ha),4*(1/plot_size_ha))) %>%
  # calc individ stem BA
  mutate(tree_ba = pi*(DBH_cm/200)^2) %>%
  group_by(Plot_code) %>%
  summarise(tph=sum(tph_adj),ba_m2=sum(tree_ba*tph_adj),
            dom_ht_m = quantile(Ht_ft_calc, probs=0.9)*0.3048,
            mean_tree_ht_m = mean(Ht_ft_calc)*0.3048)

trees.large = trees.in %>%
  # add tph multiplier based on plot size
  mutate(tph_adj = ifelse(Plot_code != "Con_BerryGlade_Unburned3",(1/plot_size_ha),4*(1/plot_size_ha))) %>%
  # filter for large trees only
  filter(DBH_cm >= 7.6) %>%
  group_by(Plot_code) %>%
  summarise(TPH_largetrees=sum(tph_adj))

# sapling density and mean height
saps.tph = saps.in %>%
  mutate(tph_adj = (1/sap_transect_size),
         tph_trans = Live_count * tph_adj) %>%
  group_by(Plot_code) %>%
  summarise(tph_saps=sum(tph_trans))

# mean sapling height
saps.ht = saps.in %>%
  dplyr::select(Plot_code, c(Height_m_1:Height_m_3)) %>%
  gather(key="Height_nbr",value="Height_m",-Plot_code, na.rm=TRUE) %>%
  group_by(Plot_code) %>%
  summarise(mean_sap_ht_m = mean(Height_m))

# canopy cover measured with densiometer
cc.densio = plots %>%
  dplyr::select(Plot_code, c(CanopyCover_T1:CanopyCover_T4)) %>%
  gather(key="CC_mean",value="CC_dots",-Plot_code) %>%
  group_by(Plot_code) %>%
  summarise(CC_pct = mean(CC_dots)*1.04)

####
# 3. compare canopy height and cover, FVS v. field, also vet densities, include deciduous
####

# read in updated data with corrected sapling count
fvs.fin = read.csv("processed_data/fuels_map_variables/fvs/outputs/trees_all_dflive_cbhcalc_6ftcutoff_outputs.csv") %>%
  # first entry only, they are identical
  filter(!duplicated(StandID)) %>%
  mutate(CFL_kg_m2 = Standing_Foliage * 0.22417,
         CBH_m = Canopy_Ht * 0.3048,
         SH_m = TopHt * 0.3048,
         TPH = Tpa * 2.47105,
         BA_m2 = (BA/10.764)*2.47105) %>%
  rename(CBD_kg_m3 = Canopy_Density, Plot_code = StandID, CC_pct = Total_Cover) %>%
  dplyr::select(Plot_code, CFL_kg_m2, CBD_kg_m3, CBH_m, SH_m, TPH, BA_m2, CC_pct, RunTitle)

# prep field data
field.comp = left_join(trees.stand,cc.densio) %>%
  left_join(saps.tph) %>%
  rename(tph_trees=tph,SH_m=dom_ht_m,BA_m2=ba_m2) %>%
  group_by(Plot_code) %>%
  mutate(TPH=sum(tph_trees,tph_saps,na.rm=TRUE),
         RunTitle="field") %>%
  dplyr::select(Plot_code,TPH,BA_m2,SH_m,CC_pct,RunTitle) %>%
  ungroup()

# prep fvs data
fvs.prep = fvs.fin %>%
  dplyr::select(Plot_code,TPH,BA_m2,SH_m,CC_pct,RunTitle) %>%
  mutate(Plot_code = as.character(Plot_code),
         CC_pct = as.numeric(CC_pct))

### combine and plot
stand.comp = rbind(field.comp,fvs.prep)

# plot boxplots
stand.comp %>%
  pivot_longer(c(TPH:CC_pct),names_to = "fuel_metric") %>%
  ggplot(aes(x=RunTitle,y=value, fill=RunTitle)) +
  facet_wrap(~fuel_metric, scales="free_y") +
  geom_boxplot() + 
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) 



# plot comparisons

# 1:1 comparisons, new function to include all plots
plot_fuelcomp2 = function(df1,df2,var_name) {
  # set equal xlim and ylim for easier comparison
  maxlim = df1 %>%
    left_join(df2, by="Plot_code") %>%
    select(c(paste0(var_name,".x"),paste0(var_name,".y"))) %>%
    max() * 1.05
  
  df1 %>%
    left_join(df2, by="Plot_code") %>%
    ggplot(aes_string(x=paste0(var_name,".x"),y=paste0(var_name,".y"))) +
    geom_point() + 
    geom_text(aes(label=Plot_code),hjust=-0.1,size=2) +
    geom_abline(col="red") +
    xlab(paste(df1$RunTitle[1],var_name)) +
    ylab(paste(df2$RunTitle[1],var_name)) +
    xlim(0,maxlim) +
    ylim(0,maxlim)+
    coord_cartesian(clip="off") +
    theme_bw()
}


plot_grid(plot_fuelcomp2(field.comp,fvs.prep,"TPH"),
          plot_fuelcomp2(field.comp,fvs.prep,"BA_m2"),
          plot_fuelcomp2(field.comp,fvs.prep,"SH_m"),
          plot_fuelcomp2(field.comp,fvs.prep,"CC_pct"),
          align="hv")


####
# 4. combine and visualize
####


### final fuels selection, see readme
fuels.out = plots %>%
  dplyr::select(Plot_code,Vegetation_type) %>%
  full_join(shrub.cover, by="Plot_code") %>%
  full_join(shrub.ht1, by="Plot_code") %>%
  full_join(quad.dead, by="Plot_code") %>%
  full_join(sage.cover, by="Plot_code") %>%
  # full_join(sage.prop, by="Plot_code") %>%
  full_join(sage.ht, by="Plot_code") %>%
  full_join(cwd.cover, by="Plot_code") %>%
  full_join(cwd.biomass, by="Plot_code") %>%
  full_join(fbfm.in, by="Plot_code") %>%
  full_join(trees.stand, by="Plot_code") %>%
  full_join(trees.large, by="Plot_code") %>%
  full_join(saps.tph, by="Plot_code") %>%
  full_join(saps.ht, by="Plot_code") %>%
  full_join(fvs.fin, by="Plot_code") %>%
  dplyr::select(-c(Sum_cover_q,Slope_deg,Slope_pct,c,CWD_bm_met,CWD_bm_eng,CWD_bm_ck,
                   dom_ht_m, mean_tree_ht_m,TPH,BA_m2,RunTitle)) %>%
  rename(Shrub_cover_pct = Shrub_cover, Shrub_ht_m = Shrub_smtree_avg_ht, Shrub_pctDead = Pct_dead,
         CWD_cover_pct = Cwd_cover,CWD_Mg_ha = CWD_bm_met_slope, TPH_trees = tph, BA_m2_ha = ba_m2, TPH_saps = tph_saps,
         Sap_ht_m = mean_sap_ht_m) %>%
  dplyr::select(c(Plot_code:Sap_ht_m,SH_m,CC_pct,CBH_m,CBD_kg_m3,CFL_kg_m2)) %>%
  # replace NAs with 0 where appropriate
  replace_na(list(Sagebrush_cover_pct=0,
                  CWD_cover_pct=0,
                  CWD_Mg_ha=0,
                  TPH_trees=0,
                  BA_m2_ha=0,
                  TPH_largetrees=0,
                  TPH_saps=0,
                  CC_pct=0))


str(fuels.out)
head(fuels.out)


fuels.out %>%
  # filter(Vegetation_type=="Conifer") %>%
  ggplot(aes(x=CC_pct, y=SH_m)) +
  ylab("Stand height (m)") +
  xlab("Canopy cover (%)") +
  geom_point(aes(color=factor(Vegetation_type))) +
  theme_bw()


fuels.out %>%
  dplyr::select(-FBFM) %>%
  pivot_longer(cols=c(Shrub_cover_pct:CFL_kg_m2),names_to="metric") %>%
  ggplot(aes(x=Vegetation_type, y=value)) +
  facet_wrap(~metric, scales="free") +
  # ylab("Stand height (m)") +
  # xlab("Canopy cover (%)") +
  geom_boxplot(aes(fill=factor(Vegetation_type))) +
  theme_bw() 


fuels.out %>%
  write.csv("processed_data/fuels_map_variables/field_plot_2019_fuels.csv", row.names=FALSE)
