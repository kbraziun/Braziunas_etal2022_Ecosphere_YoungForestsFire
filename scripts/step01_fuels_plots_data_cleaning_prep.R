#####
#
## fuels plot data inspection, cleaning, and prep
#
#####

rm(list=ls())

# wd inherited from project

# load libraries
library(openxlsx) # version 4.1.3
library(tidyr) # version 1.0.2
library(dplyr) # version 0.8.3
library(ggplot2) # version 3.2.1
library(ggpubr) # version 0.2.4
library(cowplot) # version 1.0.0
library(qpcR) # version  1.4-1
library(RSQLite) # version 2.1.2

####
# 1. inspect gps data
####

gps.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx", 
                   sheet="GPS",
                   colNames=TRUE)

str(gps.in)

# compare field gps to postprocessed

gps_plot = function(x.name, y.name) {
  gps.in %>%
    ggplot(aes_string(x=x.name, y=y.name)) +
    facet_wrap(~Field_GPS) +
    geom_point() +
    geom_abline(slope=1, col="red") +
    theme_bw()
}

gps_plot("Field_easting","Postprocessed_easting")
gps_plot("Field_northing","Postprocessed_northing")
gps_plot("Field_elevation_m","Postprocessed_elevation_m")

# all look good
# elevation includes +1.26 m for where GPS was positioned above the ground
# but this is not an issue since elevation is not a predictor and the overall
# difference is negligible. this explains why it's consistently higher than
# points taken with the Oregon, though

####
# 2. general plot measurements
####

plots.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                     sheet="General_plot_measurements",colNames=TRUE)

# check that gps is matched up correctly
plots.gps = plots.in %>%
  left_join(gps.in, by=c("Plot_code","Postprocessed_easting","Postprocessed_northing","Postprocessed_elevation_m"))  # no warnings

# vegetation type
summary(as.factor(plots.in$Vegetation_type))

# flags
summary(as.factor(plots.in$Flag))  # 8 flagged
plots.in %>%
  filter(Flag=="Y") %>%
  mutate(Vegetation_type = as.factor(Vegetation_type)) %>%
  summary()  # mostly conifer

# canopy cover should be no more than 96, because out of 96
summary(plots.in)  # looks good

####
# 3. Surface FBFMs
####

# these will need to be revisited and recoded based on plot photos

# take random sample of plots for testing
plots.in %>%
  filter(Vegetation_type=="Conifer") %>%
  sample_n(4)

plots.in %>%
  filter(Vegetation_type=="Deciduous") %>%
  sample_n(1)

plots.in %>%
  filter(Vegetation_type=="Shrub") %>%
  sample_n(3)

####
# 4. shrubs
####

### load three shrub tables
ocular = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                 sheet="Ocular_shrubs_discontinued", colNames=TRUE)

quadrats = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                     sheet="Shrub_quadrats", colNames=TRUE)

summary(quadrats)

shrub.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                     sheet="Shrub_intercepts", colNames=TRUE) %>%
  arrange(Plot_code, T, Intersect_begin)

str(shrub.in)

### inspect shrub intercepts for 15 cm rule
# if errors identified in data recording, these
# were corrected in original xlsx data and noted on
# hard copy files
# this is consistent with data cleaning steps taken during
# data entry

# look at species
summary(as.factor(shrub.in$Species))


shrub.ck = shrub.in %>%
  # remove multi-intersects
  filter(!is.na(Intersect_end))

# go one species at a time
# iterate through until catch all errors

spec_list = unique(shrub.in$Species)

for(s in spec_list) {
  sp = s
  print(sp)
  
  # subset by species
  # and sort by plot, transect, intersect start location
  shrub.sub = shrub.ck %>%
    filter(Species == sp)
  
  # skip if only 1 entry
  if(dim(shrub.sub)[1] < 2) {
    next
  }
  
  # loop through each row, except last one because na
  for(i in 1:(dim(shrub.sub)[1]-1)) {
    # extract end of intersect in given row plus 15 cm
    # and beginning of intersect in next row
    i.end = round(shrub.sub[i,"Intersect_end"] + 0.15,2)
    i.next = round(shrub.sub[i+1,"Intersect_begin"],2)
    
    # only compare when in same plot and transect
    # if intersect + 15 cm overlaps, print error message
    if(shrub.sub[i,"T"] == shrub.sub[i+1,"T"] & 
       shrub.sub[i,"Plot_code"] == shrub.sub[i+1,"Plot_code"] &
       i.end > i.next) {
      print(paste("data record error at",sp,shrub.sub[i,"Plot_code"],shrub.sub[i,"T"],
                  shrub.sub[i,"Intersect_begin"]))
    }
  }
}

# go through each intercept with multiple entries
# and ensure have the same starting point

for(p in unique(shrub.in$Plot_code)) {
  
  print(p)
  # subset by species
  # and sort by plot, transect, intersect start location
  shrub.sub = shrub.in %>%
    filter(Plot_code==p) 
  
  int_list = unique(shrub.sub$Multi_intersect)
  int_list = int_list[int_list != "N"]
  
  if(length(int_list) == 0) {
    next
  }
  
  # loop through each multi-intersect
  for(i in 1:(length(int_list))) {
    
    # subset shrub.sub
    shrub.int = shrub.sub %>%
      filter(Multi_intersect==i)
  
    
    if(abs(max(shrub.int$Intersect_begin) - min(shrub.int$Intersect_begin)) != 0) {
      print(paste("data record error at",shrub.sub[i,"Plot_code"],"Multi_intersect",i))
    }
    
  }
}



### shrub intercept total cover
# have to create a version with all shrubs regardless of species
# to get total cover

# all shrubs, remove multi-intersects
# remove species ID
shrub.all = shrub.in %>%
  filter(!is.na(Intersect_end)) %>%
  dplyr::select(-c(Species,Multi_intersect,Height_m))

plot_list = unique(shrub.all$Plot_code)

# this loop will create a reduced dataframe
# merging intersections of different species
shrub.red = data.frame()

# loop through each row, except last one because na
for(p in 1:length(plot_list)) {
  
  shrub.sub = shrub.all[shrub.all$Plot_code==plot_list[p],]
  
  # loop through each transect
  for(t in unique(shrub.sub$T)) {
    
    shrub.tsub = shrub.sub[shrub.sub$T==t,]
    
    for(j in 1:10) {
      
      # reset end value
      i.newend = 0
      
      # get transect length
      t.length = dim(shrub.tsub)[1]
      
      # blank dataframe
      shrub.trans = data.frame()
      
      # use T/F to indicate if need to continue
      # looping through
      continue=FALSE
      
      for(i in 1:t.length) {
        
        # extract information from previous,
        # current, and next intersect
        i.prevend = ifelse(i>1,round(shrub.tsub[i-1,"Intersect_end"],2),0)
        i.beg = round(shrub.tsub[i,"Intersect_begin"],2)
        i.end = round(shrub.tsub[i,"Intersect_end"],2)
        i.next = ifelse(i<t.length,round(shrub.tsub[i+1,"Intersect_begin"],2),26)
        i.nextend = ifelse(i<t.length,round(shrub.tsub[i+1,"Intersect_end"],2),0)
        
        # does current intersect begin before previous one ends
        overlap.prev = ifelse(i>1 & i.beg < i.prevend, TRUE, FALSE)
        
        # catch possible problem: does current intersect end after
        # the new ending point (updated end of each loop interation)
        keep = ifelse(i.end > i.newend, TRUE, FALSE)
        
        # check if overlaps next intersect
        overlap.next = ifelse(i.end > i.next, TRUE, FALSE)
        
        # update continue to TRUE if at least one intersect overlaps
        if(i.end > i.next) {
          continue=TRUE
        }
        
        # if overlap.prev is true, go to next intersect
        # however, retain intersect if it would extend the previously
        # recorded intersect. situation that is likely if there are
        # lots of lots of intersections recorded in succession that
        # continually lengthen one another
        if(isTRUE(overlap.prev) & isFALSE(keep)) {
          next
        }
        
        # update shrub.red with correct info
        # if overlap,
        # use maximum of the intersect end of this and subsequent intersect
        i.newend = ifelse(isTRUE(overlap.next),max(i.end,i.nextend),i.end)
        
        # output updated intersect
        shrub.trans = rbind(shrub.trans, 
                          data.frame(Plot_code = plot_list[p],
                                     T=t,
                                     Intersect_begin = i.beg,        
                                     Intersect_end = i.newend))
  
      }
      
      # update df for rerunning
      shrub.tsub = shrub.trans
      
      # exit loop if no overlaps were found
      if(isFALSE(continue)) {break}
    }
    

    # update output dataframe
    shrub.red = rbind(shrub.red,shrub.trans)
    
    # all overlaps caught?
    if(isTRUE(continue)) {
      print("some overlaps remain")
    }
  }
}

head(shrub.red)
summary(shrub.red)

# compare with manual edits to data
shrub.red %>%
  filter(Plot_code=="Con_24_1")
shrub.red %>%
  filter(Plot_code=="Con_93_1")
# looks good!

### compare 3 methods of estimating shrub cover
# just compare total cover and average shrub height
# recall that for ocular and quadrat data, avg ht is 70% of max ht

# ocular cover
head(ocular)
ocular.tot = ocular %>%
  # add ht in proportion to cover for weighted avg
  mutate(Ht_prop = Avg_ht_m * Pct_cover) %>%
  group_by(Plot_code) %>%
  # sum to total cover, avg height
  summarise(Tot_cover_o = sum(Pct_cover),
            Avg_ht_o = sum(Ht_prop)/Tot_cover_o)

# quadrat cover
head(quadrats)
# total of 8 .25-m2 quadrats for a total of 2-m2 per plot
# but at this point, just averaging is good
quad.tot = quadrats %>%
  # add ht in proportion to cover for weighted avg
  mutate(Ht_prop = Avg_ht_m * Pct_cover) %>%
  group_by(Plot_code) %>%
  # Total cover is mean of cover values across quadrats
  # Avg height is weighted average of height divided by sum of total cover
  summarise(Tot_cover_q = mean(Pct_cover),
            Sum_cover_q = sum(Pct_cover),
            Avg_ht_q = sum(Ht_prop, na.rm=TRUE)/Sum_cover_q)


# shrub intercept cover
shrub.cover = shrub.red %>%
  # calc each intersect length
  mutate(Intersect_length = Intersect_end - Intersect_begin) %>%
  group_by(Plot_code) %>%
  # sum by intersect and divide by total intersect length
  # which is 50 m
  summarise(Tot_cover_int = (sum(Intersect_length)/50) * 100)

# shrub intercept average ht
shrub.ht = shrub.in %>%
  group_by(Plot_code) %>%
  summarise(Avg_ht_int = quantile(Height_m, probs=c(0.7)))

# join ocular, quadrat, intercept
shrub.comp = shrub.cover %>%
  left_join(shrub.ht, by="Plot_code") %>%
  left_join(quad.tot, by="Plot_code") %>%
  left_join(ocular.tot, by="Plot_code") 

head (shrub.comp) # NAs for lots of ocular data

# generic shrub comparison plot
shrub_plot = function(y.val, x.val,y.lab,x.lab,low.in,hi.in) {
  shrub.comp %>%
    ggplot(aes_string(y=y.val, x=x.val)) +
    geom_point(size=1,na.rm=TRUE) +
    # 1:1 line
    geom_abline(slope=1, col="red") +
    # pearson's correlation coefficient
    stat_cor(method="pearson", na.rm=TRUE) +
    xlim(low.in,hi.in) +
    ylim(low.in,hi.in) +
    ylab(y.lab) +
    xlab(x.lab) +
    theme_bw() +
    theme(plot.margin = margin(c(5,10,10,10)))
  
}

# how well does ocular do?
s1 = shrub_plot("Tot_cover_o","Tot_cover_q","Ocular shrub cover (%)","Quadrat shrub cover (%)",0,80)
s2 = shrub_plot("Tot_cover_o","Tot_cover_int","Ocular shrub cover (%)","Intercept shrub cover (%)",0,80)
s3 = shrub_plot("Avg_ht_o","Avg_ht_q","Ocular avg shrub height (m)","Quadrat avg shrub height (m)",0,1.6)
s4 = shrub_plot("Avg_ht_o","Avg_ht_int","Ocular avg shrub ht (m)","Intercept avg shrub height (m)",0,1.6)

# how do quadrat and intercept compare?
s5 = shrub_plot("Tot_cover_q","Tot_cover_int","Quadrat shrub cover (%)","Intercept shrub cover (%)", 0,80)
s6 = shrub_plot("Avg_ht_q","Avg_ht_int","Quadrat avg shrub height (m)","Intercept avg shrub height (m)",0,1.6)

### compare shrub cover methods

plot_grid(s1,s2,s5,s4,s3,s6, nrow=2)


### write output
# write modified shrub cover as csv

write.csv(shrub.red, file="data/Field_plots_2019/Cleaned_data/shrub_intercept_cover.csv", row.names=FALSE)

####
# 5. CWD
####

# load table
cwd.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                   sheet="CWD", colNames=TRUE) %>%
  arrange(Plot_code, T, Intersect_begin)

head(cwd.in)
summary(cwd.in)

# ensure no missed intersection overlaps

cwd.ck = cwd.in %>%
  # remove multi-intersects
  filter(!is.na(Intersect_end))


# loop through each row, except last one because na
for(i in 1:(dim(cwd.ck)[1]-1)) {
  # extract end of intersect in given row plus 15 cm
  # and beginning of intersect in next row
  i.end = round(cwd.ck[i,"Intersect_end"],2)
  i.next = round(cwd.ck[i+1,"Intersect_begin"],2)
  
  # only compare when in same plot and transect
  # if intersect + 15 cm overlaps, print error message
  if(cwd.ck[i,"T"] == cwd.ck[i+1,"T"] & 
     cwd.ck[i,"Plot_code"] == cwd.ck[i+1,"Plot_code"] &
     i.end > i.next) {
    print(paste("data record error at",cwd.ck[i,"Plot_code"],cwd.ck[i,"T"],
                cwd.ck[i,"Intersect_begin"]))
  }
}  # good!

# go through each intercept with multiple entries
# and ensure have the same starting point

for(p in unique(cwd.in$Plot_code)) {
  
  print(p)
  # subset by plot
  cwd.sub = cwd.in %>%
    filter(Plot_code==p) 
  
  int_list = unique(cwd.sub$Multi_intersect)
  int_list = int_list[int_list != "N"]
  
  if(length(int_list) == 0) {
    next
  }
  
  # loop through each multi-intersect
  for(i in 1:(length(int_list))) {
    
    # subset shrub.sub
    cwd.int = cwd.sub %>%
      filter(Multi_intersect==i)
    
    
    if(abs(max(cwd.int$Intersect_begin) - min(cwd.int$Intersect_begin)) != 0) {
      print(paste("data record error at",cwd.sub[i,"Plot_code"],"Multi_intersect",i))
    }
    
  }
} # good!

# examine whether cover estimates might be low due to 
# measurement methods issues
# i.e., when intercept is less than diameter
# only using first diameter if multi-intersect
head(cwd.ck)

cwd.dia = cwd.ck %>%
  mutate(Intersect_length = (Intersect_end - Intersect_begin) * 100) %>%
  mutate(Intersect_ratio = Intersect_length/Diam_cm)

summary(cwd.dia[cwd.dia$Intersect_ratio<0.8,])  # 17 with especially problematic measurements

cwd.plots = cwd.dia %>%
  group_by(Plot_code) %>%
  summarise(Intersect_ratio_mean = mean(Intersect_ratio),
          Intersect_ratio_median = median(Intersect_ratio))

cwd.plots
summary(cwd.plots)
# when aggregated to plot level, initial glance is that looks like
# nothing is completely off at entire plot level
# so moving on for now

####
# 6. Trees
####

trees.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                  sheet="Trees", colNames = TRUE)

head(trees.in)
summary(trees.in)

### step 1: ID issues with data

unique(trees.in$Plot_code)  # looks good

# trees with unknown status
unique(trees.in$Status)
trees.in %>%
  filter(Status=="Unknown")  # handful of trees where not recorded
# none of these will be used in fitting regression
# RULE 1: set trees with unknown status to live (only a few)

# trees with unknown structural stage
unique(trees.in$Struct_stage)
trees.in %>%
  filter(is.na(Struct_stage), Status != "D")
# all DF trees need structural stage
trees.in %>%
  filter(is.na(Struct_stage), Status == "L")
trees.in %>%
  filter(Struct_stage=="Unknown", Status == "L")
# will need to assign structural stage for some live trees
# RULE 2: assign trees with missing structural stage to structural
# stage of most similarly-sized tree based on DBH. only important
# for fuelcalc, not used in FVS

# split stems
trees.in %>%
  filter(is.na(DBH_cm)) # one tree with no DBH measurement
# this is a split stem above DBH, will have to address
# RULE 3: if step is split above dbh, consolidate to one record by choosing
# tallest height (largest dbh) and dropping shorter (smaller) one

unique(trees.in$Multi_stem)
# how many trees
multis = trees.in %>%
  group_by(Plot_code) %>%
  summarise(multi_stems = max(as.numeric(Multi_stem),na.rm=TRUE)) %>% 
  filter(!is.infinite(multi_stems))

sum(multis$multi_stems)  # 113, minus 2 will be dropped because split above dbh
# no need to change any info for multistems split below dbh, just treat them as trees
# NO ACTION NEEDED: count multi-stems that split below DBH as individual trees

# any trees with cbh > ht, indicates problem with data recording or entry
trees.in %>%
  filter(CBH_m>=Height_m) 
# for all of these, reviewing raw data sheets
# if carries through to data sheets, assume these were switched by
# accident while recording in field and fix here
# all of these from data sheets
# fix below
# RULE 4a: Fix these but still omit these trees from regression calculations
# RULE 4b: For one tree where noted that likely height estimated, switch to
# expected height value rather than switching height and CBH

unique(trees.in$Species)  # some unknown
trees.in %>%
  filter(Species=="Unknown")  # these are dead
# RULE 5: set unknown species (all dead) to most common species in that stand
# or in case of two that are "maybe POBA", set to POBA


trees.in %>%
  filter(is.na(Height_m), DBH_cm >= 7.6)  # trees where height should have been recorded
# one tree exactly 7.6, height and cbh not recorded, just include in regression
# some dead, OK for now
# some bent or leaning, but only a couple live, leaving these in the mix
# NO ACTION NEEDED: if CBH or height missing, assign by regression

trees.in %>%
  filter(is.na(CBH_m), DBH_cm >=7.6, Status != "D")
trees.in %>%
  filter(Status == "DF")
# DF trees will need CBH calc'd
# will need to examine flagged trees on case by case basis
# RULE 14: consider DF trees as live during time of LiDAR flight
# and determine CBH, although exclude these trees from fitting regressions


unique(trees.in$CR_field)  # if assigned in field, make sure used
trees.in %>%
  filter(CR_field==0)  # one tree live but with "no foliage"
# RULE 6: trees with CR_field recorded should be assigned a damage code of 97 (dead top)
# and then use CR_field to calculate HtTopK (height to top kill)


unique(trees.in$Exclude)  # Y, N, and Flag
trees.in %>%
  filter(Exclude=="Y") %>%
  tally()
unique(trees.in$Exclude_reason)  # these will be excluded from height regression
# NO ACTION NEEDED: if leaning or bent, leave be

unique(trees.in$Flag)  # Y, N
trees.in %>%
  filter(Flag=="Y") %>%
  tally()
flagged = trees.in %>%
  filter(Flag=="Y") 
unique(trees.in$Notes)  # review for flagged trees
# NO ACTION NEEDED: half in half out; rooted outside but crown in; etc leave in
# NO ACTION NEEDED: irregular crown but no CR_field, leave be
# NO ACTION NEEDED: DBH incorrect, this one is dead so leave be
# RULE 7: if rooted in but crown mostly or all out, remove (2 entries)
# RULE 8: if no CBH, possibly 0; no CBH recorded because bent over; set to 0 (2 entries)
# RULE 9: if broken top, set damage code to 96 (broken/missing top), set HtTopK to recorded
# height or leave HtTopK blank if not recorded (FVS will calc this then as 80% of total height,
# which would be calculated from the height regression)
# RULE 10: CBH appears incorrect, erased; exclude from regression but otherwise leave be
# RULE 11: CBH recorded but status is dead: remove CBH (1 entry)
# RULE 12: Set POBA to POTR
# RULE 13: height regression

### RULE and NO ACTION NEEDED summary
# RULE 1: set trees with unknown status to live (only a few)
# RULE 2: assign trees with missing structural stage to structural
# stage of most similarly-sized tree based on DBH. only important
# for fuelcalc, not used in FVS
# RULE 3: if step is split above dbh, consolidate to one record by choosing
# tallest height (largest dbh) and dropping shorter (smaller) one
# RULE 4a: Fix these but still omit these trees from regression calculations
# RULE 4b: For one tree where noted that likely height estimated, switch to
# expected height value rather than switching height and CBH
# RULE 5: set unknown species (all dead) to most common species in that stand
# or in case of two that are "maybe POBA", set to POBA
# RULE 6: trees with CR_field recorded should be assigned a damage code of 97 (dead top)
# and then use CR_field to calculate HtTopK (height to top kill)
# RULE 7: if rooted in but crown mostly or all out, remove (2 entries)
# RULE 8: if no CBH, possibly 0; no CBH recorded because bent over; set to 0 (2 entries)
# RULE 9: if broken top, set damage code to 96 (broken/missing top), set HtTopK to recorded
# height or leave HtTopK blank if not recorded (FVS will calc this then as 80% of total height,
# which would be calculated from the height regression)
# RULE 10: CBH appears incorrect, erased; exclude from regression but otherwise leave be
# RULE 11: CBH recorded but status is dead: remove CBH (1 entry)
# RULE 12: Set POBA to POTR prior to height regression, carry through to FVS
# RULE 13: height regression by species, only for trees that are not excluded from
# regression, use to set height and CBH for all trees in which these were not recorded
# including live and dead trees
# RULE 14: experiment with DF trees. could consider DF trees as live during time of LiDAR flight
# or dead to better reflect current fuel conditions. determine what to go with based
# on fuels regression fit

# NO ACTION NEEDED: count multi-stems that split below DBH as individual trees
# NO ACTION NEEDED: if CBH or height missing, assign by regression
# NO ACTION NEEDED: if leaning or bent, leave be
# NO ACTION NEEDED: half in half out; rooted outside but crown in; etc leave in
# NO ACTION NEEDED: irregular crown but no CR_field, leave be
# NO ACTION NEEDED: DBH possibly incorrect flag, this one is dead so leave be

### implement rules and fit regression for assigning height and CBH to
# small trees

# RULE 1: trees with missing status
# assume live if not recorded
trees.in %>%
  filter(Status=="Unknown") %>%
  tally()  # this affects 9 trees

trees.out = trees.in
trees.out$Status = ifelse(trees.in$Status=="Unknown","L",trees.in$Status)

summary(as.factor(trees.in$Status))
summary(as.factor(trees.out$Status))

# RULE 2: trees with missing structural stage
# includes some live trees with error in data recording
# and dead with foliage trees where this was NA
# NOTE: do not need structural stage for FVS

trees.out %>%
  filter(Struct_stage == "Unknown" | is.na(Struct_stage), Status != "D") %>%
  group_by(Status) %>%
  tally()  # 63 DF, 10 L

stage.lookup = trees.out %>%
  filter(Status=="L") %>%
  group_by(Plot_code, Struct_stage) %>%
  summarise(mean_dbh = mean(DBH_cm, na.rm=TRUE))

# loop through all trees
for(i in 1:dim(trees.out)[1]) {
  
  if(trees.out[i,]$Struct_stage %in% c("C","D","I","O","S") | trees.out[i,]$Status == "D") {
    next
  } else {
    
    # get lookup for this plot
    stage.sub = stage.lookup %>%
      filter(Plot_code==trees.out[i,]$Plot_code, Struct_stage %in% c("C","D","I","O","S"))
    
    rownum = which(abs(stage.sub$mean_dbh-trees.out[i,]$DBH_cm)==min(abs(stage.sub$mean_dbh-trees.out[i,]$DBH_cm)))
    
    # structural stage
    trees.out[i,]$Struct_stage = stage.sub[rownum,]$Struct_stage
    
    # print(trees.out[i,]$Struct_stage)
    # print(stage.sub[rownum,])
  }
  

  }

trees.out %>%
  filter(Struct_stage == "Unknown" | is.na(Struct_stage), Status != "D") %>%
  group_by(Status) %>%
  tally()  # all taken care of

summary(as.factor(trees.in$Struct_stage))
summary(as.factor(trees.out$Struct_stage))

# RULE 3: split trees
# if step is split above dbh, consolidate to one record by choosing
# tallest height (largest dbh) and dropping shorter (smaller) one
trees.out %>%
  filter(Multi_stem != "N", !is.na(Notes))

# use tallest height
trees.out[trees.out$Plot_code=="Con_13_2" & trees.out$Q==2 & trees.out$Multi_stem=="3" & !is.na(trees.out$DBH_cm),]$Height_m = 16.5

# drop 2 extra stems
trees.backup = trees.out
trees.temp = trees.backup %>%
  filter(!(Plot_code=="Con_13_2" & Q==2 & Multi_stem=="3" & is.na(DBH_cm)))
trees.out = trees.temp %>%
  filter(!(Plot_code=="Con_31_1" & Q==1 & Multi_stem=="1" & DBH_cm==0.7))

trees.out %>%
  filter(Multi_stem != "N", !is.na(Notes)) # looks good

# RULE 4a,b: if CBH is greater than height, assume
# recorded incorrectly and fix
# but also exclude these from regression calcs
trees.out %>%
  filter(CBH_m >= Height_m) 

summary(trees.out$Height_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.70    7.30   10.15   11.32   14.00   43.00    3143 

trees.out %>%
  group_by(Exclude) %>%
  tally()

# 149, 4184, 28

# first exclude all these trees
trees.out$Exclude = ifelse(!is.na(trees.out$CBH_m) & trees.out$CBH_m>trees.out$Height_m, "Y", trees.out$Exclude) # good

# backup height and cbh
trees.out$Height_backup = trees.out$Height_m
trees.out$CBH_backup = trees.out$CBH_m

# first fix trees where it is probably flipped
trees.out$Height_m = ifelse(!is.na(trees.out$CBH_m) & trees.out$CBH_m>trees.out$Height_m, trees.out$CBH_backup, trees.out$Height_backup)
trees.out$CBH_m = ifelse(!is.na(trees.out$CBH_m) & trees.out$CBH_backup>trees.out$Height_backup, trees.out$Height_backup, trees.out$CBH_backup)

# then fix tree where it is probably 21.5
trees.out$Height_m = ifelse(!is.na(trees.out$CBH_m) & trees.out$CBH_backup >= trees.out$Height_backup & trees.out$Plot_code=="Con_42_1", 21.5, trees.out$Height_m)

trees.out %>%
  filter(CBH_m >= Height_m)  # addressed
trees.out %>%
  filter(CBH_backup >= Height_backup)  # addressed

summary(trees.out$Height_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.70    7.30   10.20   11.34   14.00   43.00    3143 


# RULE 5: set unknown species (all dead) to most common species in that stand
# or in case of two that are "maybe POBA", set to POBA
trees.out %>%
  filter(Species=="Unknown")  # these are dead

# Create the function to get mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# get most common species in lookup table
# "maybe POBA" trees are in stand dominated by POBA
domspec.lookup = trees.out %>%
  group_by(Plot_code) %>%
  summarise(dom_species = getmode(Species))

# join, substitute unknown species
trees.backup2 = trees.out

trees.out = trees.backup2 %>%
  left_join(domspec.lookup, by="Plot_code") %>%
  mutate(Species = ifelse(Species=="Unknown",dom_species,Species))

trees.backup2 %>%
  filter(Species=="Unknown")  
trees.out %>%
  filter(Species=="Unknown")  # addressed



# RULE 6: trees with CR_field recorded should be assigned a damage code of 97 (dead top)
# and then use CR_field to calculate HtTopK (height to top kill)
unique(trees.out$CR_field)  # if assigned in field, make sure used
trees.out %>%
  filter(CR_field != "N")
trees.out %>%
  filter(CR_field == "N")

trees.out$damage_code = ifelse(trees.out$CR_field !="N",97,NA)
trees.out$CR_m = ifelse(trees.out$CR_field !="N",trees.out$Height_m * as.numeric(trees.out$CR_field)/100,NA)
trees.out$Height_topkill_m = ifelse(trees.out$CR_field !="N",round(trees.out$CBH_m + trees.out$CR_m,1),NA)


# RULE 7: if rooted in but crown mostly or all out, remove (2 entries)
unique(trees.out$Notes)  # review for flagged trees

trees.backup3 = trees.out
trees.out %>%
  group_by(Status) %>%
  tally() # 565 D, 63 DF, 3733 L

trees.out = trees.out %>%
  filter(!Notes %in% c("Rooted in, but crown out","Half out, crown mostly out"))

trees.out %>%
  group_by(Status) %>%
  tally() # 563 D, 63 DF, 3733 L


# RULE 8: if no CBH, possibly 0; no CBH recorded because bent over; set to 0 (2 entries)
summary(trees.out$CBH_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.90    2.60    3.64    5.60   18.30    3276 

trees.out$CBH_m = ifelse(trees.out$Notes %in% c("No CBH recorded, possibly 0", "No CBH recorded, probably because bent over"), 0, trees.out$CBH_m)

summary(trees.out$CBH_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.900   2.600   3.633   5.600  18.300    3274


# RULE 9: if broken top, set damage code to 96 (broken/missing top), set HtTopK to recorded
# height or leave HtTopK blank if not recorded (FVS will calc this then as 80% of total height,
# which would be calculated from the height regression)
unique(trees.out$Exclude_reason)

trees.out %>%
  group_by(damage_code) %>%
  tally()  # 14 97, 4345 NA

trees.out$damage_code = ifelse(trees.out$Exclude_reason %in% c("Broken top", "Broken top (inferred)","Snapped, bent to ground","Broken top, still attached","Possible broken top"), 96, trees.out$damage_code)

trees.out %>%
  group_by(damage_code) %>%
  tally()  # 48 96, 14 97, 4297 NA

summary(trees.out$Height_topkill_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.800   4.400  12.100   9.429  13.250  16.500    4345 

trees.out$Height_topkill_m = ifelse(trees.out$damage_code==96 & !is.na(trees.out$Height_m),trees.out$Height_m,trees.out$Height_topkill_m)

summary(trees.out$Height_topkill_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   1.400   2.350   4.000   6.281  10.350  16.500    4316


# RULE 10: CBH appears incorrect, erased; exclude from regression but otherwise leave be
trees.out %>%
  filter(Notes=="CBH may be incorrect, was erased")

trees.out %>%
  group_by(Exclude) %>%
  tally() # 148 F, 4179 N, 32 Y

trees.out$Exclude = ifelse(trees.out$Notes %in% c("CBH may be incorrect, was erased"),"Y",trees.out$Exclude)

trees.out %>%
  group_by(Exclude) %>%
  tally() # 148 F, 4178 N, 33 Y


# RULE 11: CBH recorded but status is dead: remove CBH (1 entry)
trees.out %>%
  filter(Notes=="CBH recorded but status is dead" )

summary(trees.out$CBH_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.900   2.600   3.633   5.600  18.300    3274

trees.out$CBH_m = ifelse(trees.out$Notes %in% c("CBH recorded but status is dead"),NA,trees.out$CBH_m)

summary(trees.out$CBH_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.900   2.550   3.634   5.600  18.300    3275 


# RULE 12: Set POBA to POTR prior to height regression, carry through to FVS
trees.out %>%
  group_by(Species) %>%
  tally()
# 1 ABLA     1132
# 2 PIAL       65
# 3 PICO     1826
# 4 PIEN      160
# 5 POBA        8
# 6 POTR     1054
# 7 PSME      114

trees.out$Species = ifelse(trees.out$Species=="POBA","POTR",trees.out$Species)

trees.out %>%
  group_by(Species) %>%
  tally()
# Species     n
# <chr>   <int>
#   1 ABLA     1132
# 2 PIAL       65
# 3 PICO     1826
# 4 PIEN      160
# 5 POTR     1062
# 6 PSME      114


# RULE 13a: height regression by species, only for trees that are not excluded from
# regression, use to set height for all trees in which these were not recorded
# including live and dead trees
# RULE 13b: use default FVS equations to set crown ratio for all trees where CBH not recorded,
# which takes stand density and species-specific crown length into account
# fit separate regressions for each species

# compare model fit using FVS defaults vs. my data
fvs.params = data.frame(Species = c("PIAL","PSME","POTR","PICO","PIEN","ABLA"),
                        fvs.b1 = c(4.1920,4.5175,4.4625,4.4625,4.5822,4.3603),
                        fvs.b2 = c(-5.1651,-6.5129,-5.2223,-5.2223,-6.4818,-5.2148))


# fit regressions for English units to enable comparison with FVS
trees.ht = trees.out %>%
  filter(Status=="L", !is.na(Height_m), Exclude == "N") %>%
  mutate(DBH_in=round(DBH_cm/2.54,2), 
         Ht_ft = round(Height_m/0.3048,2),
         CBH_ft = round(CBH_m/0.3048,2)) %>%
  left_join(fvs.params,by="Species")
# n = 1049

models.out = data.frame()

# using logistic equation of same form as fvs
# use similar form for cbh for consistency
# however, RMSE is improved by fitting to field measured data rather than using
# default FVS equations
for(i in unique(trees.ht$Species)) {

  trees.sp = trees.ht %>%
    filter(Species==i)

  model.ht = nls(Ht_ft~4.5+exp(b1+b2/(DBH_in+1)), data=trees.sp, start=list(b1=1,b2=-1))
  fit.b1 = summary(model.ht)$coefficients["b1",1]
  fit.b2 = summary(model.ht)$coefficients["b2",1]
  fit.rmse = RMSE(model.ht)

  # ck rmse against rmse if using FVS stock equation, improvement
  trees.sp$pred = 4.5+exp(trees.sp$fvs.b1+(trees.sp$fvs.b2)/(trees.sp$DBH_in+1))
  trees.sp$error = trees.sp$Ht_ft - trees.sp$pred
  trees.sp$err2 = trees.sp$error^2
  fvs.rmse = sqrt(mean(trees.sp$err2))

  # # ck rmse. looks good.
  # trees.sp$pred = 4.5+exp(fit.b1+fit.b2/(trees.sp$DBH_in+1))
  # trees.sp$error = trees.sp$Ht_ft - trees.sp$pred
  # trees.sp$err2 = trees.sp$error^2
  # sqrt(mean(trees.sp$err2))

  model.cbh = nls(CBH_ft~exp(b1+b2/(DBH_in+1)), data=trees.sp, start=list(b1=1,b2=-1))
  cbh.b1 = summary(model.cbh)$coefficients["b1",1]
  cbh.b2 = summary(model.cbh)$coefficients["b2",1]
  cbh.rmse = RMSE(model.cbh)

  model.n = length(trees.sp[,1])

  models.out = rbind(models.out, data.frame(Species = i,b1_ht = fit.b1,b2_ht=fit.b2,
                                            rmse_ht=fit.rmse,rmse_fvsht=fvs.rmse,
                                            b1_cbh = cbh.b1, b2_cbh=cbh.b2, rmse_cbh = cbh.rmse, n_trees=model.n))
}


models.out  # RMSE similar between ht and cbh. plots below make height look like better fit
str(models.out)

trees.ht %>%
  left_join(models.out, by="Species") %>%
  ggplot(aes(x=DBH_in, y=Ht_ft)) +
  facet_wrap(~Species) +
  geom_point() +
  geom_line(aes(y=4.5+exp(b1_ht + b2_ht/(DBH_in+1))), color="red",size=1) +
  geom_line(aes(y=4.5+exp(fvs.b1 + (fvs.b2)/(DBH_in+1))), color="blue",size=1) +
  theme_bw() 

trees.ht %>%
  left_join(models.out, by="Species") %>%
  ggplot(aes(x=DBH_in, y=CBH_ft)) +
  facet_wrap(~Species) +
  geom_point() +
  geom_line(aes(y=exp(b1_cbh + b2_cbh/(DBH_in+1))), color="red",size=1) +
  # geom_line(aes(y=4.5+exp(fvs.b1 + (fvs.b2)/(DBH_in+1))), color="blue",size=1) +
  theme_bw() 

# calculate height for all trees
# leave out CBH calculations, allow FVS to fill in
trees.calc = trees.out %>%
  dplyr::select(-c(Height_backup:dom_species,CR_m)) %>%
  left_join(models.out, by="Species") %>%
  mutate(DBH_in=round(DBH_cm/2.54,2), 
         Ht_ft = round(Height_m/0.3048,2),
         CBH_ft = round(CBH_m/0.3048,2)) %>%  
  mutate(Predicted = ifelse(is.na(Ht_ft),"yes","no")) %>%
  mutate(Ht_ft_calc = ifelse(is.na(Ht_ft), round(4.5+exp(b1_ht + b2_ht/(DBH_in+1)),2), Ht_ft)) %>%
  mutate(CBH_ft_calc = ifelse(is.na(CBH_ft) & Status != "D", round(exp(b1_cbh + b2_cbh/(DBH_in+1)),2), CBH_ft))


# examine figures
trees.calc %>%
  ggplot(aes(x=DBH_cm, y=Ht_ft_calc*0.3048, col=Predicted)) +
  facet_wrap(~Species) +
  geom_point() +
  geom_line(aes(y=(4.5+exp(b1_ht + b2_ht/(DBH_in+1)))*0.3048), color="red", size=1) +
  scale_color_manual(values=c("black","blue")) +
  theme_bw() 

trees.calc %>%
  ggplot(aes(x=DBH_cm, y=CBH_ft_calc*0.3048, col=Predicted)) +
  facet_wrap(~Species) +
  geom_point() +
  geom_line(aes(y=exp(b1_cbh + b2_cbh/(DBH_in+1))*0.3048), color="red",size=1) +
  scale_color_manual(values=c("black","blue")) +
  theme_bw() 

# export final fit data
write.csv(trees.calc, file="data/Field_plots_2019/Cleaned_data/trees_cleaned_with_heights_grouped.csv",row.names=FALSE)

# RULE 14: assume DF trees are live. could consider DF trees as live during time of LiDAR flight
# or dead to better reflect current fuel conditions. determine what to go with based
# on fuels regression fit
# this will be incorporated into fvs and fuelcalc prep below

####
# 7. saplings
####

saps.in = read.xlsx("data/Field_plots_2019/Raw_data/GRTE_LiDAR_field_data_2019.xlsx",
                    sheet="Saplings", colNames = TRUE)

head(saps.in)
summary(saps.in)  # should be all < 1.4 m height

saps.in %>%
  filter(Height_m_1>1.4)  # exclude this 1 live tree that should not have been counted

# just one with NA in there
saps.in %>%
  filter(is.na(Q)) # 1 with record that 0 live

saps.in %>%
  filter(is.na(Height_m_1))  # just dead ones

saps.in %>%
  filter(Species=="Unknown")  # just dead ones

# subset to just live saplings, drop nas, drop > 1.4 m height
saps.out = saps.in %>%
  filter(Live_count>0,
         Height_m_1 <=1.4)

summary(saps.out)
unique(saps.out$Species)

write.csv(saps.out, "data/Field_plots_2019/Cleaned_data/saplings_cleaned.csv", row.names=FALSE)

# also create sapling data for input to FVS
# plot size multiplier for saps, only on T1 and T3
sap_transect_size = ((2 * 26))/10000

# sapling density and mean height
saps.tph = saps.out %>%
  dplyr::select(-Dead_count) %>%  # exclude dead
  mutate(tph_adj = (1/sap_transect_size),
         tph_trans = Live_count * tph_adj) %>%
  group_by(Plot_code,Q,Species) %>%
  summarise(tph_saps=sum(tph_trans))

saps.tph %>%
  group_by(Species) %>%
  tally()

saps.ht = saps.out %>%
  dplyr::select(-c(Live_count,Dead_count)) %>%
  gather(key="Height_nbr",value="Height_m",-c(Plot_code:Species), na.rm=TRUE) %>%
  group_by(Plot_code,Q,Species) %>%
  summarise(mean_sap_ht_m=mean(Height_m)) %>%
  left_join(saps.tph, by=c("Plot_code","Q","Species")) %>%
  # also set POBA and POAN to POTR
  mutate(Species = ifelse(Species %in% c("POAN","POBA"),"POTR",Species))

summary(saps.ht)
saps.ht %>%
  group_by(Species) %>%
  tally()

####
# 8. Prep for canopy fuels calculations in FVS 
####

# this will include all live and dead trees, DF trees counted as live,
# for trees where CBH not recorded, use CBH calculated from regression fit to live, large trees

tree.plots = unique(trees.calc$Plot_code)
length(tree.plots)

# does this include all sapling plots?
!unique(saps.ht$Plot_code) %in% tree.plots # yes

# read in blank database
conn=DBI::dbConnect(RSQLite::SQLite(), dbname = paste0("processed_data/fuels_map_variables/fvs/BlankDatabase.db")) # connect to the db
dbListTables(conn)

fvs.groupin = dbReadTable(conn,"FVS_GroupAddFilesAndKeywords")
fvs.plotin = dbReadTable(conn,"FVS_PlotInit")
fvs.standin = dbReadTable(conn, "FVS_StandInit")
fvs.treesin = dbReadTable(conn, "FVS_TreeInit")  

dbDisconnect(conn)

# create stand list
fvs.stand = fvs.standin %>%
  # remove na row
  filter(!is.na(Stand_ID)) %>%
  # add row for each stand (plot)
  add_row(Stand_ID=tree.plots) %>%
  # match up with information from plots
  left_join(plots.in, by=c("Stand_ID" = "Plot_code")) %>%
  mutate(Variant = "TT", 
         Inv_Year = 2019,
         Groups = paste("All_Stands",Vegetation_type),
         # tetons nf
         Location = 416,
         Aspect = Aspect_deg,
         Slope = round(tan(Slope_deg*pi/180)*100,2),
         ElevFt = round(Postprocessed_elevation_m * 3.28,2),
         # these from donato example
         Basal_Area_Factor = 0,
         Inv_Plot_Size = 1,
         Brk_DBH = 999,
         Num_Plots = 1) %>%
  dplyr::select(-c(Photo_code:CanopyCover_T4))

# get size of tree plot
plot_size_m2 = pi * 13^2
plot_size_acres = plot_size_m2/10000 * 2.471
tree_exp_factor = round((1/plot_size_acres),2)
# exception!!! Berry Glade unburned, multiply this expansion factor x 4 (trees only)

# create species lookup
sp.lookup = data.frame(Species_code = c("ABLA","PIAL","PICO","PIEN","POTR","PSME"),
                       Species_FVS = c("AF","WB","LP","ES","AS","DF"))

# create tree list
fvs.treescbh = fvs.treesin %>%
  # remove na row
  filter(!is.na(Stand_ID)) %>%
  # add row for each stand (plot)
  add_row(Stand_ID=tree.plots) %>%
  # join with tree data
  right_join(trees.calc, by=c("Stand_ID" = "Plot_code")) %>%
  # join with species codes
  left_join(sp.lookup, by=c("Species.y" = "Species_code")) %>%
  rename(Species = Species.x) %>%
  # calculate crown ratio, use value set in field if there
  # otherwise calculate from ht and cbh only if recorded
  mutate(CrRatio_calc = as.numeric(ifelse(CR_field != "N",CR_field,round(100*(Ht_ft_calc - CBH_ft_calc)/Ht_ft_calc,0)))) %>%
  mutate(Plot_ID = 1,
         StandPlot_ID = paste(Stand_ID,Plot_ID,sep="_"),
         Species = Species_FVS,
         DBH=DBH_in, 
         Ht = Ht_ft_calc, 
         # temporary status column
         Status_placeholder = Status,
         # damage code and ht to top kill if applicable
         Damage1 = damage_code,
         HtTopK = round(Height_topkill_m/0.3048,2),
         # if crown ratio is below 10, set to 1 per FVS
         # which will set it to 5%
         # leave NAs as is
         CrRatio = ifelse(CrRatio_calc<10, 1,CrRatio_calc),
         # expansion factor is set for all plots except BerryGlade where only
         # 1/4 of plot was sampled
         Tree_Count =ifelse(Stand_ID != "Con_BerryGlade_Unburned3",tree_exp_factor,tree_exp_factor*4)) %>%
  dplyr::select(-c(Q:CrRatio_calc))

# then create corresponding list for sapling entries
fvs.sapstemp = fvs.treesin %>%
  # remove na row
  filter(!is.na(Stand_ID)) %>%
  # add row for each stand (plot)
  add_row(Stand_ID=tree.plots) %>%
  # join with tree data
  right_join(saps.ht, by=c("Stand_ID" = "Plot_code")) %>%
  # join with species codes
  left_join(sp.lookup, by=c("Species.y" = "Species_code")) %>%
  rename(Species = Species.x) %>%
  mutate(Plot_ID = 1,
         StandPlot_ID = paste(Stand_ID,Plot_ID,sep="_"),
         Species = Species_FVS,
         # trees below bh get assigned small, non-zero dbh
         DBH=0.1, 
         Ht = round(mean_sap_ht_m/0.3048,2), 
         # temporary status column
         Status_placeholder = "L",
         # tree count, convert to acres
         Tree_Count =round(tph_saps/2.471,2)) %>%
  dplyr::select(-c(Q:Species_FVS))

# combine these two
fvs.treesout = rbind(as.data.frame(fvs.treescbh),as.data.frame(fvs.sapstemp)) %>%
  # add History based on L & DF = live, D = dead
  mutate(History = ifelse(Status_placeholder %in% c("L","DF"),1,8)) %>%
  # sort by Stand_ID to integrate saps with trees
  arrange(Stand_ID) %>%
  # add tree number
  group_by(Stand_ID) %>%
  mutate(Tree_ID = row_number()) %>%
  dplyr::select(-Status_placeholder)

summary(fvs.treesout)

fvs.plotin = fvs.plotin %>%
  filter(!is.na(Stand_ID))

# write out db
# note need to manually copy to FVS_Group section
fvs.out = dbConnect(RSQLite::SQLite(), "processed_data/fuels_map_variables/fvs/trees_all_dflive_cbhcalc.db")
dbWriteTable(fvs.out,"FVS_GroupAddFilesAndKeywords",fvs.groupin)  # needs to be manually updated
dbWriteTable(fvs.out, "FVS_PlotInit", fvs.plotin)
dbWriteTable(fvs.out, "FVS_StandInit", fvs.stand)
dbWriteTable(fvs.out, "FVS_TreeInit", fvs.treesout)
dbListTables(fvs.out)
dbDisconnect(fvs.out)
