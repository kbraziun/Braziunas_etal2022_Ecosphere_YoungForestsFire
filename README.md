# readme for Braziunas et al. in review Young Forests and Fire manuscript

## see also

This github repository includes data and code associated with Braziunas et al. In review, Young forests
and fire: Using lidar-imagery fusion to explore fuels and burn severity in a subalpine forest reburn.
Data and code are provided here for peer review, and a permanent data deposit will be uploaded to
the Environmental Data Initiative (EDI) upon publication.

## purpose

This readme gives an overview of directory structure, files, and steps for recreating
outputs and analyses associated with Braziunas et al. in review.

**Manuscript citation:** Braziunas, K.H., D.C. Abendroth, and M.G. Turner. In review. Young
forests and fire: Using lidar-imagery fusion to explore fuels and burn severity in a subalpine
forest reburn

## organization and file descriptions

Directory structure and files:

- `analysis/`: Results data.
  - `fuels_prediction_map/`: Summary statistics for field data, fuels models, and final fuels map.
      Includes field data plots and fuels summary statistics (`plot_summary_stats.csv` and
      `fuels_summary_stats.csv`), final linear models for predicting forest and shrubland fuels
      (`model_final_selection.csv` and `model_final_preditors.csv`), predicted versus observed
      values for final model fits (`forest_predicted_observed.csv` and `shrubland_predicted_observed.csv`),
      and comparison of final lidar-imagery fusion fuels map and LANDFIRE with field plot data
      (`fuels_map_landfire_comparison.csv`).
  - `q1_young_old_forest_comparison`: Data (`young_old_fuels_severity.csv`) and summary statistics
      (`young_old_forest_medians.csv`) used to answer question 1, How do pre-fire fuels and burn
      severity copmare between young (~30-year-old) and mature (> 125-year-old) forests that burned
      under similar fire weather conditions?
  - `q2_fuels_fire_severity_models`: Data (`weather_fuels_severity.csv`) used to answer question 2,
      How well do pre-fire fuels and forest structure predict burn severity under extreme versus
      moderate fire weather?
- `data/`: Input field and weather data.
  - `Field_plots_2019`: Field data collected for this study in 2019 in Grand Teton National Park. Includes
      `Raw_data` folder with `GRTE_LiDAR_field_data_2019.xls` that includes associated metadata on sheet 1.
      Also includes `Cleaned_data` folder, with `.csv` outputs from R data cleaning script.
  - `GRTE_RAWS`: RAWS weather station data during the Berry Fire.
- `processed_data`: Data created during analyses.
  - `Berry_Fire_sample_polys`: Shapefiles and associated lidar and imagery (NAIP) predictors for random
      samples of grid cells from the Berry Fire used to answer Q1 (files with `age` in title) and Q2
      (files with `wx` in title).
  - `field_plot_selection`: Shapefiles for field plot footprints used in extracting lidar and imagery
      predictors for create linear regression models to predict fuels.
  - `fuels_map_variables`: Fuels used as dependence variables in linear regression models
      (`field_plot_2019_fuels.csv`). Lidar (`lidar_metrics.csv`)
      and imagery (`naip_metrics.csv`) predictors for field plot footprints used in linear regression models.
      Also includes `fvs` folder with inputs and outputs from estimating fuels in Forest Vegetation
      Simulator software.
  - `GRTE_rasters`: Final fuels map rasters (in `final_fuels_map` folder) and rasters used in creating
      final fuels maps or performing manuscript analyses. These include a digital elevation model, lidar
      and imagery predictors at 30 m resolution, and vegetation map at 30 m resolution. Also includes a
      `fuels_map_comparison` folder with unmasked version of canopy fuels for LANDFIRE comparison. One
      raster `grte_30m_naip_masked.tif` is not uploaded because it exceeds Github's file size limits,
      but this can be created from data included here by running `step07_create_fuels_map.R`.
  - `GRTE_shps`: Outline of study region in which fuels were mapped.
- `scripts`: R scripts used to recreate data and analyses. Scripts are numbered in order they were run.
  Scripts with `.R` can be rerun with data provided in deposit. Scripts with `.Rmd` require additional
  data not included in deposit such as NAIP imagery files, raw lidar point cloud, LANDFIRE fuels
  maps, MTBS burn severity rasters, and fire progression maps. Examples of some outputs of these scripts
  and included in knitted `.html` files.
  - `step01_fuels_plots_data_cleaning_prep.R`: Quality check and cleaning of field plot data.
  - `step02_field_plot_fuels_calculations.R`: Calculate fuel loads for field plots.
  - `step03_berry_fire_point_selection.Rmd`: Random selection of points for Q1 and Q2 analyses. Also
      calculates summary information about the Berry Fire on high and low spread days.
  - `step04_las_statistics_extraction.Rmd`: Extract lidar metrics from lidar point cloud for fitting
      fuels regressions, predicting fuel loads in Berry Fire sample points (Q1 and Q2), and creating
      the final fuels map. This script also creates the common grid and masks used for final fuels map.
  - `step05_naip_statistics_extraction.Rmd`: Extract  imagery metrics from NAIP data for fitting fuels
      regressions, predicting fuel loads in Berry Fire sample points (Q1 and Q2), and creating
      the final fuels map.
  - `step06_fit_fuels_regressions.R`: Fit linear regression models to predict field-measured fuels from
      lidar and imagery predictors.
  - `step07_create_fuels_map.R`: Create lidar-imagery fusion fuels maps.
  - `step08_fuels_fit_analyses.R`: Use final linear regression models to predict fuels in field data
      plots and generate summary data on model fits (predicted versus observed values).
  - `step09_fuels_map_landfire_comparison.Rmd`: Compare predicted fuels from final lidar-imagery fusion fuels
      and LANDFIRE fuels maps (rasters) with field data observations (using centroids from field plots
      to extract corresponding raster values).
  - `step10_fuels_severity_analyses.R`: Analyses for Q1 (fuels and burn severity in young versus mature
      forests) and Q2 (how well fuels predict burn severity under extreme versus moderate fire weather).

## platforms

- Operating systems and software used for development and implementation
  - OS: Windows 10 Pro
  - R version: 3.6.1
  - Forest Vegetation Simulator Complete Package Software Version: 2020.03.11
  - ArcGIS Desktop 10.6
  