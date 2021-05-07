# Local climate impacts of tropical deforestation

# Scripts for processing and analysing satellite data
Author: Callum Smith (ee13c2s@leeds.ac.uk)

The scripts for processing the raw data used in the following paper:

Smith, C., Baker, J.C.A. & Spracklen, D.V. In review. Transition from shifting to commodity agriculture doubles warming from tropical forest loss (AGU Advances).

Analysis in this study was carried out using satellite data from a range of sources. We employed tree canopy cover data from the Global Forest Change group (Hansen et al., 2013) and a range of climate and land surface variables from MODIS, as well as Preicpitation data from CHIRPS and elevation data from GLOBE.


# List of scripts

## Scripts used to process raw satellite data to harmonised netcdf files:

1. LAI_processing.py
2. ET_processing.py - this includes putting into folders, etc as below.
3. etc.


## Scripts with additional functions called during harmonisation and data processing

1. cube_funcs.py
2. harmonise.py


# categorise the forest loss hansen data

1. canopy_cover.py


## Scripts to analyse the satellite data

1. Figure_3_all_vars.py
2. Figure_2.py




Dataset stored here: "Leeds data repo". Describe the dataset. Full details of all source datasets are provided in 'ETransition from shifting to commodity agriculture doubles warming from tropical forest loss' in review in AGU Advances.











# JESS
## List of scripts
1. Scripts to process raw satellite, ERA5 and CMIP model data to harmonised netcdfs
2. chirps_p_data_processing.py
3. clara-a1_radiation_data_processing.py
gleam_et_data_processing.py
grace_tws_data_processing.py
p-lsh_et_data_processing.py
era5_et_data_processing.py
era5_lai_high_data_processing.py
era5_pr_data_processing.py
era5_rdn_data_processing.py
cmip5_et_data_processing.py
cmip6_et_data_processing.py
## Scripts with additional functions called during harmonisation and data processing
cube_funcs.py
harmonise.py
## Scripts for processing 500-m sinusoidal MODIS ET tiles
modis_et_make_dirs.sh # this bash script sorts modis raw data into directories by year and day of year
modis_et_hdf2tiff_merge.sh # this bash script calls modis_et_merge (must change path to this in modis_et_hdf2tiff_merge.py), modis_et_merge in turn calls modis_et_update_tif.py so this must be in the same directory
modis_et_final_processing.py # this python script loops over files from each year, and for each 8-day file adds date, converts units, converts to monthly data
harmonise_modis_et.py # this script harmonises MODIS ET data to single data cube of desired resolution
## Script to process raw river runoff data from ANA
amazon_river_station_data_processing.py
## Script to process raw flux tower data from LBA
lba_eco_data_processing.py
## Scripts to extract basin monthly-mean values of ET, leaf area index, precipitation and radiation
get_catchment_balance_ET_all_basins.py # this script calculates catchment-balance ET for Amazon and sub-basins
get_satellite_et_all_basins_interannual.py
get_satellite_lai_all_basins_interannual.py
get_satellite_pre_all_basins_interannual.py
get_satellite_rdn_all_basins_interannual.py
get_era5_et_all_basins_interannual.py
get_era5_lai_high_all_basins_interannual.py
get_era5_pr_all_basins_interannual.py
get_era5_rdn_all_basins_interannual.py
get_cmip5_et_all_basins_interannual.py
get_cmip5_lai_all_basins_interannual.py
get_cmip5_pr_all_basins_interannual.py
get_cmip5_rdn_all_basins_interannual.py
get_cmip6_et_all_basins_interannual.py
get_cmip6_lai_all_basins_interannual.py
get_cmip6_pr_all_basins_interannual.py
get_cmip6_rdn_all_basins_interannual.py
## Script to collate Amazon monthly mean ET, P, RDN and LAI values
Dataset stored here: https://zenodo.org/record/4271331#.YD-emS2l30o Amazon basin-mean monthly estimates of evapotranspiration estimated from catchment-balance analysis, satellites (MODIS, P-LSH and GLEAM), reanalysis (ERA5) and the CMIP5 and CMIP6 climate models. Catchment level estimates of climate variables that influence ET are also included (precipitation, radiation and leaf area index). Full details of all source datasets are provided in 'Evapotranspiration in the Amazon: spatial patterns, seasonality and recent trends in observations, reanalysis and CMIP models' in review in Hydrology and Earth System Sciences.

all_data_catchment_level.py
## Script to estimate uncertainty in Amazon catchment-balance ET
quantify_uncertainty_in_catchment_et.py
