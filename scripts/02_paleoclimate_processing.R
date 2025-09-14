# ==============================================================================
# Paleoclimatic Data Processing for Bioclimatic Analysis
# ==============================================================================
# Description: Download and process paleoclimatic data from WorldClim/PaleoClim
#              for multiple time periods and prepare geodiversity layers
# 
# Input:  Spatial polygon (satz_sf) defining region of interest
# Output: Processed bioclimatic rasters in data/processed/
#         - bio_{1,4,12,15}_{period}.tif for each time period
#         - geodiv_{variable}.tif for geodiversity layers
#
# Time periods processed:
# - Present (current)
# - Pliocene (M2, ~3.3 Ma)
# - Last Glacial Maximum (LGM, ~21 ka)
# - Last Interglacial (LIG, ~130 ka) 
# - Late Holocene (LH, 4.2-0.3 ka)
#
# Bioclimatic variables extracted:
# - bio_1: Annual Mean Temperature
# - bio_4: Temperature Seasonality (Standard Deviation)
# - bio_12: Annual Precipitation
# - bio_15: Precipitation Seasonality (Coefficient of Variation)
# =============================================================================

## 1. Load required packages ----

# Prepare environmental layers
library(tidyverse)
library(stars)
library(terra)
library(rpaleoclim)

satz_sf<- st_transform(satz_sf, crs = 4326)

raster_prep <- function(x, y){
  filename<- word(x, 6,sep = '/')
  
  read_stars(x) %>%
    st_warp(y) %>%
    st_crop(y)%>% as(., "SpatRaster") %>%
    writeRaster(paste0('data/processed/',filename),
                overwrite=TRUE)}


process_bioclim <- function(raster_obj, period_name) {
  # Seleccionar bioclim variables de interés
  selected_vars <- raster_obj[[names(raster_obj) %in% 
                                 c('bio_1', 'bio_4', 'bio_12', 'bio_15')]]
  
  # Recortar y enmascarar usando la región
  processed <- crop(selected_vars, vect(satz_sf), mask = TRUE)
  
  # Escribir los archivos a disco
  writeRaster(processed, 
              filename = paste0('data/processed/',  
                                names(processed),'_',
                                period_name,'.tif'),
              overwrite = TRUE)}

# Current
present <- paleoclim(period = 'cur',region = satz_sf,
                     resolution = '2_5m')

process_bioclim(present, 'present')


# Pliocene: M2
pliocene <- paleoclim(period = 'm2',region = satz_sf, skip_cache = T,
                     resolution = '2_5m')

process_bioclim(pliocene, 'pliocene')

# Last Glacial Maximum
lgm <- paleoclim(period = 'lgm',region = satz_sf, skip_cache = T,
                      resolution = '2_5m')

process_bioclim(lgm, 'lgm') 

# Last interglacial
lig <- paleoclim(period = 'lig',region = satz_sf,
                      resolution = '2_5m')
  
process_bioclim(lig, 'lig')
  
# Late holocene 
lh <- paleoclim(period = 'lh',region = satz_sf,
                   resolution = '2_5m')
  
process_bioclim(lh, 'lh')



# Geodiversity
y <- read_stars('data/processed/bio_1_present.tif')

geodiv <- list.files('data/raw', pattern = '.tif$', full.names = T) %>%
  lapply(., read_stars) %>%
  lapply(., st_warp, read_stars('data/processed/bio_1_present.tif')) %>% 
  lapply(., as, "SpatRaster") %>% rast() %>% 
  crop(vect(satz_sf), mask = T) 

writeRaster(geodiv, paste0('data/processed/geodiv_', names(geodiv), '.tif'),
            overwrite=TRUE)
  

