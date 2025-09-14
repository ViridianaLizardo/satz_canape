# Reproducible repository for: [Add title here]

This repo has the files and scripts to recreate phylogenetic diversity calculations, CANAPE, and statistical analysis. 

## Data
- env_vars: raster files of the enviromental variables analized
- processed: comunity matrix, phylogenetic tree and, occurrence data.

## Scripts: With the key parts for recreate this study.

## Results: metrics_final_100.gpkg
This file is a **GeoPackage** containing the results of the CANAPE (Categorical Analysis of Neo- and Paleo-Endemism) analysis at a spatial resolution of 100 km grid cells. Each record corresponds to a grid cell and includes taxonomic, phylogenetic, and environmental metrics used in the study.  

#### Diversity metrics
- `rich`: Species richness per grid cell.  
- `redundancy`: Proportion of species shared with other cells.  
- `pd_obs`, `pd_rand_mean`, `pd_rand_sd`, `pd_obs_z`, … :  
  Observed PD, null expectations, and standardized effect sizes.  
- `pd_alt_*`: Alternative PD metric values.  
- `rpd_*`: Relative PD metrics.  
- `pe_obs`, `pe_rand_mean`, `pe_obs_z`, … :  
  Observed PE, null model expectations, and significance tests.  
- `pe_alt_*`: Alternative PE metric values.  
- `rpe_*`: Relative PE metrics.  

#### CANAPE results
- `endem_type`: Endemism category (e.g., paleo-, neo-, mixed, not significant).  
- `pd_signif`, `rpd_signif`, `pe_signif`, `rpe_signif`: Significance results for each metric.  

#### Phylogenetic ages
- `mean_sp_age`: Mean species age in the cell.  
- `mean_crown_age_fam`: Mean crown age of families present.  
- `fam_age`: Age of dominant families.  

#### Biogeographic context
- `Subregion`, `Provincias`, `perc_area`, `NAME`, `FROMAGE`:  
  Regional and provincial classifications, and coverage of the cell in each unit.  

#### Paleoclimate and current climate variables (WorldClim)
Variables are provided for five periods:  
- `lgm` = Last Glacial Maximum  
- `lh` = Late Holocene  
- `lig` = Last Interglacial  
- `pliocene` = Mid-Pliocene Warm Period  
- `present` = Current climate  

Variables:  
- `bio_1`: Mean Annual Temperature  
- `bio_4`: Temperature Seasonality  
- `bio_12`: Annual Precipitation  
- `bio_15`: Precipitation Seasonality
