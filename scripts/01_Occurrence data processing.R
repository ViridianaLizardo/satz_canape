# ==============================================================================
# SUPPLEMENTARY MATERIAL: GBIF Occurrence Processing & Native Species Filtering
# 
# Description: Data preparation pipeline for CANAPE analysis including:
#              - GBIF occurrence data retrieval with error handling
#              - Native species identification using Plants of the World Online
#              - Phylogenetic tree pruning to native species only
#
# Inputs:
#   - GBIF species keys (from taxonomic backbone)
#   - Study region polygon (WKT format)
#   - Full phylogenetic tree (GBOTB.extended.TPL)
#   - TDWG Level 3 geographic regions
#
# Outputs:
#   - data/processed/occurrences_final.csv: Filtered occurrence data (native species only)
#   - data/processed/final_tree.nex: Pruned phylogenetic tree (native species only)
#   - gbif_ping.csv: Raw GBIF occurrence counts per species
#
#==============================================================================


# 1. LOADING LIBRARIES----

library(dplyr)    # For data manipulation
library(purrr)    # For functional programming
library(rgbif)    # For GBIF data access
library(kewr)     # For POWO data access
library(sf)       # For spatial operations
library(ape)      # For phylogenetic tree manipulation


# 2. GBIF OCCURRENCE DATA DOWNLOAD FUNCTION----


#' Safe GBIF Occurrence Data Retrieval
#'
#' @param species_keys Vector of GBIF species keys
#' @param polygon_wkt WKT polygon defining search area
#' @param output_file Name of output CSV file
#' @param batch_size Number of species to process at once
#' @param max_attempts Maximum API query attempts per species
#' @param delay Initial delay between API calls (in seconds)
#'
#' @return Data frame with occurrence counts per species

occ_explore_safe <- function(species_keys, 
                             polygon_wkt, 
                             output_file = "gbif_results.csv",
                             batch_size = 100,
                             max_attempts = 3,
                             delay = 1) {
  
  # Initialize results file if it doesn't exist
  if (!file.exists(output_file)) {
    data.frame(speciesKey = character(),
               occ = integer(),
               stringsAsFactors = FALSE) %>%
      write.csv(output_file, row.names = FALSE)
  }
  
  # Internal function with error handling for single species query
  get_occ_data <- function(x) {
    attempt <- 1
    result <- NULL
    
    while (attempt <= max_attempts && is.null(result)) {
      try({
        res_occ <- occ_search(
          taxonKey = x,
          geometry = polygon_wkt,
          hasCoordinate = TRUE,
          fields = "minimal",
          limit = 300
        )$data
        
        if (is.null(res_occ)) {
          result <- data.frame(speciesKey = x, occ = 0)
        } else {
          result <- res_occ %>%
            reframe(occ = n()) %>%
            mutate(speciesKey = x, .before = occ)
        }
      }, silent = TRUE)
      
      if (is.null(result)) {
        message("\nError with speciesKey: ", x, " (Attempt ", attempt, "/", max_attempts, ")")
        Sys.sleep(delay * attempt)
        attempt <- attempt + 1
      }
    }
    
    if (is.null(result)) {
      message("\nFailed to query speciesKey: ", x, " after ", max_attempts, " attempts")
      result <- data.frame(speciesKey = x, occ = NA)
    }
    
    return(result)
  }
  
  # Process species in batches with progress tracking
  total_species <- length(species_keys)
  processed <- 0
  
  # Check which species have already been processed
  if (file.exists(output_file)) {
    completed <- read.csv(output_file)$speciesKey
    species_keys <- setdiff(species_keys, completed)
  }
  
  for (i in seq(1, length(species_keys), by = batch_size)) {
    batch <- species_keys[i:min(i + batch_size - 1, length(species_keys))]
    
    batch_results <- map(batch, function(x) {
      res <- get_occ_data(x)
      cat(".")  # Progress indicator
      processed <<- processed + 1
      if (processed %% 50 == 0) {
        cat(" ", round(processed/total_species*100, 1), "%\n")
      }
      return(res)
    }) %>% bind_rows()
    
    # Append batch results to output file
    write.table(batch_results, 
                file = output_file,
                append = TRUE,
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
    
    # Delay between batches to avoid API overload
    Sys.sleep(delay * 2)
  }
  
  # Process and return final results
  final_results <- read.csv(output_file) %>%
    left_join(gbif_keys_all, by = "speciesKey") %>%
    filter(occ > 0)
  
  return(final_results)
}

# 3. DATA PROCESSING PIPELINE ----

# Example usage:
res_gbif <- occ_explore_safe(
   species_keys = gbif_keys_all_clean$speciesKey,
   polygon_wkt = wkt_polygon,
   output_file = "gbif_ping.csv")

# Process occurrence coordinates
gbif_result_coords <- cbind(gbif_result, 
                            st_coordinates(gbif_result))

# Clean and filter occurrence data
gbif_result2 <- gbif_result_coords %>% 
  select(-X, -Y) %>%
  st_drop_geometry() %>%
  left_join(., res_gbif) %>% 
  select(class, family, label, lat, lon, perc_zt) %>%
  filter(class %in% c('Magnoliopsida', 'Liliopsida')) %>%
  distinct()  


# 4. PHYLOGENETIC TREE PROCESSING----


# Prune tree to include only species with occurrence data
to_prune <- tree_all$tip.label[!tree_all$tip.label %in% labs_in_db]
final_tree <- drop.tip(tree_all, to_prune)

# Save pruned tree
write.nexus(final_tree, file = 'data/processed/final_tree.nex')

# Filter occurrences to match pruned tree
gbif_result3 <- gbif_result2 %>%
  filter(label %in% final_tree$tip.label) 


# 5. NATIVE SPECIES IDENTIFICATION----

# Match species with Plants of the World Online (POWO)
POWO <- unique(gbif_result3$label) %>%
  kewr::match_knms(.) %>% 
  tidy()

# Load TDWG Level 3 polygons for study region
tdwg_lev3_satz <- st_read('data/wgsrpd-master/level3') %>%
  st_make_valid() %>%
  st_intersection(., st_read('data/SATZ.gpkg')) %>% 
  pull(LEVEL3_COD)

#' Get Native Distribution Status from POWO
#'
#' @param x Species name (IPNI ID)
#' @return Data frame with native status

get_distribution <- function(x) {
  cat(paste0('starting: ', x))  # Print current species name
  
  # Try to retrieve distribution data with error handling
  res_powo_prev <- tryCatch(
    {
      lookup_powo(x, distribution = TRUE) %>%
        tidy()
    },
    error = function(e) {
      return(NULL)
    }
  )
  
  # Determine distribution status
  if (is.null(res_powo_prev) || !'distribution' %in% names(res_powo_prev)) {
    res <- 'unknown'
  } else {
    res_powo <- res_powo_prev %>%
      select(distribution) %>%
      unnest(cols = distribution)
    
    # Check for introduced distributions
    if ('introduced' %in% names(res_powo)) {
      intro <- res_powo %>%
        select(introduced) %>%
        unnest(cols = introduced) %>%
        select(tdwgCode) %>%
        mutate(presence = 'introduced')
      
      native <- res_powo %>%
        select(natives) %>%
        unnest(cols = natives) %>%
        select(tdwgCode) %>%
        mutate(presence = 'native')
      
      distribution <- bind_rows(native, intro) %>%
        rename(LEVEL3_COD = tdwgCode)
      
      # Check native status in study region
      keep <- distribution %>%
        filter(LEVEL3_COD %in% tdwg_lev3_satz) %>%
        group_by(presence) %>%
        tally() %>%
        filter(presence == 'native')
      
      res <- if (nrow(keep) == 1) 'native' else 'introduced'
    } else {
      res <- 'native'
    }
  }
  
  cat(paste0(': ', res, ' .... done!\n'))
  return(data.frame(ipni_id = x, distribution = res))
}

# Wrapper function with error handling for batch processing
try_get_distr <- function(x) {
  tryCatch(get_distribution(x),
           error = function(e) {
             return(data.frame(ipni_id = x, distribution = 'unknown'))
           })
}

# Query native status for all species
native_test <- lapply(POWO$ipni_id, try_get_distr)

# Identify native species to retain
to_keep <- native_test %>% 
  bind_rows() %>% 
  na.omit() %>%
  left_join(POWO, .) %>%
  rename(label = submitted) %>%
  filter(distribution != 'introduced') %>%
  pull(label) %>% unique()


# 6. FINAL DATA PRODUCTS----

# Load complete phylogenetic tree
tree_all <- V.PhyloMaker2::GBOTB.extended.TPL

# Prune tree to native species only
to_prune <- tree_all$tip.label[!tree_all$tip.label %in% to_keep]
final_tree <- drop.tip(tree_all, to_prune)

# Save final tree
write.nexus(final_tree, file = 'data/processed/final_tree.nex')

# Save final occurrence data (native species only)
gbif_result3 %>% 
  filter(label %in% to_keep) %>% 
  write.csv('data/processed/occurrences_final.csv', row.names = TRUE)