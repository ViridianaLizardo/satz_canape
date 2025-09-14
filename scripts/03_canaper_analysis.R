# ==============================================================================
# SUPPLEMENTARY MATERIAL SCRIPT
# CANAPE (Categorical Analysis of Neo- and Paleo-Endemism) Workflow
#
# Description: Complete workflow for phylogenetic endemism analysis using CANAPE
#              methodology to identify neoendemic, paleoendemic, and mixed regions
#
# Inputs:
#   - data/processed/occurrences_final.csv: Species occurrence data
#   - data/processed/final_tree.nex: Phylogenetic tree in Nexus format
#
# Outputs:
#   - canaper_sf: Spatial dataframe with endemism classification results
#   - rand_res: Randomization test results for phylogenetic diversity metrics
#
# Analysis Steps:
# 1. Spatial data preparation and coordinate system management
# 2. Hexagonal grid creation (100km resolution) for community matrix
# 3. Community matrix construction from occurrence data
# 4. Phylogenetic tree processing and name validation
# 5. CANAPE randomization testing with 9999 replicates
# 6. Endemism type classification and significance testing
# ==============================================================================

# 1. LOADING LIBRARIES----


# Core spatial and data manipulation packages
library(tidyverse)
library(sf)
library(terra)

# Phylogenetic analysis packages
library(ape)        # For handling phylogenies

# CANAPE-specific packages
library(canaper)    # For CANAPE analysis



# 2. DATA PREPARATION----


# Define coordinate reference system (SAD_1969_Albers_South_America)
crs_102033 <- "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs +type=crs"

# Load base spatial data
data("World")

# Define South American Transition Zone (SATZ) as all of South America excluding Falklands
satz_sf <- filter(World, continent == 'South America' & iso_a3 !='FLK') %>% 
  tally() %>% 
  st_transform(., crs = crs_102033) %>% 
  mutate(Region = 'Sudamerica')


# Load occurrence data
occ_sf <- read.csv('data/processed/occurrences_final.csv') %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>% 
  filter(!is.na(label)) %>% 
  st_transform(., crs = crs_102033) %>% 
  st_join(satz_sf) %>% 
  filter(!is.na(Region))

# Load phylogenetic tree
tree <- read.nexus('data/processed/final_tree.nex')


# 3. VISUALIZATION SETUP----

# Create base map template for visualizations
gg_basemap <- ggplot(outline) + 
  geom_sf(fill = 'grey70', col = NA, alpha = 0.4) +
  theme_light() + 
  theme(legend.position = "top",
        legend.title.position = 'top',
        legend.key.width = unit(35, "points"),
        legend.key.height = unit(7, "points"))

# Preview base map
gg_basemap + 
  coord_sf(xlim = c(-2500000 ,2600000),
           ylim = c(-2500000 , 4700000))


# 4. COMMUNITY MATRIX CONSTRUCTION----


# Define hexagonal grid cell size (100km)
cell_size <- 100000

# Create hexagonal grid over study area
hex_grid <- st_make_grid(x = satz_sf,
                         crs = satz_sf,
                         cellsize = cell_size,
                         square = F) %>% 
  st_as_sf() %>%
  rownames_to_column(var = 'cell') %>%
  st_join(satz_sf[,1]) %>% 
  na.omit()

# Construct community matrix from occurrence data
comm_sf <- hex_grid %>%            
  st_join(occ_sf, .) %>%      
  st_drop_geometry() %>%         
  group_by(cell, label) %>%
  summarise(records = n()) %>%     
  group_by(cell) %>%
  mutate(rich = n_distinct(label),
         redundancy = 1- (rich/sum(records))) %>% 
  pivot_wider(id_cols = c(cell,rich, redundancy),
              names_from = label,
              values_from = records,values_fill = 0) %>%
  inner_join(hex_grid, ., by = 'cell') %>% 
  select(-n)

# Create non-spatial version for analysis
comm_df <- st_drop_geometry(comm_sf)


# 5. PHYLOGENETIC DATA PREPARATION----

# Prepare community matrix for phylogenetic analysis
comm <- comm_df %>% 
  select(-rich, -redundancy) %>%
  column_to_rownames('cell')

# Check species names between tree and community matrix
names_comm <- data.frame(sp = colnames(comm),
                         data = rep(1, ncol(comm))) %>%
  pull(data,name = sp)

geiger::name.check(tree, names_comm)

# Check for species without locations
sites_per_sp <- colSums(comm)
sites_per_sp[sites_per_sp == 0]


# 7. CANAPE ANALYSIS PREPARATION----


# Prepare binary community matrix
com_matrix <- comm_df %>% 
  select(-rich, -redundancy) %>%
  column_to_rownames('cell') %>%           
  mutate_all(~ ifelse(. > 0, 1, 0))


# 8. RANDOMIZATION TEST (CANAPE) ----

# Set random seed for reproducibility
set.seed(12345)


# Run CANAPE randomization test
rand_res <- cpr_rand_test(
  com_matrix, tree,
  null_model = "curveball",
  n_reps = 9999, 
  n_iterations = 8000,
  tbl_out = TRUE
)



# 9. RESULTS PROCESSING

# Define significance levels
clasif_levels <- names(cpr_signif_cols)

# Classify endemism types and significance
canaper_sf <- cpr_classify_endem(rand_res) %>%
  cpr_classify_signif("pd") %>%
  cpr_classify_signif("rpd") %>%
  cpr_classify_signif("pe") %>%
  cpr_classify_signif("rpe") %>%
  rename(., cell = site) %>%
  left_join(comm_sf[,1:3], .) %>%
  mutate(
    endem_type = factor(endem_type, 
                        levels = c('neo', 'paleo', 'not significant', 'mixed', 'super')),
    pd_signif = factor(pd_signif, levels = clasif_levels),
    rpd_signif = factor(rpd_signif, levels = clasif_levels),
    pe_signif = factor(pe_signif, levels = clasif_levels),
    rpe_signif = factor(rpe_signif, levels = clasif_levels)
  )

