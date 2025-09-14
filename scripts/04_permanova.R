# SUPPLEMENTARY MATERIAL SCRIPT
# PERMANOVA and Kernels 



# This script performs a Categorical Analysis of Neo- and Paleo-Endemism (CANAPE)
# across the South American Transition Zone (SATZ) and correlates these
# classifications with environmental variables.




# LOAD REQUIRED LIBRARIES ----
library(MASS)       # Modern applied statistics functions
library(tidyverse)  # Data manipulation and transformation
library(tidyr)      # Data tidying and reshaping
library(sf)         # Simple features for spatial data handling
library(tmap)       # Thematic mapping and spatial visualization
library(vegan)      # Community ecology and multivariate analysis
library(ggplot2)    # Grammar of graphics for data visualization
library(ggrepel)    # Prevent overlapping labels in ggplot2



library(patchwork)


# READ SPATIAL DATA----
metrics_sf <- st_read("results/metrics_final_100.gpkg")%>%  
st_drop_geometry() 
  

# FILTER FOR ENV SPACE ANALYSIS----


df_present <- metrics_sf %>%
                            # Convert sf to dataframe
  dplyr::select(
    cell,
    endem_type,                                     # Endemic type classification
    bio_4_present,                                  # Temperature seasonality (standard deviation)
    bio_15_present,                                 # Precipitation seasonality (coefficient of variation)
    bio_1_present,                                  # Annual mean temperature
    bio_12_present                                  # Annual precipitation
  ) %>%
  filter(!is.na(endem_type)) %>%                    # Remove rows with missing endemic type
  filter(complete.cases(.))                         # Remove rows with any NA values


df_all<-metrics_sf %>%
  dplyr::select(
    cell,
    endem_type,                                     # Endemic type classification
    bio_1_lgm:bio_4_present) %>%
  pivot_longer(cols = bio_1_lgm:bio_4_present, 
               names_sep = '_', names_to = c('bio', 'n', 'time')) %>%
  mutate(bio = paste(bio, n, sep = '_'),
         var_type = if_else(bio %in% c('bio_4', 'bio_15'),
                            'Seasonality' ,'Annual Conditions'),
         var = if_else(bio %in% c('bio_1', 'bio_4'),
                            'Temperature' ,'Precipitation')) %>%
  
  filter(!is.na(endem_type)) %>%                    # Remove rows with missing endemic type
  filter(complete.cases(.)) %>%                         # Remove rows with any NA values
  dplyr::select(-n) %>%

  pivot_wider(id_cols = c(cell, endem_type, time, var_type),
              names_from = var, values_from = value) %>%
  group_by(var_type, time) %>% group_split()

# GLOBAL PERMANOVA + DISPERSION ANALYSIS ----

# Single function to handle both global and pairwise PERMANOVA
run_complete_analysis <- function(data, env_vars, group_var, analysis_name) {
  
  # Extract data
  env_data <- data[, env_vars]
  groups <- data[[group_var]]
  
  # Global PERMANOVA
  global <- adonis2(env_data ~ groups, method = "mahalanobis", permutations = 999)
  
  # Dispersion test
  dispersion <- permutest(betadisper(dist(env_data), groups))
  
  # Pairwise comparisons
  group_levels <- unique(groups)
  pairwise <- combn(group_levels, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      subset_data <- data[data[[group_var]] %in% pair, ]
      ad <- adonis2(subset_data[, env_vars] ~ subset_data[[group_var]],
                    method = "mahalanobis", permutations = 999)
      data.frame(
        Group1 = pair[1],
        Group2 = pair[2],
        F_value = ad$F[1],
        R2 = ad$R2[1],
        p_value = ad$`Pr(>F)`[1]
      )
    })
  
  # Return everything in a list
  list(
    analysis_name = analysis_name,
    global = global,
    dispersion = dispersion,
    pairwise = pairwise
  )
}


# EXECUTE ALL ANALYSES IN 2 LINES ----


set.seed(123)

# Run both analyses
results_seasonality <- run_complete_analysis(df_present, 
                                             c("bio_4_present", "bio_15_present"), 
                                             "endem_type", 
                                             "Seasonality")

results_annual <- run_complete_analysis(df_present,
                                        c("bio_1_present", "bio_12_present"),
                                        "endem_type", 
                                        "Annual")


# VIEW RESULTS
# seasonality
results_seasonality$global
results_seasonality$pairwise

# annual
results_annual$global
results_annual$pairwise


# KDE MODE CENTROIDS FOR PLOTTING ----
get_kde_peak2 <- function(x) {
  y <- na.omit(x)
  d <- kde2d(y$Temperature, y$Precipitation, n = 200)
  
  max_idx <- which(d$z == max(d$z), arr.ind = TRUE)
  
  data.frame(
    endem_type = unique(y$endem_type),
    mode_temp = d$x[max_idx[1, 1]],     
    mode_precip = d$y[max_idx[1, 2]],     
    var_type = unique(y$var_type),
    time = unique(y$time),
    n_points = nrow(y),                 
    stringsAsFactors = FALSE
  )
}

# Aplicar con manejo de errores
centroids_mode <- df_all %>%
  
  map(safely(get_kde_peak2)) %>%    # error handling
  map("result") %>%                # extract results
  bind_rows()

centroids_mode <-bind_rows(df_all) %>% group_by( time,var_type, endem_type) %>%
  group_split() %>% lapply(., get_kde_peak2) %>% bind_rows()%>%
  mutate(time = factor(time, 
                       levels = c('present', 'lh', 'lgm', 'lig', 'pliocene')),
         var_type = factor(var_type, 
                           levels = c('Annual Conditions', 'Seasonality')))



# SEGMENTS & SIGNIFICANCE ANNOTATION ----
# This section creates line segments between centroid pairs and annotates them
# with statistical significance based on pairwise PERMANOVA results

## Seasonality

pairwise_df_season <- results_seasonality$pairwise %>%
  mutate(pair_id = paste0(pmin(Group1, Group2), "_", pmax(Group1, Group2)))

centroid_segments_mode_season <- centroids_mode %>%
  filter(var_type == 'Seasonality' ) %>%
  dplyr::select(endem_type, mode_temp, mode_precip) %>%
  cross_join(., ., suffix = c("1", "2")) %>%
  filter(endem_type1 != endem_type2) %>%
  mutate(pair_id = paste(pmin(endem_type1, endem_type2), 
                         pmax(endem_type1, endem_type2), sep = "_")) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  rename(Group1 = endem_type1, Group2 = endem_type2,
         x1 = mode_temp1, y1 = mode_precip1,
         x2 = mode_temp2, y2 = mode_precip2) %>%
  left_join(dplyr::select(pairwise_df_season, pair_id, p_value), by = "pair_id") %>%
  mutate(significance = ifelse(!is.na(p_value) & p_value < 0.05, 
                               "significant", "non-significant"))%>%
  mutate(var_type = 'Seasonality')


## Annual

pairwise_df_annual <- results_annual$pairwise %>%
  mutate(pair_id = paste0(pmin(Group1, Group2), "_", pmax(Group1, Group2)))

centroid_segments_mode_annual <- centroids_mode %>%
  filter(var_type == 'Annual Conditions') %>%
  dplyr::select(endem_type, mode_temp, mode_precip) %>%
  cross_join(., ., suffix = c("1", "2")) %>%
  filter(endem_type1 != endem_type2) %>%
  mutate(pair_id = paste(pmin(endem_type1, endem_type2), 
                         pmax(endem_type1, endem_type2), sep = "_")) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  rename(Group1 = endem_type1, Group2 = endem_type2,
         x1 = mode_temp1, y1 = mode_precip1,
         x2 = mode_temp2, y2 = mode_precip2) %>%
  left_join(dplyr::select(pairwise_df_annual, pair_id, p_value), by = "pair_id") %>%
  mutate(significance = ifelse(!is.na(p_value) & p_value < 0.05, 
                               "significant", "non-significant")) %>%
  mutate(var_type = 'Annual Conditions')




# COLOR PALETTE

canape_colors  <- c('neo' = "#DC085F",
                 'paleo' = "#359EDB",
                 'not significant' = "grey60", 
                 'mixed'= "#D89EFA",
                 'super' = "#6C02A6")



# FINAL PLOT with names and axis identified for individual plots

label_type <- c(
  "Annual" = "Annual Conditions",
  "Seasonality" = "Seasonality"
)

label_time <- c(
  "present" = "Present",
  "lh" = "Late Holocene\n(4.2-0.3 ka)",
  "lgm" = "Last Glacial\nMaximum\n(ca. 21 ka)",
  "lig" = "Last Interglacial\n(ca. 130 ka)",
  "pliocene" = "Pliocene\n(ca. 3.3 Ma)") 



df_plots <- bind_rows(df_all) %>% 
  mutate(endem_type = factor(endem_type, levels = c('neo','paleo',
                                        'not significant',
                                        'mixed', 'super')))


plot_annual<- df_plots %>%
  filter(Precipitation<6000) %>%
  mutate(time = factor(time, 
                       levels = c('present', 'lh', 'lgm', 'lig', 'pliocene'))) %>%
  filter(var_type == 'Annual Conditions') %>%
  ggplot(aes(x = Temperature, y = Precipitation, 
             fill = endem_type, col = endem_type)) +
  geom_point(alpha = 0.5, size = 0.2) +
  stat_density_2d(geom = "polygon", alpha = 0.25, 
                  contour = TRUE,lwd = 0.1,
                  aes(alpha = after_stat(level),
                  )) +
  scale_color_manual(values = canape_colors) +
  scale_fill_manual(values = canape_colors) +
  theme_light()  +
  facet_grid(var_type~  time,   
             labeller = labeller(type = label_type, time = label_time, 7)) +
  
  # Replace X with open circle (shape 21) for KDE centroids
  geom_point(data = filter(centroids_mode,var_type == 'Annual Conditions'),
             aes(x = mode_temp, y = mode_precip,),
             color = "black", shape = 21, size = 2, stroke = 1.5, alpha = 0.75) +
  
  #coord_cartesian(expand = F) +
  # theme(legend.position = 'none') +
  labs(color = "CANAPE group",
       fill = "CANAPE group",
       x = 'Mean Annual Temperature (BIO 1)',
       y = 'Annual Precipitation (BIO 12)') +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 4) ) +
  coord_cartesian(expand = F)+
  theme(strip.text = element_text(size = 10))



plot_season<- df_plots %>%
  mutate(time = factor(time, 
                       levels = c('present', 'lh', 'lgm', 'lig', 'pliocene'))) %>%
  filter(var_type == 'Seasonality') %>%
  ggplot(aes(x = Temperature, y = Precipitation, 
             fill = endem_type, col = endem_type)) +
  geom_point(alpha = 0.5, size = 0.2) +
  stat_density_2d(geom = "polygon", alpha = 0.25, 
                  contour = TRUE,lwd = 0.1,
                  aes(alpha = after_stat(level),
                  )) +
  scale_color_manual(values = canape_colors) +
  scale_fill_manual(values = canape_colors) +
  theme_light()  +
  facet_grid(var_type~  time,    
             labeller = labeller(type = label_type, time = label_time,7))+
  
  # Replace X with open circle (shape 21) for KDE centroids
  geom_point(data = filter(centroids_mode,var_type == 'Seasonality'),
             aes(x = mode_temp, y = mode_precip,),
             color = "black", shape = 21, size = 2, stroke = 1.5, alpha = 0.75) +
  
  #coord_cartesian(expand = F) +
  # theme(legend.position = 'none') +
  labs(color = "CANAPE group",
       fill = "CANAPE group",
       x = 'Temperature Seasonality (BIO 4)',
       y = 'Precipitation Seasonality (BIO 15)') +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 4) ) +
  coord_cartesian(expand = F) +
  theme(strip.text = element_text(size = 10))


((plot_annual / plot_season) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') + 
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))  %>%
  ggsave(filename = 'figure_all_wide.tiff', plot = ., device = 'tiff',
         dpi = 300, width = 168,height = 100,units = 'mm', scale = 2, bg = 'white')
