# Code to generate main figure for manuscript
# Figure N
# Panel A: World map with strains from all 3 species with high-impact variants
# Panel B: C. elegans tree with strains colored by high-impact variant
# Panel C: C. briggsae tree with strains colored by high-impact variant
# Panel D: C. tropicalis tree with strains colored by high-impact variant

#### Load packages ####
source("scripts/phylogeographic_distribution/plot_maps.R")

#### Inputs ####
isotype_folder_id <- "20250128"

#### Outputs ####


ben1_maps <- list(
  png = "figures/figure_6/figure_6a.png",
  eps = "figures/figure_6/figure_6a.eps"
)

tbb2_maps <- list(
  png = "figures/figure_S10/figure_S10a.png",
  eps = "figures/figure_S10/figure_S10a.eps"
)

tbb1_maps <- list(
  png = "figures/figure_S11/figure_S11a.png",
  eps = "figures/figure_S11/figure_S11a.eps"
)

#### Load data ####

# Load isotype variant summary data
ce_var_file <-
  glue::glue(
    "data/isotype_variant_table/c_elegans/{isotype_folder_id}/isotype_variant_summary.tsv"
  )
cb_var_file <-
  glue::glue(
    "data/isotype_variant_table/c_briggsae/{isotype_folder_id}/isotype_variant_summary.tsv"
  )

ct_var_file <-
  glue::glue(
    "data/isotype_variant_table/c_tropicalis/{isotype_folder_id}/isotype_variant_summary.tsv"
  )

ce_df <- data.table::fread(ce_var_file) %>%
  # clean ben-1 calls added from concordance check
  dplyr::mutate(
    `ben-1_clean_call` = case_when(
      strain == "JU3125" ~ "Transposon insertion",
      strain == "ED3011" ~ "Splice Donor",
      TRUE ~ `ben-1_clean_call`
    )
  )

cb_df <- data.table::fread(cb_var_file)
ct_df <- data.table::fread(ct_var_file)

#### Filter dataframes for each tubulin gene ####

# ben-1
ce_ben1_var_df <- filter_high_impact_variants(ce_df, "ben-1_clean_call", "C. elegans")
cb_ben1_var_df <- filter_high_impact_variants(cb_df, "ben-1_clean_call", "C. briggsae")
ct_ben1_var_df <- filter_high_impact_variants(ct_df, "ben-1_clean_call", "C. tropicalis")

# tbb-1
ce_tbb1_var_df <- filter_high_impact_variants(ce_df, "tbb-1_clean_call", "C. elegans")
cb_tbb1_var_df <- filter_high_impact_variants(cb_df, "tbb-1_clean_call", "C. briggsae")
ct_tbb1_var_df <- filter_high_impact_variants(ct_df, "tbb-1_clean_call", "C. tropicalis")

# tbb-2
ce_tbb2_var_df <- filter_high_impact_variants(ce_df, "tbb-2_clean_call", "C. elegans")
cb_tbb2_var_df <- filter_high_impact_variants(cb_df, "tbb-2_clean_call", "C. briggsae")
ct_tbb2_var_df <- filter_high_impact_variants(ct_df, "tbb-2_clean_call", "C. tropicalis")

#### Combine species data ####

## ben-1 ##
ben1_var_df <- combine_filtered_data_frames(
  ce_ben1_var_df,
  cb_ben1_var_df,
  ct_ben1_var_df
)

tbb1_var_df <- combine_filtered_data_frames(
  ce_tbb1_var_df,
  cb_tbb1_var_df,
  ct_tbb1_var_df
)

tbb2_var_df <- combine_filtered_data_frames(
  ce_tbb2_var_df,
  cb_tbb2_var_df,
  ct_tbb2_var_df
)

#### Create sf obj. from combined data ####

## ben-1
ben1_var_sf <- convert_to_sf(ben1_var_df)

## tbb-1
tbb1_var_sf <- convert_to_sf(tbb1_var_df)

## tbb-2
tbb2_var_sf <- convert_to_sf(tbb2_var_df)




#### Generate map plots ####
ben1_map <- plot_high_impact_variants_map(
  world = ben1_var_sf$world,
  all_var_sf = ben1_var_sf$all_var_sf
)

tbb1_map <- plot_high_impact_variants_map(
  world = tbb1_var_sf$world,
  all_var_sf = tbb1_var_sf$all_var_sf
)

tbb2_map <- plot_high_impact_variants_map(
  world = tbb2_var_sf$world,
  all_var_sf = tbb2_var_sf$all_var_sf
)

#### Save plots ####

print(
  'saving ben-1 map'
)

# Ben-1
save_plot(
  fn_list = ben1_maps,
  plot = ben1_map
)

print(
  'saving tbb-2 map'
)

# tbb-2
save_plot(
  fn_list = tbb2_maps,
  plot = tbb2_map
)

print(
  'saving tbb-1 map'
)

# tbb-1
save_plot(
  fn_list = tbb1_maps,
  plot = tbb1_map
)

