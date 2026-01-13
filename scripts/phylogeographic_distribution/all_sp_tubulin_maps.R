# Code to generate main figure for manuscript
# Figure N
# Panel A: World map with strains from all 3 species with high-impact variants
# Panel B: C. elegans tree with strains colored by high-impact variant
# Panel C: C. briggsae tree with strains colored by high-impact variant
# Panel D: C. tropicalis tree with strains colored by high-impact variant

#### Load packages ####
source("scripts/phylogeographic_distribution/plot_maps.R")

# Load standardized output functions and figure mappings
source("bin/outs.R")
source("bin/materials_key.R")

#### Inputs ####
isotype_folder_id <- "20250128"

#### Outputs ####
# Output file paths are now managed by materials_key.R via panel_fns

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

# mec-7
ce_mec7_var_df <- filter_high_impact_variants(ce_df, "mec-7_clean_call", "C. elegans")
cb_mec7_var_df <- filter_high_impact_variants(cb_df, "mec-7_clean_call", "C. briggsae")
ct_mec7_var_df <- filter_high_impact_variants(ct_df, "mec-7_clean_call", "C. tropicalis")

# tbb-4
ce_tbb4_var_df <- filter_high_impact_variants(ce_df, "tbb-4_clean_call", "C. elegans")
cb_tbb4_var_df <- filter_high_impact_variants(cb_df, "tbb-4_clean_call", "C. briggsae")
ct_tbb4_var_df <- filter_high_impact_variants(ct_df, "tbb-4_clean_call", "C. tropicalis")

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

mec7_var_df <- combine_filtered_data_frames(
  ce_mec7_var_df,
  cb_mec7_var_df,
  ct_mec7_var_df
)

tbb4_var_df <- combine_filtered_data_frames(
  ce_tbb4_var_df,
  cb_tbb4_var_df,
  ct_tbb4_var_df
)

#### Create sf obj. from combined data ####

## ben-1
ben1_var_sf <- convert_to_sf(ben1_var_df)

## tbb-1
tbb1_var_sf <- convert_to_sf(tbb1_var_df)

## tbb-2
tbb2_var_sf <- convert_to_sf(tbb2_var_df)

## mec-7
mec7_var_sf <- convert_to_sf(mec7_var_df)

## tbb-4
tbb4_var_sf <- convert_to_sf(tbb4_var_df)

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

mec7_map <- plot_high_impact_variants_map(
  world = mec7_var_sf$world,
  all_var_sf = mec7_var_sf$all_var_sf
)

tbb4_map <- plot_high_impact_variants_map(
  world = tbb4_var_sf$world,
  all_var_sf = tbb4_var_sf$all_var_sf
)

#### Save plots ####

# ben-1 map (Figure 5, panel a)
save_plot(
  tplot = ben1_map,
  fn_list = panel_fns$ben1_map,
  w_in = 7.5,
  h_in = 4
)

# tbb-1 map (S26_FIG, panel a)
save_plot(
  tplot = tbb1_map,
  fn_list = panel_fns$tbb1_map,
  w_in = 7.5,
  h_in = 4
)

# tbb-2 map (S27_FIG, panel a)
save_plot(
  tplot = tbb2_map,
  fn_list = panel_fns$tbb2_map,
  w_in = 7.5,
  h_in = 4
)

# mec-7 map (S28_FIG, panel a)
save_plot(
  tplot = mec7_map,
  fn_list = panel_fns$mec7_map,
  w_in = 7.5,
  h_in = 4
)

# tbb-4 map (S29_FIG, panel a)
save_plot(
  tplot = tbb4_map,
  fn_list = panel_fns$tbb4_map,
  w_in = 7.5,
  h_in = 4
)
