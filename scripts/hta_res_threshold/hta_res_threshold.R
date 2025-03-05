library(dplyr)
library(data.table)

set.seed(123)

source("bin/get_res_threshold.R")

#### Output ####

# output directory named with the summary
data_out_dir <- glue::glue(
  "data/hta_res_threshold"
)

# if it does exist, delete the contents
if (dir.exists(data_out_dir)) {
  unlink(data_out_dir, recursive = TRUE)
}
# if it doesn't exist, create it
dir.create(data_out_dir, recursive = TRUE)

#### Input ####

ref_threshold <- 0.75

ce_hta_file <-
    "data/hta_summaries/elegans_df_ABZ_summary_mean_median_wormlength.csv"

cb_hta_file <-
    "data/hta_summaries/briggsae_df_ABZ_summary_mean_median_wormlength.csv"

ct_hta_file <-
    "data/hta_summaries//tropicalis_df_ABZ_summary_mean_median_wormlength.csv"

ce_hta <- fread(ce_hta_file)
cb_hta <- fread(cb_hta_file)
ct_hta <- fread(ct_hta_file)

#### Calculate Resistance Thresholds ####

ce_res <- get_res_strains_ref(
    exp_summary_df = ce_hta,
    pheno_col = "mean_median_wormlength_um_delta_reg",
    ref_strain = "N2",
    threshold_per = ref_threshold
)

# C. briggsae
cb_res <- get_res_strains_ref(
    exp_summary_df = cb_hta,
    pheno_col = "mean_median_wormlength_um_delta_reg",
    ref_strain = "AF16",
    threshold_per = ref_threshold
)

# C. tropicalis
ct_res <- get_res_strains_ref(
    exp_summary_df = ct_hta,
    pheno_col = "mean_median_wormlength_um_reg_delta",
    ref_strain = "NIC58",
    threshold_per = ref_threshold
)

#### Save the data ####
write_res_strains_csv <- function(ce_output, cb_output, ct_output, file_path) {
  # Create data frames for each species
  ce_df <- data.frame(
    species = "c_elegans",
    hta_threshold = ce_output$threshold
  )
  
  cb_df <- data.frame(
    species = "c_briggsae",
    hta_threshold = cb_output$threshold
  )
  
  ct_df <- data.frame(
    species = "c_tropicalis",
    hta_threshold = ct_output$threshold
  )
  
  # Combine the data frames
  combined_df <- rbind(ce_df, cb_df, ct_df)
  
  # Write the combined data frame to a .csv file
  write.csv(combined_df, file = file_path, row.names = FALSE)
}

# Assuming `ce_res`, `cb_res`, and `ct_res` are the outputs from `get_res_strains_ref` function
write_res_strains_csv(
  ce_res, cb_res, ct_res,
  file.path(data_out_dir, "res_thresholds.csv"))

