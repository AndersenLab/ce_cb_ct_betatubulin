library(dplyr)
library(data.table)

set.seed(123)

#### Functions ####

# Function also used in `normalize_phenotypes.R`
# to flag strains as resistant or sensitive
#' Get Resistant Strains
#'
#' This function calculates the resistance threshold and identifies strains with resistance based on a specified phenotype column.
#'
#' @param exp_summary_df A data frame containing the experimental summary with columns `strain` and a dynamic phenotype column.
#' @param pheno_col A string specifying the name of the phenotype column to be used for resistance calculation.
#' @param threshold_per A numeric value representing the threshold percentage to classify strains as resistant. Default is 0.50.
#' @param ref_strain Reference strain to calculate the resistance threshold.
#'
#' @return A list containing:
#' \item{min}{The minimum value of the specified phenotype column.}
#' \item{max}{The maximum value of the specified phenotype column.}
#' \item{res_strains}{A vector of strains identified as resistant based on the threshold.}
#' \item{threshold}{The calculated threshold value used for resistance classification.}
#'
# Set a resistance threshold using the normalized phenotypes
# Function to calculate the resistance threshold and identify strains with resistance
# function expects a data frame with columns `strain` and a dynamic `pheno_col`
get_res_strains_ref <- function(exp_summary_df, pheno_col, ref_strain, threshold_per = 0.50) {
  # exp_summary_df: a data frame with the columns `strain` and a dynamic `pheno_col`
  # max <- exp_summary_df %>%
  #   dplyr::pull(!!sym(pheno_col)) %>%
  #   max()

  min <- exp_summary_df %>%
    dplyr::filter(strain == ref_strain) %>%
    dplyr::pull(!!sym(pheno_col)) %>%
    min()

  res_thres <- (min * -threshold_per) + min  

  res_strains <- exp_summary_df %>%
    dplyr::filter(!!sym(pheno_col) > res_thres) %>%
    dplyr::pull(strain)

  out <- list(
    min = min,
    max = max,
    res_strains = res_strains,
    threshold = res_thres
  )

  return(out)
}

#### Output ####

# output directory named with the summary
data_out_dir <- glue::glue(
  "data/proc/hta_res_threshold"
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
    "data/raw/species_hta/elegans_df_ABZ_summary_mean_median_wormlength.csv"

cb_hta_file <-
    "data/raw/species_hta/briggsae_df_ABZ_summary_mean_median_wormlength.csv"

ct_hta_file <-
    "data/raw/species_hta/tropicalis_df_ABZ_summary_mean_median_wormlength.csv"

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

