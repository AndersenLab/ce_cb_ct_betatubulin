# Script where function to calculate resistance threshold for strains in an 
# assay is defined

# used in several scripts in the project
# - `normalize_phenos.R`
# - `beta_tub_expression_plots.R`
# - `hta_res_thresholds.R`


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

# #' Get Resistant Strains
# #'
# #' This function calculates the resistance threshold and identifies strains with resistance based on a specified phenotype column.
# #'
# #' @param exp_summary_df A data frame containing the experimental summary with columns `strain` and a dynamic phenotype column.
# #' @param pheno_col A string specifying the name of the phenotype column to be used for resistance calculation.
# #' @param threshold_per A numeric value representing the threshold percentage to classify strains as resistant. Default is 0.50.
# #'
# #' @return A list containing:
# #' \item{min}{The minimum value of the specified phenotype column.}
# #' \item{max}{The maximum value of the specified phenotype column.}
# #' \item{res_strains}{A vector of strains identified as resistant based on the threshold.}
# #' \item{threshold}{The calculated threshold value used for resistance classification.}
# #'
# # Set a resistance threshold using the normalized phenotypes
# # Function to calculate the resistance threshold and identify strains with resistance
# # function expects a data frame with columns `strain` and a dynamic `pheno_col`
# get_res_strains <- function(exp_summary_df, pheno_col, threshold_per = 0.50) {
#   # exp_summary_df: a data frame with the columns `strain` and a dynamic `pheno_col`
#   max <- exp_summary_df %>%
#     dplyr::pull(!!sym(pheno_col)) %>%
#     max()

#   min <- exp_summary_df %>%
#     dplyr::pull(!!sym(pheno_col)) %>%
#     min()

#   threshold <- abs(max - min) * threshold_per # .50 threshold only returns 7 strains

#   res_strains <- exp_summary_df %>%
#     dplyr::filter(!!sym(pheno_col) > threshold) %>%
#     dplyr::pull(strain)

#   out <- list(
#     min = min,
#     max = max,
#     res_strains = res_strains,
#     threshold = threshold
#   )

#   return(out)
# }

