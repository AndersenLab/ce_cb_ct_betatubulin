# Materials Key
# This script maps script variables to figure/table IDs
# and creates the output directories and file names with the correct syntax

library(here)

# Tables list with variable name mappings
tables <- list(
  # Add table mappings here as needed
  # Example: blossum_grantham_table = "table_S13"
)

# Main figures list with variable name mappings
figures <- list(
  # figure_9 panels (ben-1 trees from all_sp_tubulin_trees.R)
  ce_ben1_tree_plot = "figure_5",
  cb_ben1_tree_plot = "figure_5",
  ct_ben1_tree_plot = "figure_5"
)

# Supplementary figures list with variable name mappings
sup_figures <- list(
  # S1_FIG - ben-1 expression x ABZ response (from beta_tub_expression_plots.R)
  main_figure = "S1_FIG",

  # figure_S4 - C. briggsae HTA strains tree (from all_sp_tubulin_trees.R)
  cb_hta_tree_plot = "figure_S4",

  # figure_S5 - C. tropicalis HTA strains tree (from all_sp_tubulin_trees.R)
  ct_hta_tree_plot = "figure_S5",

  # S19_FIG - other beta tubulin expression plots (from beta_tub_expression_plots.R)
  final_plot_patchwork = "S19_FIG",

  # figure_S20 - tubulin alignments MSA (from align_tub.R)
  combined_plot = "figure_S20",

  # figure_S21 - C. elegans missense BLOSUM/Grantham (from pull_bg.R)
  ce_plot = "figure_S21",

  # figure_S22 - tbb-1 trees (from all_sp_tubulin_trees.R)
  ce_tbb1_tree_plot = "figure_S22",
  cb_tbb1_tree_plot = "figure_S22",
  ct_tbb1_tree_plot = "figure_S22",

  # figure_S23 - tbb-2 trees (from all_sp_tubulin_trees.R)
  ce_tbb2_tree_plot = "figure_S23",
  cb_tbb2_tree_plot = "figure_S23",
  ct_tbb2_tree_plot = "figure_S23",

  # figure_S24 - mec-7 trees (from all_sp_tubulin_trees.R)
  ce_mec7_tree_plot = "figure_S24",
  cb_mec7_tree_plot = "figure_S24",
  ct_mec7_tree_plot = "figure_S24",

  # figure_S25 - tbb-4 trees (from all_sp_tubulin_trees.R)
  ce_tbb4_tree_plot = "figure_S25",
  cb_tbb4_tree_plot = "figure_S25",
  ct_tbb4_tree_plot = "figure_S25",

  # figure_S26 - combined missense panels (from combine_missense_panels.R)
  missense_combined = "figure_S26"
)

# Generate file names for the tables
generate_table_fns <- function(tables) {
  file_names <- list()
  for (name in names(tables)) {
    folder_path <- here("tables", tables[[name]])
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    file_names[[name]] <- list(
      html = here("tables", tables[[name]], paste0(tables[[name]], ".html")),
      tsv = here("tables", tables[[name]], paste0(tables[[name]], ".tsv.zip")),
      docx = here("tables", tables[[name]], paste0(tables[[name]], ".docx")),
      csv = here("tables", tables[[name]], paste0(tables[[name]], ".csv"))
    )
  }
  return(file_names)
}

table_fns <- generate_table_fns(tables)

# Generate file names for the figures
generate_figure_fns <- function(figures) {
  file_names <- list()
  for (name in names(figures)) {
    folder_path <- here("figures", figures[[name]])
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    file_names[[name]] <- list(
      png = here("figures", figures[[name]], paste0(figures[[name]], ".png")),
      jpg = here("figures", figures[[name]], paste0(figures[[name]], ".jpg")),
      tiff = here("figures", figures[[name]], paste0(figures[[name]], ".tiff"))
    )
  }
  return(file_names)
}

figure_fns <- generate_figure_fns(figures)

sup_figure_fns <- generate_figure_fns(sup_figures)
