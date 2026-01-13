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
# For standalone figures (not panels)
figures <- list(
  # Add main figures here as needed
)

# Supplementary figures list with variable name mappings
# For standalone figures (not panels)
sup_figures <- list(
  # S1_FIG - ben-1 expression x ABZ response (from beta_tub_expression_plots.R)
  main_figure = "S1_FIG",

  # S19_FIG - other beta tubulin expression plots (from beta_tub_expression_plots.R)
  final_plot_patchwork = "S19_FIG",

  # figure_S20 - tubulin alignments MSA (from align_tub.R)
  combined_plot = "figure_S20",

  # figure_S21 - C. elegans missense BLOSUM/Grantham (from pull_bg.R)
  ce_plot = "figure_S21",

  # figure_S26 - combined missense panels (from combine_missense_panels.R)
  missense_combined = "figure_S26"
)

# Panel figures - figures composed of multiple panels (map + trees)
# Panel naming convention:
#   - Panel A: map
#   - Panel B: C. elegans (ce) tree
#   - Panel C: C. briggsae (cb) tree
#   - Panel D: C. tropicalis (ct) tree
# Format: list(figure_id, panel_letter)

figure_panels <- list(

  # figure_5 - ben-1 trees and map (main figure)
  ce_ben1_tree_plot = list(figure_id = "figure_5", panel = "b"),
  cb_ben1_tree_plot = list(figure_id = "figure_5", panel = "c"),
  ct_ben1_tree_plot = list(figure_id = "figure_5", panel = "d"),

  # figure_S4 - C. briggsae HTA strains tree (standalone, no panels)
  cb_hta_tree_plot = list(figure_id = "figure_S4", panel = NULL),

  # figure_S5 - C. tropicalis HTA strains tree (standalone, no panels)
  ct_hta_tree_plot = list(figure_id = "figure_S5", panel = NULL),

  # figure_S22 - tbb-1 trees and map
  ce_tbb1_tree_plot = list(figure_id = "figure_S22", panel = "b"),
  cb_tbb1_tree_plot = list(figure_id = "figure_S22", panel = "c"),
  ct_tbb1_tree_plot = list(figure_id = "figure_S22", panel = "d"),

  # figure_S23 - tbb-2 trees and map
  ce_tbb2_tree_plot = list(figure_id = "figure_S23", panel = "b"),
  cb_tbb2_tree_plot = list(figure_id = "figure_S23", panel = "c"),
  ct_tbb2_tree_plot = list(figure_id = "figure_S23", panel = "d"),

  # figure_S24 - mec-7 trees and map
  ce_mec7_tree_plot = list(figure_id = "figure_S24", panel = "b"),
  cb_mec7_tree_plot = list(figure_id = "figure_S24", panel = "c"),
  ct_mec7_tree_plot = list(figure_id = "figure_S24", panel = "d"),

  # figure_S25 - tbb-4 trees and map
  ce_tbb4_tree_plot = list(figure_id = "figure_S25", panel = "b"),
  cb_tbb4_tree_plot = list(figure_id = "figure_S25", panel = "c"),
  ct_tbb4_tree_plot = list(figure_id = "figure_S25", panel = "d")
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

# Generate file names for standalone figures (not panels)
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

# Generate file names for panel figures
# Panels are saved in subdirectories: figures/figure_X/figure_Xa/figure_Xa.{png,jpg,tiff}
# If panel is NULL, saves directly to figure directory (standalone)
generate_panel_fns <- function(panels) {
  file_names <- list()
  for (name in names(panels)) {
    figure_id <- panels[[name]]$figure_id
    panel_letter <- panels[[name]]$panel

    if (is.null(panel_letter)) {
      # Standalone figure (no panel subdirectory)
      folder_path <- here("figures", figure_id)
      file_base <- figure_id
    } else {
      # Panel figure - create subdirectory
      panel_name <- paste0(figure_id, panel_letter)
      folder_path <- here("figures", figure_id, panel_name)
      file_base <- panel_name
    }

    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }

    file_names[[name]] <- list(
      png = here(folder_path, paste0(file_base, ".png")),
      jpg = here(folder_path, paste0(file_base, ".jpg")),
      tiff = here(folder_path, paste0(file_base, ".tiff"))
    )
  }
  return(file_names)
}

panel_fns <- generate_panel_fns(figure_panels)
