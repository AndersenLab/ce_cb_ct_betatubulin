# Materials Key
# This script maps script variables to figure/table IDs
# and creates the output directories and file names with the correct syntax
#
# Naming conventions:
#   Main figures:
#     - Folder: Figure_X (e.g., Figure_1, Figure_5)
#     - Filename: FIGX.{tiff,png,jpg} (e.g., FIG1.tiff, FIG5.png)
#   Supplementary figures:
#     - Folder: S#_FIG (e.g., S1_FIG, S19_FIG)
#     - Filename: S#_FIG.{tiff,png,jpg} (e.g., S1_FIG.tiff, S19_FIG.png)
#   Panel figures:
#     - Main: Folder Figure_X/Figure_Xa, Filename FIGXa.{tiff,png,jpg}
#     - Supp: Folder S#_FIG/S#_FIGa, Filename S#_FIGa.{tiff,png,jpg}

library(here)


# Main figures list with variable name mappings
# Format: variable_name = figure_number (integer)
# For standalone figures (not panels)
figures <- list(
  # Add main figures here as needed
)

# Tables list with variable name mappings
tables <- list()

# Supplementary figures list with variable name mappings
# Format: variable_name = supplementary_figure_number (integer)
# For standalone figures (not panels)
sup_figures <- list(
  # S1_FIG - ben-1 expression x ABZ response (from beta_tub_expression_plots.R)
  main_figure = 1,

  # S2_FIG - C. briggsae HTA strains tree (from all_sp_tubulin_trees.R)
  cb_hta_tree_final = 2,

  # S3_FIG - C. tropicalis HTA strains tree (from all_sp_tubulin_trees.R)
  ct_hta_tree_final = 3,

  # S19_FIG - other beta tubulin expression plots (from beta_tub_expression_plots.R)
  final_plot_patchwork = 19,

  # S20_FIG - tubulin alignments MSA (from align_tub.R)
  combined_plot = 20,

  # S21_FIG - C. elegans missense BLOSUM/Grantham (from pull_bg.R)
  ce_plot = 21,

  # S26_FIG - combined missense panels (from combine_missense_panels.R)
  missense_combined = 26
)

# Panel figures - figures composed of multiple panels (map + trees)
# Panel naming convention:
#   - Panel A: map
#   - Panel B: C. elegans (ce) tree
#   - Panel C: C. briggsae (cb) tree
#   - Panel D: C. tropicalis (ct) tree
# Format: list(figure_num, panel_letter, is_supplementary)

figure_panels <- list(

  # Figure_5 - ben-1 trees and map (main figure)
  ce_ben1_tree_plot = list(figure_num = 5, panel = "b", is_supplementary = FALSE),
  cb_ben1_tree_plot = list(figure_num = 5, panel = "c", is_supplementary = FALSE),
  ct_ben1_tree_plot = list(figure_num = 5, panel = "d", is_supplementary = FALSE),

  # S22_FIG - tbb-1 trees and map
  ce_tbb1_tree_plot = list(figure_num = 22, panel = "b", is_supplementary = TRUE),
  cb_tbb1_tree_plot = list(figure_num = 22, panel = "c", is_supplementary = TRUE),
  ct_tbb1_tree_plot = list(figure_num = 22, panel = "d", is_supplementary = TRUE),

  # S23_FIG - tbb-2 trees and map
  ce_tbb2_tree_plot = list(figure_num = 23, panel = "b", is_supplementary = TRUE),
  cb_tbb2_tree_plot = list(figure_num = 23, panel = "c", is_supplementary = TRUE),
  ct_tbb2_tree_plot = list(figure_num = 23, panel = "d", is_supplementary = TRUE),

  # S24_FIG - mec-7 trees and map
  ce_mec7_tree_plot = list(figure_num = 24, panel = "b", is_supplementary = TRUE),
  cb_mec7_tree_plot = list(figure_num = 24, panel = "c", is_supplementary = TRUE),
  ct_mec7_tree_plot = list(figure_num = 24, panel = "d", is_supplementary = TRUE),

  # S25_FIG - tbb-4 trees and map
  ce_tbb4_tree_plot = list(figure_num = 25, panel = "b", is_supplementary = TRUE),
  cb_tbb4_tree_plot = list(figure_num = 25, panel = "c", is_supplementary = TRUE),
  ct_tbb4_tree_plot = list(figure_num = 25, panel = "d", is_supplementary = TRUE)
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

# Generate file names for main figures
# Folder: Figure_X, Filename: FIGX.{ext}
generate_main_figure_fns <- function(figures) {
  file_names <- list()
  for (name in names(figures)) {
    fig_num <- figures[[name]]
    folder_name <- paste0("Figure_", fig_num)
    file_base <- paste0("FIG", fig_num)
    folder_path <- here("figures", folder_name)
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    file_names[[name]] <- list(
      png = here("figures", folder_name, paste0(file_base, ".png")),
      jpg = here("figures", folder_name, paste0(file_base, ".jpg")),
      tiff = here("figures", folder_name, paste0(file_base, ".tiff"))
    )
  }
  return(file_names)
}

# Generate file names for supplementary figures
# Folder: S#_FIG, Filename: S#_FIG.{ext}
generate_sup_figure_fns <- function(sup_figures) {
  file_names <- list()
  for (name in names(sup_figures)) {
    fig_num <- sup_figures[[name]]
    folder_name <- paste0("S", fig_num, "_FIG")
    file_base <- folder_name
    folder_path <- here("figures", folder_name)
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    file_names[[name]] <- list(
      png = here("figures", folder_name, paste0(file_base, ".png")),
      jpg = here("figures", folder_name, paste0(file_base, ".jpg")),
      tiff = here("figures", folder_name, paste0(file_base, ".tiff"))
    )
  }
  return(file_names)
}

figure_fns <- generate_main_figure_fns(figures)
sup_figure_fns <- generate_sup_figure_fns(sup_figures)

# Generate file names for panel figures
# Main panels: Folder Figure_X/Figure_Xa, Filename FIGXa.{ext}
# Supp panels: Folder S#_FIG/S#_FIGa, Filename S#_FIGa.{ext}
# If panel is NULL, saves directly to figure directory (standalone)
generate_panel_fns <- function(panels) {
  file_names <- list()
  for (name in names(panels)) {
    fig_num <- panels[[name]]$figure_num
    panel_letter <- panels[[name]]$panel
    is_sup <- panels[[name]]$is_supplementary

    if (is_sup) {
      # Supplementary figure
      parent_folder <- paste0("S", fig_num, "_FIG")
      if (is.null(panel_letter)) {
        # Standalone (no panel)
        folder_path <- here("figures", parent_folder)
        file_base <- parent_folder
      } else {
        # Panel figure
        panel_folder <- paste0(parent_folder, panel_letter)
        folder_path <- here("figures", parent_folder, panel_folder)
        file_base <- panel_folder
      }
    } else {
      # Main figure
      parent_folder <- paste0("Figure_", fig_num)
      if (is.null(panel_letter)) {
        # Standalone (no panel)
        folder_path <- here("figures", parent_folder)
        file_base <- paste0("FIG", fig_num)
      } else {
        # Panel figure
        panel_folder <- paste0("Figure_", fig_num, panel_letter)
        folder_path <- here("figures", parent_folder, panel_folder)
        file_base <- paste0("FIG", fig_num, panel_letter)
      }
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
