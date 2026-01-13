# Code to generate trees for all species denoting high-impact beta-tubulin variants
set.seed(123)

#### Load required libraries ####
library(ggtree)

# Load standardized output functions and figure mappings
source("bin/outs.R")
source("bin/materials_key.R")


#### Define parameters ####
isotype_folder_id <- "20250128"

ce_tree_file <-
  "data/caendr_species_trees/c_elegans/20231213/WI.20231213.hard-filter.isotype.min4.tree"

cb_tree_file <-
  "data/caendr_species_trees/c_briggsae/20240129/WI.20240129.hard-filter.isotype.min4.tree"

ct_tree_file <-
  "data/caendr_species_trees/c_tropicalis/20231201/WI.20231201.hard-filter.isotype.min4.tree"

#### Define outputs ####
data_out_dir <- glue::glue(
  "data/proc/phylogeographic_distribution/{isotype_folder_id}"
)

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

## Figure output paths are now managed by materials_key.R
# Tree plots use: figure_fns (for main figures) and sup_figure_fns (for supplementary figures)


#### Load the data ####

ce_iso_var_sum_df <- data.table::fread(ce_var_file) %>%
  dplyr::mutate(
    `ben-1_clean_call` = dplyr::case_when(
      strain == "JU3125" ~ "Transposon insertion",
      strain == "ED3011" ~ "Splice Donor",
      TRUE ~ `ben-1_clean_call`
    )
  )


cb_iso_var_sum_df <- data.table::fread(cb_var_file)

ct_iso_var_sum_df <- data.table::fread(ct_var_file)


ce_tree <- read.tree(ce_tree_file)

cb_tree <- read.tree(cb_tree_file)

ct_tree <- read.tree(ct_tree_file)


create_tree_anno_df <- function(iso_sum_df, col_id) {
  fil_df <- iso_sum_df %>%
    dplyr::select(label = strain, !!rlang::sym(col_id)) %>%
    dplyr::rename(!!gsub("-", "", col_id) := !!rlang::sym(col_id))

  return(fil_df)
}

source("bin/var_color_scale.R")

add_var_to_tree <- function(tree, tree_anno_df, col_scale, input_col, xpos = NULL, ypos = NULL) {
  # load the tree
  base_tree <- ggtree(
    tree,
    layout = "equal_angle",
    linewidth = 0.15
  )

  # add data to base tree
  base_tree_anno <- base_tree %<+% tree_anno_df +
    ggtree::geom_point(
      data = ggtree::td_filter(isTip & !!rlang::sym(input_col) != "No variant"),
      ggplot2::aes_string(fill = input_col),
      shape = 21,
      size = 1.5
    ) +
    ggplot2::scale_fill_manual(values = col_scale) +
    ggtree::geom_treescale(
      offset = 0.005,
      x = xpos,
      y = ypos,
      linesize = 0.25,
      fontsize = 5,
      family = "Helvetica"
    ) # Adjust offset, x, and y as needed
  return(base_tree_anno)
}

#### Save the trees ####
# save_tree function is now loaded from bin/outs.R


#### ben-1 Trees ####

## Ce ben-1 high impact tree 5a ##

# create data frame for ben-1
ce_ben1_anno_df <- create_tree_anno_df(
  ce_iso_var_sum_df,
  "ben-1_clean_call"
)

# add the variant data to the tree
ce_ben1_tree_plot <- add_var_to_tree(
  ce_tree,
  ce_ben1_anno_df,
  strain_var_colors,
  "ben1_clean_call",
  xpos = -0.05,
  ypos = -0.09
)

# save the tree
save_tree(
  ce_ben1_tree_plot + theme(legend.position = "none"),
  figure_fns$ce_ben1_tree_plot,
  2.5,
  4
)

## Cb ben-1 high impact tree 5b ##

# create data frame for ben-1
cb_ben1_anno_df <- create_tree_anno_df(
  cb_iso_var_sum_df,
  "ben-1_clean_call"
)

# add the variant data to the tree
cb_ben1_tree_plot <- add_var_to_tree(
  cb_tree,
  cb_ben1_anno_df,
  strain_var_colors,
  "ben1_clean_call",
  xpos = -0.01,
  ypos = -0.135
)

# save the tree
save_tree(
  cb_ben1_tree_plot + theme(legend.position = "none"),
  figure_fns$cb_ben1_tree_plot,
  2.5,
  4
)

## Ct ben-1 high impact tree 5c ##
# create data frame for ben-1
ct_ben1_anno_df <- create_tree_anno_df(
  ct_iso_var_sum_df,
  "ben-1_clean_call"
)

# add the variant data to the tree
ct_ben1_tree_plot <- add_var_to_tree(
  ct_tree,
  ct_ben1_anno_df,
  strain_var_colors,
  "ben1_clean_call",
  xpos = -0.01,
  ypos = -0.125
)

# save the tree
save_tree(
  ct_ben1_tree_plot + theme(legend.position = "none"),
  figure_fns$ct_ben1_tree_plot,
  2.5,
  4
)

#### tbb-2 trees ####

## Ce tbb-2 high impact tree S10a ##
# create data frame for tbb-2
ce_tbb2_anno_df <- create_tree_anno_df(
  ce_iso_var_sum_df,
  "tbb-2_clean_call"
)

# add the variant data to the tree
ce_tbb2_tree_plot <- add_var_to_tree(
  ce_tree,
  ce_tbb2_anno_df,
  strain_var_colors,
  "tbb2_clean_call",
  xpos = -0.05,
  ypos = -0.09
)

# save the tree
save_tree(
  ce_tbb2_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ce_tbb2_tree_plot,
  2.5,
  4
)

## Cb tbb-2 high impact tree S10b ##
# create data frame for tbb-2
cb_tbb2_anno_df <- create_tree_anno_df(
  cb_iso_var_sum_df,
  "tbb-2_clean_call"
)

# add the variant data to the tree
cb_tbb2_tree_plot <- add_var_to_tree(
  cb_tree,
  cb_tbb2_anno_df,
  strain_var_colors,
  "tbb2_clean_call",
  xpos = -0.01,
  ypos = -0.135
)

# save the tree
save_tree(
  cb_tbb2_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$cb_tbb2_tree_plot,
  2.5,
  4
)

## Ct tbb-2 high impact tree S10c ##
# create data frame for tbb-2
ct_tbb2_anno_df <- create_tree_anno_df(
  ct_iso_var_sum_df,
  "tbb-2_clean_call"
)

# add the variant data to the tree
ct_tbb2_tree_plot <- add_var_to_tree(
  ct_tree,
  ct_tbb2_anno_df,
  strain_var_colors,
  "tbb2_clean_call",
  xpos = -0.01,
  ypos = -0.125
)

# save the tree
save_tree(
  ct_tbb2_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ct_tbb2_tree_plot,
  2.5,
  4
)

#### tbb-1 trees ####

## Ce tbb-1  ##

ce_tbb1_anno_df <- create_tree_anno_df(
  ce_iso_var_sum_df,
  "tbb-1_clean_call"
)

# add the variant data to the tree
ce_tbb1_tree_plot <- add_var_to_tree(
  ce_tree,
  ce_tbb1_anno_df,
  strain_var_colors,
  "tbb1_clean_call",
  xpos = -0.05,
  ypos = -0.09
)

# save the tree
save_tree(
  ce_tbb1_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ce_tbb1_tree_plot,
  2.5,
  4
)


## Cb tbb-1 tree ##

cb_tbb1_anno_df <- create_tree_anno_df(
  cb_iso_var_sum_df,
  "tbb-1_clean_call"
)


# add the variant data to the tree
cb_tbb1_tree_plot <- add_var_to_tree(
  cb_tree,
  cb_tbb1_anno_df,
  strain_var_colors,
  "tbb1_clean_call",
  xpos = -0.01,
  ypos = -0.135
)

# save the tree
save_tree(
  cb_tbb1_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$cb_tbb1_tree_plot,
  2.5,
  4
)

## Ct tbb-2 high impact tree S11c ##
# create data frame for tbb-2
ct_tbb1_anno_df <- create_tree_anno_df(
  ct_iso_var_sum_df,
  "tbb-1_clean_call"
)

# add the variant data to the tree
ct_tbb1_tree_plot <- add_var_to_tree(
  ct_tree,
  ct_tbb1_anno_df,
  strain_var_colors,
  "tbb1_clean_call",
  xpos = -0.01,
  ypos = -0.125
)

# save the tree
save_tree(
  ct_tbb1_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ct_tbb1_tree_plot,
  2.5,
  4
)

#### mec-7 trees ####

ce_mec7_anno_df <- create_tree_anno_df(
  ce_iso_var_sum_df,
  "mec-7_clean_call"
)

ce_mec7_tree_plot <- add_var_to_tree(
  ce_tree,
  ce_mec7_anno_df,
  strain_var_colors,
  "mec7_clean_call",
  xpos = -0.05,
  ypos = -0.09
)

# save the tree
save_tree(
  ce_mec7_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ce_mec7_tree_plot,
  2.5,
  4
)

cb_mec7_anno_df <- create_tree_anno_df(
  cb_iso_var_sum_df,
  "mec-7_clean_call"
)

cb_mec7_tree_plot <- add_var_to_tree(
  cb_tree,
  cb_mec7_anno_df,
  strain_var_colors,
  "mec7_clean_call",
  xpos = -0.01,
  ypos = -0.135
)

# save the tree
save_tree(
  cb_mec7_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$cb_mec7_tree_plot,
  2.5,
  4
)

ct_mec7_anno_df <- create_tree_anno_df(
  ct_iso_var_sum_df,
  "mec-7_clean_call"
)

ct_mec7_tree_plot <- add_var_to_tree(
  ct_tree,
  ct_mec7_anno_df,
  strain_var_colors,
  "mec7_clean_call",
  xpos = -0.01,
  ypos = -0.125
)

# save the tree
save_tree(
  ct_mec7_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ct_mec7_tree_plot,
  2.5,
  4
)

### tbb-4 trees ####

ce_tbb4_anno_df <- create_tree_anno_df(
  ce_iso_var_sum_df,
  "tbb-4_clean_call"
)

ce_tbb4_tree_plot <- add_var_to_tree(
  ce_tree,
  ce_tbb4_anno_df,
  strain_var_colors,
  "tbb4_clean_call",
  xpos = -0.05,
  ypos = -0.09
)

# save the tree
save_tree(
  ce_tbb4_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ce_tbb4_tree_plot,
  2.5,
  4
)

cb_tbb4_anno_df <- create_tree_anno_df(
  cb_iso_var_sum_df,
  "tbb-4_clean_call"
)

cb_tbb4_tree_plot <- add_var_to_tree(
  cb_tree,
  cb_tbb4_anno_df,
  strain_var_colors,
  "tbb4_clean_call",
  xpos = -0.01,
  ypos = -0.135
)

# save the tree
save_tree(
  cb_tbb4_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$cb_tbb4_tree_plot,
  2.5,
  4
)

ct_tbb4_anno_df <- create_tree_anno_df(
  ct_iso_var_sum_df,
  "tbb-4_clean_call"
)

ct_tbb4_tree_plot <- add_var_to_tree(
  ct_tree,
  ct_tbb4_anno_df,
  strain_var_colors,
  "tbb4_clean_call",
  xpos = -0.01,
  ypos = -0.125
)

# save the tree
save_tree(
  ct_tbb4_tree_plot + theme(legend.position = "none"),
  sup_figure_fns$ct_tbb4_tree_plot,
  2.5,
  4
)


# head(ce_mec7_anno_df)

# #### Wrapper function for generating tubulin trees ####

# generate_tubulin_tree <- function(species, tubulin_id) {
#   # Define species parameters
#   if (species == "ce") {
#     tree <- ce_tree
#     iso_var_sum_df <- ce_iso_var_sum_df
#     xpos <- -0.05
#     ypos <- -0.09
#   } else if (species == "cb") {
#     tree <- cb_tree
#     iso_var_sum_df <- cb_iso_var_sum_df
#     xpos <- -0.01
#     ypos <- -0.135
#   } else if (species == "ct") {
#     tree <- ct_tree
#     iso_var_sum_df <- ct_iso_var_sum_df
#     xpos <- -0.01
#     ypos <- -0.125
#   } else {
#     stop("Species must be 'ce', 'cb', or 'ct'")
#   }

#   # Create file paths variable name
#   fp_var_name <- paste0(species, "_", gsub("-", "", tubulin_id), "_fp")
#   fp <- get(fp_var_name)

#   # Create annotation dataframe
#   anno_df <- create_tree_anno_df(
#     iso_var_sum_df,
#     paste0(tubulin_id, "_clean_call")
#   )

#   # Add variant data to tree
#   tree_plot <- add_var_to_tree(
#     tree,
#     anno_df,
#     strain_var_colors,
#     paste0(gsub("-", "", tubulin_id), "_clean_call"),
#     xpos = xpos,
#     ypos = ypos
#   )

#   # Save the tree
#   save_tree(
#     tree_plot + theme(legend.position = "none"),
#     fp,
#     2.5,
#     4
#   )
# }
# #### tbb-1 trees ####
# generate_tubulin_tree("ce", "tbb-1")
# generate_tubulin_tree("cb", "tbb-1")
# generate_tubulin_tree("ct", "tbb-1")

# #### tbb-2 trees ####
# generate_tubulin_tree("ce", "tbb-2")
# generate_tubulin_tree("cb", "tbb-2")
# generate_tubulin_tree("ct", "tbb-2")

# #### mec-7 trees ####

# ## Generate mec-7 trees for all species ##
# generate_tubulin_tree("ce", "mec-7")
# generate_tubulin_tree("cb", "mec-7")
# generate_tubulin_tree("ct", "mec-7")

# #### tbb-4 trees ####

# ## Generate tbb-4 trees for all species ##
# generate_tubulin_tree("ce", "tbb-4")
# generate_tubulin_tree("cb", "tbb-4")
# generate_tubulin_tree("ct", "tbb-4")

#### Plot HTA strain tree for CB and CT ####
# Define function to flag HTA strains and add
# a variant status column to the df for plotting
add_hta_var_stat <- function(iso_sum_df, hta_strains) {
  hta_iso_sum_df <- iso_sum_df %>%
    dplyr::mutate(
      is_hta_strain = ifelse(strain %in% hta_strains, TRUE, FALSE),
      has_var = ifelse(
        `ben-1_clean_call` != "No variant" |
          `tbb-1_clean_call` != "No variant" |
          `tbb-2_clean_call` != "No variant" |
          `mec-7_clean_call` != "No variant" |
          `tbb-4_clean_call` != "No variant",
        "True",
        "False"
      )
    )
  return(hta_iso_sum_df)
}

# Define function to create the tree
# Slightly modified from the previous function
# Adds strain names as labels
plot_hta_strains_tree <- function(tree, tree_anno_df, col_scale, input_col, xpos = NULL, ypos = NULL) {
  base_tree <- ggtree(tree, layout = "equal_angle")

  # add data to base tree
  base_tree_anno <- base_tree %<+% tree_anno_df +
    ggtree::geom_point(
      data = ggtree::td_filter(isTip & is_hta_strain == TRUE),
      ggplot2::aes_string(fill = input_col),
      size = 2,
      shape = 21
    ) +
    ggplot2::scale_fill_manual(
      values = col_scale,
      name = "Beta-tubulin variant",
      guide = ggplot2::guide_legend(override.aes = list(shape = 21))
    ) +
    ggplot2::scale_color_manual(
      values = col_scale,
      name = "Beta-tubulin variant",
      guide = "none"
    ) +
    ggtree::geom_treescale(
      offset = 0.005,
      x = xpos,
      y = ypos,
      linesize = 0.25,
      fontsize = 5,
      family = "Arial"
    ) # Adjust offset, x, and y as needed
  return(base_tree_anno)
}

# Load the HTA strain data
hta_strains_file <-
  "bin/hta_strains.tsv"

hta_strains <- data.table::fread(hta_strains_file)

# Create a list of HTA Cb and CT strains
cb_hta_strains <- hta_strains %>%
  dplyr::filter(species == "c_briggsae") %>%
  dplyr::pull(strain)

ct_hta_strains <- hta_strains %>%
  dplyr::filter(species == "c_tropicalis") %>%
  dplyr::pull(strain)

# use the function to filter the isotype summary and add the variant status
cb_hta_iso_sum_df <- add_hta_var_stat(cb_iso_var_sum_df, cb_hta_strains)

ct_hta_iso_sum_df <- add_hta_var_stat(ct_iso_var_sum_df, ct_hta_strains)

# Create the tree annotation dataframes
cb_hta_anno_df <- cb_hta_iso_sum_df %>%
  dplyr::select(label = strain, is_hta_strain, has_var)

ct_hta_anno_df <- ct_hta_iso_sum_df %>%
  dplyr::select(label = strain, is_hta_strain, has_var)

# Plot the trees
cb_hta_tree_plot <- plot_hta_strains_tree(
  tree = cb_tree,
  tree_anno_df = cb_hta_anno_df,
  col_scale = c(
    "False" = "lightgrey",
    "True" = "red"
  ),
  input_col = "has_var",
  xpos = -0.01,
  ypos = -0.135
)


ct_hta_tree_plot <- plot_hta_strains_tree(
  tree = ct_tree,
  tree_anno_df = ct_hta_anno_df,
  col_scale = c(
    "False" = "lightgrey",
    "True" = "red"
  ),
  input_col = "has_var",
  xpos = -0.01,
  ypos = -0.125
)

# Save the trees
save_tree(
  cb_hta_tree_plot +
    theme(
      legend.position = "top",
      legend.text = ggplot2::element_text(
        size = 11,
        family = "Arial"
      ),
      legend.title = ggplot2::element_text(face = "bold")
    ),
  sup_figure_fns$cb_hta_tree_plot,
  w_in = 7.5,
  h_in = 4
)

save_tree(
  ct_hta_tree_plot +
    theme(
      legend.position = "top",
      legend.text = ggplot2::element_text(
        size = 11,
        family = "Arial"
      ),
      legend.title = ggplot2::element_text(face = "bold")
    ),
  sup_figure_fns$ct_hta_tree_plot,
  w_in = 7.5,
  h_in = 4
)
