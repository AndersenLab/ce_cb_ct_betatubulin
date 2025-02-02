# Code to generate trees for all species denoting high-impact beta-tubulin variants
set.seed(123)

#### Load required libraries ####
library(ggtree)


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

## Figure 5 Trees ##

# Ce ben-1 high impact tree
figure_5b_out_paths <- c(
  png = "figures/figure_5/figure_5a.png",
  eps = "figures/figure_5/figure_5a.eps"
)

# Cb ben-1 high impact tree
figure_5c_out_paths <- c(
  png = "figures/figure_5/figure_5b.png",
  eps = "figures/figure_5/figure_5b.eps"
)

# Ct ben-1 high impact tree
figure_5d_out_paths <- c(
  png = "figures/figure_5/figure_5c.png",
  eps = "figures/figure_5/figure_5c.eps"
)

## Figure S10 Trees ##

# Ce tbb-1 high impact tree
figure_S10b_out_paths <- c(
  png = "figures/figure_S10/figure_S10a.png",
  eps = "figures/figure_S10/figure_S10a.eps"
)

# Cb tbb-1 high impact tree
figure_S10c_out_paths <- c(
  png = "figures/figure_S10/figure_S10b.png",
  eps = "figures/figure_S10/figure_S10b.eps"
)

# Ct tbb-1 high impact tree
figure_S10d_out_paths <- c(
  png = "figures/figure_S10/figure_S10c.png",
  eps = "figures/figure_S10/figure_S10c.eps"
)

## Figure S11 Trees ##
# Ce tbb-2 high impact tree
figure_S11b_out_paths <- c(
  png = "figures/figure_S11/figure_S11a.png",
  eps = "figures/figure_S11/figure_S11a.eps"
)

# Cb tbb-2 high impact tree
figure_S11c_out_paths <- c(
  png = "figures/figure_S11/figure_S11b.png",
  eps = "figures/figure_S11/figure_S11b.eps"
)

# Ct tbb-2 high impact tree
figure_S11d_out_paths <- c(
  png = "figures/figure_S11/figure_S11c.png",
  eps = "figures/figure_S11/figure_S11c.eps"
)

## Figure S1 CB HTA strains ##

figure_S1_paths <- c(
  png = "figures/figure_S1/figure_S1.png",
  eps = "figures/figure_S1/figure_S1.eps"
)

## Figure S2 CT HTA strains ##
figure_S2_paths <- c(
  png = "figures/figure_S2/figure_S2.png",
  eps = "figures/figure_S2/figure_S2.eps"
)

# # if it does exist, delete the contents
# if (dir.exists(data_out_dir)) {
#   unlink(data_out_dir, recursive = TRUE)
# }
# # if it doesn't exist, create it
# dir.create(data_out_dir, recursive = TRUE)

# figure_out_dir <- glue::glue(
#   "figures/phylogeographic_distribution/{isotype_folder_id}"
# )

# # create the output directory if it doesn't exist
# if (!dir.exists(figure_out_dir)) {
#   dir.create(figure_out_dir, recursive = TRUE)
# }


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
  base_tree <- ggtree(tree, layout = "equal_angle")

  # add data to base tree
  base_tree_anno <- base_tree %<+% tree_anno_df +
    ggtree::geom_point(
      data = ggtree::td_filter(isTip & !!rlang::sym(input_col) != "No variant"),
      ggplot2::aes_string(fill = input_col),
      shape = 21,
      size = 2
    ) +
    ggplot2::scale_fill_manual(values = col_scale) +
    ggtree::geom_treescale(
      offset = 0.001,
      x = xpos,
      y = ypos,
      fontsize = 5,
      family = "Arial"
    ) # Adjust offset, x, and y as needed
  return(base_tree_anno)
}

#### Save the trees ####
save_tree <- function(tree_plot, fn_list, w_in, h_in) {
  # get the folder name from the first file name
  folder <- dirname(fn_list[1])
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }

  fn_png <- fn_list[1]
  fn_eps <- fn_list[2]


  # save eps plot
  ggplot2::ggsave(
    filename = fn_eps,
    plot = tree_plot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 300
  )
  # save png plot
  ggplot2::ggsave(
    filename = fn_png,
    plot = tree_plot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 300
  )
  # Check if either fn_eps or fn_png exist
  # if so delete
}

# # Function to process all genes
# process_gene <- function(
#   gene_id, 
#   ce_iso_var_sum_df, 
#   cb_iso_var_sum_df, 
#   ct_iso_var_sum_df, 
#   ce_tree, 
#   cb_tree, 
#   ct_tree, 
#   figure_out_dir
# ) {
#   clean_call_col <- paste0(gsub("-", "", gene_id), "_clean_call")
#   ce_anno_df <- create_tree_anno_df(ce_iso_var_sum_df, paste0(gene_id, "_clean_call"))
#   cb_anno_df <- create_tree_anno_df(cb_iso_var_sum_df, paste0(gene_id, "_clean_call"))
#   ct_anno_df <- create_tree_anno_df(ct_iso_var_sum_df, paste0(gene_id, "_clean_call"))

#   ce_tree_plot <- add_var_to_tree(
#     ce_tree, 
#     ce_anno_df, 
#     strain_var_colors, 
#     clean_call_col, 
#     xpos = -0.05, 
#     ypos = -0.09
#     )
#   cb_tree_plot <- add_var_to_tree(cb_tree, cb_anno_df, strain_var_colors, clean_call_col, xpos = -0.01, ypos = -0.135)
#   ct_tree_plot <- add_var_to_tree(ct_tree, ct_anno_df, strain_var_colors, clean_call_col, xpos = -0.01, ypos = -0.125)

#   save_tree(
#     ce_tree_plot + theme(legend.position = "none"),
#     gene_id,
#     "c_elegans",
#     f_name = "var_tree",
#     figure_out_dir,
#     w_in = 2.5,
#     h_in = 4
#     )
   
#   save_tree(
#     cb_tree_plot + theme(legend.position = "none"),
#     gene_id,
#     "c_briggsae",
#     f_name = "var_tree",
#     figure_out_dir,
#     w_in = 2.5,
#     h_in = 4
#     )
  
#   save_tree(
#     ct_tree_plot + theme(legend.position = "none"),
#     gene_id,
#     "c_tropicalis",
#     f_name = "var_tree",
#     figure_out_dir,
#     w_in = 2.5,
#     h_in = 4
#     )
# }

# # Process all genes
# genes <- c("ben-1", "tbb-1", "tbb-2")
# for (gene in genes) {
#   process_gene(gene, ce_iso_var_sum_df, cb_iso_var_sum_df, ct_iso_var_sum_df, ce_tree, cb_tree, ct_tree, figure_out_dir)
# }

#### Figure 5 Trees ####

## Ce ben-1 high impact tree 5a ##

#create data frame for ben-1
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
  figure_5b_out_paths,
  2.5,
  4
  )

## Cb ben-1 high impact tree 5b ##

#create data frame for ben-1
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
  figure_5c_out_paths,
  2.5,
  4
  )

## Ct ben-1 high impact tree 5c ##
#create data frame for ben-1
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
  figure_5d_out_paths,
  2.5,
  4
  )

#### Figure S10 Trees ####

## Ce tbb-2 high impact tree S10a ##
#create data frame for tbb-2
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
  figure_S10b_out_paths,
  2.5,
  4
  )

## Cb tbb-2 high impact tree S10b ##
#create data frame for tbb-2
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
  figure_S10c_out_paths,
  2.5,
  4
  )

## Ct tbb-2 high impact tree S10c ##
#create data frame for tbb-2
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
  figure_S10d_out_paths,
  2.5,
  4
  )

#### Figure S11 Trees ####

## Ce tbb-2 high impact tree S11a ##
#create data frame for tbb-2
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
  figure_S11b_out_paths,
  2.5,
  4
  )

## Cb tbb-2 high impact tree S11b ##
#create data frame for tbb-2

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
  figure_S11c_out_paths,
  2.5,
  4
  )

## Ct tbb-2 high impact tree S11c ##
#create data frame for tbb-2
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
  figure_S11d_out_paths,
  2.5,
  4
  )


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
        `tbb-2_clean_call` != "No variant",
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
      name = "Beta-tubulin Variant",
      guide = ggplot2::guide_legend(override.aes = list(shape = 21))
      ) +
    ggplot2::scale_color_manual(
      values = col_scale,
      name = "Beta-tubulin Variant",
      guide = "none"
      ) +
    ggtree::geom_treescale(
      offset = 0.01,
      x = xpos,
      y = ypos,
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
    "False" = "#747070",
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
    "False" = "#747070",
    "True" = "red"
    ), 
  input_col = "has_var", 
  xpos = -0.01,
  ypos = -0.125
  )

# Save the trees
save_tree(
  cb_hta_tree_plot + theme(legend.position = "top"),
  figure_S1_paths,
  w_in = 7.5,
  h_in = 4
  )

save_tree(
  ct_hta_tree_plot + theme(legend.position = "top"),
  figure_S2_paths,
  w_in = 7.5,
  h_in = 4
  )
