library(tidyverse)
library(ggpubr)
library(patchwork)
set.seed(123)

# Load the variant color scale
source("bin/var_color_scale.R")
source("bin/get_res_threshold.R")

# Load standardized output functions and figure mappings
source("bin/outs.R")
source("bin/materials_key.R")

source("scripts/expression_variation/exp_plots.R")

### Inputs ###
isotype_folder_id <- "20250128"
ref_threshold <- 0.75



#### Load data ####

print("loading data")

# define path with folder ID
isotype_var_file <-
  glue::glue(
    "data/isotype_variant_table/c_elegans/{isotype_folder_id}/isotype_variant_summary.tsv"
  )

isotype_summary <- data.table::fread(isotype_var_file)

print("data loaded")

#### Fix Discordant strains ####

print("Fixing discordant strains")

# If a strain has a discordant variant call for ben-1,
# adjust `ben-1_var` to correct id
isotype_summary <- isotype_summary %>%
  dplyr::mutate(
    `ben-1_clean_call` = case_when(
      strain == "JU3125" ~ "Transposon insertion",
      strain == "ED3011" ~ "Splice Donor",
      TRUE ~ `ben-1_clean_call`
    )
  )


print("Discordant strains fixed")

#### Add expression data and variant data to Normalized phenos ####

print("formatting data for expression analysis")

exp_summary <- isotype_summary %>%
  dplyr::select(
    strain,
    abz_hta_norm_pheno,
    abz_hta_norm_res,
    `tbb-1_exp`,
    `tbb-2_exp`,
    `ben-1_exp`,
    `mec-7_exp`,
    `tbb-4_exp`,
    `ben-1_clean_call`
  ) %>%
  # remove any strains that do not have expression data
  dplyr::filter(
    !is.na(`tbb-1_exp`) |
      !is.na(`tbb-2_exp`) |
      !is.na(`ben-1_exp`) |
      !is.na(`mec-7_exp`) |
      !is.na(`tbb-4_exp`)
  ) %>%
  # column to indicate ben-1 variant status
  dplyr::mutate(
    ben1_var_stat = case_when(
      `ben-1_clean_call` == "No variant" ~ FALSE,
      TRUE ~ TRUE
    )
  )

print("data formatted")

# Check how many strains with exp. data have a prev. phenotype valye
# `abz_hta_norm_pheno` is not NA
n_strains_exp_phenoed <- exp_summary %>%
  dplyr::filter(!is.na(abz_hta_norm_pheno)) %>%
  dplyr::pull(strain) %>%
  length()

print(
  glue::glue(
    "Number of strains with expression data and a previous phenotype value: {n_strains_exp_phenoed}"
  )
)

# Check number of strains with low ben-1 expression
n_strains_low_exp <- exp_summary %>%
  dplyr::filter(`ben-1_exp` < 3.75) %>%
  dplyr::pull(strain) %>%
  length()

print(
  glue::glue(
    "Number of strains with low ben-1 expression: {n_strains_low_exp}"
  )
)



# 23 but 4 have low-exp calls

#### Define the resistance thresholds ####

print("Calculating resistance thresholds")

# filter to just previously phenotyped strains
all_phenotyped_iso_var_summary <- isotype_summary %>%
  dplyr::filter(abz_hta_2018 == TRUE | abz_hta_2024 == TRUE)


ce_res <- get_res_strains_ref(
  exp_summary_df = all_phenotyped_iso_var_summary,
  pheno_col = "abz_hta_norm_pheno",
  ref_strain = "N2",
  threshold_per = ref_threshold
)


all_phenotyped_iso_threshold <- ce_res$threshold

print(
  glue::glue(
    "Resistance threshold for all phenotyped strains: {all_phenotyped_iso_threshold}"
  )
)

#### Add catagory groups to ben-1 variants ####

print("Adding ben-1 variant categories")

ben1_meta <- add_ben1_metadata(
  exp_summary
)

print("Ben-1 variant categories added")

#### Figure S4a ben-1 expr. x BZ response with variant categories ####

print("Creating figure ben-1 exp x BZ response with variant categories")

# use the function to generate the plot and get the p-value and r-squared
ben1_bz_var_cat_exp <- create_expression_scatter_plot(
  exp_data = ben1_meta,
  x_column_id = "ben-1_exp",
  y_column_id = "abz_hta_norm_pheno",
  fill_column_id = "ben1_var_cat_meta",
  fill_scale = meta_cat_cols,
  x_label = expression(bolditalic("ben-1") * bold(" expression (TPM)")),
  y_label = expression(bold("Normalized ABZ Response")),
  fill_label = expression(bold("BEN-1 Variation")),
  fill_labels = custom_meta_cat_labels
)


ben1_bz_var_cat_exp_plot <- ben1_bz_var_cat_exp$plot


# print p-value and r-squared
p <- ben1_bz_var_cat_exp$p_value
print(
  glue::glue(
    "P-value for ben-1 expression and BZ response: {p}"
  )
)

r <- ben1_bz_var_cat_exp$r_squared
print(
  glue::glue(
    "R-squared for ben-1 expression and BZ response: {ben1_bz_var_cat_exp_plot$r_squared}"
  )
)

print("Figure ben-1 exp x BZ response with variant categories created")

#### Figure S4b ben-1 var x ben-1 exp ####

print("Creating figure ben-1 variant catagory x ben-1 exp")

# ben1_var_dat_no_low <- ben1_meta %>%
#   dplyr::filter("ben1_var_cat_meta" != "Low ben-1 expression")

# unique(ben1_var_dat_no_low$ben1_var_cat_meta)

# # group by the ben-1 variant category and the number of strains in each category
# ben1_meta %>%
#   dplyr::group_by(ben1_var_cat_meta) %>%
#   dplyr::summarize(
#     n = n()
#   )

# perfrom a wilcox test to see if the expression of ben-1 is different between the variant categories
ben1_var_exp_wilcox <- rstatix::wilcox_test(
  data = ben1_meta,
  formula = `ben-1_exp` ~ ben1_var_cat_meta,
  paired = FALSE
)

# clean up the wilcox test results to add to the plot
ben1_var_exp_wilcox_df <- ben1_var_exp_wilcox %>%
  rstatix::add_significance() %>%
  dplyr::rename(
    p.signif = p.adj.signif
  ) %>%
  # remove ns comparisons
  dplyr::filter(
    p.signif != "ns"
  ) %>%
  rstatix::add_xy_position(
    x = "ben1_var_cat_meta",
    step.increase = 0.03
  )


# get list of significant comparisons from the wilcox test
g1 <- ben1_var_exp_wilcox_df %>%
  dplyr::pull(group1)
g2 <- ben1_var_exp_wilcox_df %>%
  dplyr::pull(group2)

# create a list of significant comparisons
sig_comparisons <- list()
for (i in 1:length(g1)) {
  sig_comparisons[[i]] <- c(g1[i], g2[i])
}



ben1_exp_var_cat_boxplot <- ggplot2::ggplot(
  data = ben1_meta,
  ggplot2::aes(
    x = ben1_var_cat_meta,
    y = `ben-1_exp`,
    fill = ben1_var_cat_meta
  )
) +
  ggplot2::geom_boxplot(
    outliers = FALSE
  ) +
  # add significant comparisons
  ggsignif::geom_signif(
    comparisons = sig_comparisons,
    y_position = ben1_var_exp_wilcox_df$y.position,
    margin_top = 0.01,
    map_signif_level = TRUE,
    annotations = ben1_var_exp_wilcox_df$p.signif
  ) +
  ggplot2::geom_jitter(
    width = 0.2,
    height = 0,
    size = 2,
    alpha = 0.8,
    shape = 21,
    fill = "grey",
    color = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = meta_cat_cols,
    labels = custom_meta_cat_labels
  ) +
  ggplot2::labs(
    x = expression(bolditalic("ben-1") * bold(" consequence")),
    y = expression(bolditalic("ben-1") * bold(" expression (TPM)"))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )



#### plot figure all other tubulin exp x ABZ response ####

# ## Create expression scatter plots ##
tbb1_bz_var_cat_exp_out <-
  create_expression_scatter_plot(
    exp_data = ben1_meta,
    x_column_id = "tbb-1_exp",
    y_column_id = "abz_hta_norm_pheno",
    fill_column_id = "ben1_var_cat_meta",
    fill_scale = meta_cat_cols,
    x_label = expression(bolditalic("tbb-1") * bold(" expression (TPM)")),
    y_label = "Normalized ABZ Response",
    fill_label = expression(bold("BEN-1 Variation")),
    fill_labels = custom_meta_cat_labels
  )

tbb2_bz_var_cat_exp_out <-
  create_expression_scatter_plot(
    exp_data = ben1_meta,
    x_column_id = "tbb-2_exp",
    y_column_id = "abz_hta_norm_pheno",
    fill_column_id = "ben1_var_cat_meta",
    fill_scale = meta_cat_cols,
    x_label = expression(bolditalic("tbb-2") * bold(" expression (TPM)")),
    y_label = "Normalized ABZ Response",
    fill_label = expression(bold("BEN-1 Variation")),
    fill_labels = custom_meta_cat_labels
  )

mec7_bz_var_cat_exp_out <-
  create_expression_scatter_plot(
    exp_data = ben1_meta,
    x_column_id = "mec-7_exp",
    y_column_id = "abz_hta_norm_pheno",
    fill_column_id = "ben1_var_cat_meta",
    fill_scale = meta_cat_cols,
    x_label = expression(bolditalic("mec-7") * bold(" expression (TPM)")),
    y_label = "Normalized ABZ Response",
    fill_label = expression(bold("BEN-1 Variation")),
    fill_labels = custom_meta_cat_labels
  )

# create tbb-4 expression scatter plot with variant categories
tbb4_bz_var_cat_exp_out <-
  create_expression_scatter_plot(
    exp_data = ben1_meta,
    x_column_id = "tbb-4_exp",
    y_column_id = "abz_hta_norm_pheno",
    fill_column_id = "ben1_var_cat_meta",
    fill_scale = meta_cat_cols,
    x_label = expression(bolditalic("tbb-4") * bold(" expression (TPM)")),
    y_label = "Normalized ABZ Response",
    fill_label = expression(bold("BEN-1 Variation")),
    fill_labels = custom_meta_cat_labels
  )

# Function to prepare individual plots (remove legend and adjust theme)
prepare_plot <- function(plot, remove_y_title = FALSE) {
  p <- plot +
    theme(
      legend.position = "none",
      axis.title.y = if (remove_y_title) {
        element_blank()
      } else {
        element_text(
          size = 11,
          face = "bold",
          color = "black",
          family = "Arial"
        )
      },
      axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
      axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", family = "Arial"),
      plot.margin = margin(5, 5, 5, 5)
    )
  return(p)
}

# plot with patchwork to get legend
p1_with_legend <- tbb1_bz_var_cat_exp_out$plot +
  theme(
    legend.position = "top",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
    axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", family = "Arial"),
    plot.margin = margin(5, 5, 5, 5),
    legend.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 11, face = "bold", family = "Arial")
  )

# Remove legends and y-axis titles from other plots
p2_no_legend <- prepare_plot(tbb2_bz_var_cat_exp_out$plot, remove_y_title = TRUE)
p3_no_legend <- prepare_plot(mec7_bz_var_cat_exp_out$plot, remove_y_title = TRUE)
p4_no_legend <- prepare_plot(tbb4_bz_var_cat_exp_out$plot, remove_y_title = TRUE)

# Combine with patchwork
combined_patchwork <- (p1_with_legend / p2_no_legend / p3_no_legend / p4_no_legend) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 12, family = "Helvetica")
  )

# Add common y-axis label
other_beta_tubulin_exp_plot <- cowplot::plot_grid(
  grid::textGrob(
    "Normalized ABZ Response",
    gp = grid::gpar(fontsize = 11, fontface = "bold", fontfamily = "Arial"),
    rot = 90
  ),
  combined_patchwork,
  ncol = 2,
  rel_widths = c(0.03, 0.97)
)


save_plot(
  tplot = other_beta_tubulin_exp_plot,
  fn_list = sup_figure_fns$other_beta_tubulin_exp_plot,
  w_in = 7.5,
  h_in = 10
)



# # create tbb-1 expression scatter plot with variant categories
# tbb1_bz_var_cat_exp_out <-
#   create_expression_scatter_plot(
#     exp_data = ben1_meta,
#     x_column_id = "tbb-1_exp",
#     y_column_id = "abz_hta_norm_pheno",
#     fill_column_id = "ben1_var_cat_meta",
#     fill_scale = meta_cat_cols,
#     x_label = expression(bolditalic("tbb-1") * bold(" expression (TPM)")),
#     y_label = "Normalized ABZ Response",
#     fill_label = expression(bold("BEN-1 Variation")),
#     res_threshold = all_phenotyped_iso_threshold,
#     fill_labels = custom_meta_cat_labels
#   )

# # create tbb-2 expression scatter plot with variant categories
# tbb2_bz_var_cat_exp_out <-
#   create_expression_scatter_plot(
#     exp_data = ben1_meta,
#     x_column_id = "tbb-2_exp",
#     y_column_id = "abz_hta_norm_pheno",
#     fill_column_id = "ben1_var_cat_meta",
#     fill_scale = meta_cat_cols,
#     x_label = expression(bolditalic("tbb-2") * bold(" expression (TPM)")),
#     y_label = "Normalized ABZ Response",
#     fill_label = expression(bold("BEN-1 Variation")),
#     res_threshold = all_phenotyped_iso_threshold,
#     fill_labels = custom_meta_cat_labels
#   )

# # create mec-7 expression scatter plot with variant categories
# mec7_bz_var_cat_exp_out <-
#   create_expression_scatter_plot(
#     exp_data = ben1_meta,
#     x_column_id = "mec-7_exp",
#     y_column_id = "abz_hta_norm_pheno",
#     fill_column_id = "ben1_var_cat_meta",
#     fill_scale = meta_cat_cols,
#     x_label = expression(bolditalic("mec-7") * bold(" expression (TPM)")),
#     y_label = "Normalized ABZ Response",
#     fill_label = expression(bold("BEN-1 Variation")),
#     res_threshold = all_phenotyped_iso_threshold,
#     fill_labels = custom_meta_cat_labels
#   )

# # create tbb-4 expression scatter plot with variant categories
# tbb4_bz_var_cat_exp_out <-
#   create_expression_scatter_plot(
#     exp_data = ben1_meta,
#     x_column_id = "tbb-4_exp",
#     y_column_id = "abz_hta_norm_pheno",
#     fill_column_id = "ben1_var_cat_meta",
#     fill_scale = meta_cat_cols,
#     x_label = expression(bolditalic("tbb-4") * bold(" expression (TPM)")),
#     y_label = "Normalized ABZ Response",
#     fill_label = expression(bold("BEN-1 Variation")),
#     res_threshold = all_phenotyped_iso_threshold,
#     fill_labels = custom_meta_cat_labels
#   )

# ## Adjust theme ##
# adjust_theme <- function(plot, top = FALSE, y_title_width = 15) {
#   p <- plot +
#     ggplot2::theme(
#       axis.title.y = element_text(
#         size = 9,
#         face = "bold",
#         color = "black",
#         family = "Arial"
#       ),
#       axis.text.y = element_text(
#         size = 10,
#         family = "Arial",
#         color = "black"
#       ),
#       axis.text.x = element_text(
#         size = 10,
#         family = "Arial",
#         color = "black"
#       ),
#       axis.title.x = element_text(
#         size = 11,
#         face = "bold",
#         family = "Arial"
#       ),
#     )

#   # Wrap y-axis title if it exists
#   if (!is.null(p$labels$y)) {
#     p$labels$y <- stringr::str_wrap(p$labels$y, width = y_title_width)
#   }

#   # if the top argument is TRUE, add a top legend
#   if (top) {
#     p <- p +
#       ggplot2::theme(
#         legend.position = "top"
#       )
#   } else {
#     p <- p +
#       ggplot2::theme(
#         legend.position = "none"
#       )
#   }
#   return(p)
# }

# plot_list <- list(
#   tbb1_bz_var_cat_exp_out$plot,
#   tbb2_bz_var_cat_exp_out$plot,
#   mec7_bz_var_cat_exp_out$plot,
#   tbb4_bz_var_cat_exp_out$plot
# )

# # adjust the theme of each plot with the
# # adjust_theme function, for the first plot
# # in the list set top = TRUE
# adjusted_plot_list <- lapply(
#   seq_along(plot_list),
#   function(i) {
#     adjust_theme(plot_list[[i]], top = (i == 1), y_title_width = 15)
#   }
# )

# # combine
# other_beta_tubulin_expression <- adjusted_plot_list[[1]] / adjusted_plot_list[[2]] / adjusted_plot_list[[3]] / adjusted_plot_list[[4]] +
#   plot_annotation(
#     tag_levels = "a",
#     tag_prefix = "",
#     tag_suffix = ""
#   ) & theme(
#   plot.tag = element_text(face = "bold", size = 12, family = "Helvetica")
# )

# save_plot(
#   tplot = other_beta_tubulin_expression,
#   fn_list = other_beta_tubulin_abz_fn,
#   w_in = 7.5,
#   h_in = 8
# )


#### Save figures ####

## Save ben-1 exp x ABZ response sactter & bp  ###

p1 <- ben1_bz_var_cat_exp_plot +
  ggplot2::labs(caption = NULL) +
  ggplot2::theme(
    legend.text = element_text(
      size = 11,
      family = "Arial"
    ),
    legend.title = element_text(
      size = 11,
      family = "Arial",
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 11,
      face = "bold",
      family = "Arial"
    ), # Bold x-axis label in Arial
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 11,
      family = "Arial",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11,
      family = "Arial",
      face = "bold"
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 2))


# Modify panel 2 for plot
p2 <- ben1_exp_var_cat_boxplot +
  theme(
    text = ggplot2::element_text(size = 10, family = "Helvetica"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = 11,
      family = "Arial",
      color = "black"
    )
  )


# create a combined plot
main_figure <- ggpubr::ggarrange(
  p1,
  p2,
  ncol = 1,
  labels = c("A", "B"),
  font.label = list(
    size = 12,
    color = "black",
    family = "Arial"
  )
)

save_plot(
  tplot = main_figure,
  fn_list = sup_figure_fns$main_figure,
  w_in = 7.5,
  h_in = 7
)
