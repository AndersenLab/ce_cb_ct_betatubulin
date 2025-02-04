library(tidyverse)
library(ggpubr)
set.seed(123)

# Load the variant color scale
source("bin/var_color_scale.R")
source("bin/get_res_threshold.R")

source("scripts/expression_variation/exp_plots.R")


### Inputs ###
isotype_folder_id <- "20250128"
ref_threshold <- 0.75
#### Output ####



figure_out_dir <- "figures/figure_S4"

# if it does exist, delete the contents
if (dir.exists(figure_out_dir)) {
  unlink(figure_out_dir, recursive = TRUE)
}

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
    `ben-1_clean_call`
  ) %>%
  # remove any strains that do not have expression data
  dplyr::filter(
    !is.na(`tbb-1_exp`) | !is.na(`tbb-2_exp`) | !is.na(`ben-1_exp`)
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



#23 but 4 have low-exp calls

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
    fill_label = expression(italic("ben-1") * " consequence"),
    res_threshold = all_phenotyped_iso_threshold
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
    step.increase = 0.01
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
  )+
  ggplot2::geom_boxplot(
    outliers = FALSE
  )+
  # add significant comparisons
  ggsignif::geom_signif(
    comparisons = sig_comparisons,
    map_signif_level = TRUE,
    y_position = ben1_var_exp_wilcox_df$y.position,
    annotations = ben1_var_exp_wilcox_df$p.signif
  ) +
  ggplot2::geom_jitter(
    width = 0.2,
    height = 0,
    size = 2
  ) +
  ggplot2::scale_fill_manual(values = meta_cat_cols) +
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






# ben1_exp_var_cat_boxplot <- create_expression_boxplot(
#   data = ben1_meta  %>% dplyr::filter(strain != "CX11254"),
#   x_col = "ben1_var_cat_meta",
#   y_col = "ben-1_exp",
#   x_label = expression(bolditalic("ben-1") * bold(" consequence")),
#   y_label = expression(bolditalic("ben-1") * bold(" expression (TPM)")),
#   comparisons_list = list(
#     c("No variant", "SV"),
#     c("No variant", "Frame altering"),
#     c("No variant", "Missense"),
#     c("No variant", "Start/Stop altering")
#   ),
#   fill_column = "ben1_var_cat_meta",
#   fill_scale = meta_cat_cols
# )

# ben1_exp_var_cat_boxplot

print("Figure ben-1 var x ben-1 exp created")

#### Combine into manuscript Main figure ####

print("Creating main figure")

# Modify panel 1 for plot
p1 <- ben1_bz_var_cat_exp_plot +
  ggplot2::labs(caption = NULL) +
  ggplot2::theme(
    axis.title.x = element_text(
      size = 10,
      face = "bold",
      family = "Helvetica"
      ),  # Bold x-axis label in Arial
    axis.title.y = element_text(
    size = 10,
    face = "bold",
    family = "Helvetica"
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
    axis.ticks.x = element_blank()
    )


# create a combined plot
main_figure <- ggpubr::ggarrange(
  p1,
  p2,
  ncol = 1,
  labels = c("A", "B"),
  font.label = list(
    size = 10,
    color = "black", 
    family = "Arial"
    )
)

# save the plot
ggsave(
  filename = glue::glue("{figure_out_dir}/figure_S4.png"),
  plot = main_figure,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 300
)


ggsave(
  filename = glue::glue("{figure_out_dir}/figure_S4.eps"),
  plot = main_figure,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 300
)

#### Combine into supplementary figure ####

# # create a supplementary figure of tbb-1 and tbb-2 results
# supplementary_other_tub_exp <- ggpubr::ggarrange(
#   tbb1_bz_var_cat_exp_plot,
#   tbb2_bz_var_cat_exp_plot,
#   tbb1_exp_var_cat_boxplot +
#     ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#   tbb2_exp_var_cat_boxplot +
#     ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#   common.legend = TRUE,
#   labels = c("a", "b", "c", "d")
# )

# # save the plot
# ggplot2::ggsave(
#   filename = glue::glue("{figure_out_dir}/supplementary_figure_tbb1_tbb2.jpg"),
#   plot = supplementary_other_tub_exp,
#   width = 7.5,
#   height = 10,
#   units = "in",
#   dpi = 300
# )

# # create supplementary figure of scatter plots x assay
# supplementary_figure <- ggpubr::ggarrange(
#   ben1_bz_var_cat_exp_assay_plot,
#   tbb1_bz_var_cat_exp_assay_plot,
#   tbb2_bz_var_cat_exp_assay_plot,
#   common.legend = TRUE,
#   labels = c("a", "b", "c")
# )

# # save the plot
# ggplot2::ggsave(
#   filename = glue::glue("{figure_out_dir}/supplementary_figure_scatter_x_assay.jpg"),
#   plot = supplementary_figure,
#   width = 7.5,
#   height = 10,
#   units = "in",
#   dpi = 300
# )


#### Archive ####
## Linear Model if ben-1 variant category is associated with ben-1 expression ##

# # create a linear model to test if ben-1 variant meta category is associated with ben-1 expression
# ben1_var_cat_exp_lm <- lm(`ben-1_exp` ~ ben1_var_cat_meta - 1, data = exp_summary)

# ## plot relationship as coefficient plot ###
# ## create a coefficient plot
# ## plot the coefficients as points with error bars

# # extract the coefficients and standard errors + clean up the data
# coef_df <- broom::tidy(ben1_var_cat_exp_lm) %>%
#   dplyr::filter(term != "(Intercept)") %>%
#   dplyr::mutate(
#     ben1_var_cat = stringr::str_replace(
#       term,
#       pattern = "^ben1_var_cat_meta",
#       replacement = ""
#     ),
#     log10_p_value = -log10(p.value)
#   )
# # # create the ben1_var_cat_meta column based on the term column
# # coef_df <- coef_df %>%
# #   mutate(
# #     ben1_var_cat_meta = case_when(
# #       str_detect(term, "SV") ~ "SV",
# #       str_detect(term, "Frame Altering") ~ "Frame Altering",
# #       str_detect(term, "Missense") ~ "Missense",
# #       str_detect(term, "Start/Stop Codon") ~ "Start/Stop Codon",
# #       TRUE ~ "No variant"
# #     )
# #   )

# # create a plot
# coef_plot <- ggplot(
#   coef_df,
#   aes(
#     x = ben1_var_cat,
#     y = estimate,
#     ymin = estimate - std.error,
#     ymax = estimate + std.error,
#     color = p.value,
#     shape = p.value < 0.05
#   )
# ) +
#   geom_pointrange() +
#   labs(
#     x = "ben-1 variant meta category",
#     y = "Coefficient",
#     title = "Association between ben-1 variant meta category and ben-1 expression"
#   ) +
#   theme_minimal() +
#   ggplot2::coord_flip() +
#   # Add color scale for significance
#   scale_color_gradient(low = "red", high = "blue", limits = c(0, 1))

# # save the plot
# ggsave(
#   filename = glue::glue("{figure_out_dir}/ben1_var_cat_x_ben1_exp_coef.jpg"),
#   plot = coef_plot,
#   width = 7.5,
#   height = 7.5,
#   units = "in",
#   dpi = 300
# )

# ### Test plotting volcano plot where x-axis is the coefficient and y-axis is the -log10 p-value
# volcano_plot <- ggplot(
#   coef_df,
#   aes(
#     x = estimate,
#     # VARIABLE IS MISSING
#     y = log_p_value,
#     color = term
#   )
# ) +
#   geom_point() +
#   labs(
#     x = "Coefficient",
#     y = "-log10 p-value",
#     title = "Volcano plot for ben-1 variant meta category association with ben-1 expression"
#   ) +
#   theme_minimal()

# # save the plot
# ggsave(
#   filename = glue::glue("{figure_out_dir}/ben1_var_cat_x_ben1_exp_volcano.jpg"),
#   plot = volcano_plot,
#   width = 7.5,
#   height = 7.5,
#   units = "in",
#   dpi = 300
# )

# ## Test if ben-1 variant category is associated with BZ response ##
# ben1_var_cat_bz_lm <- lm(strain_norm_abz_response ~ ben1_var_cat_meta, data = exp_summary)

# # summary of the linear model
# summary(ben1_var_cat_bz_lm)

# ## Test if ben-1 variant category is associated with tbb-2 expression ##
# ben1_var_cat_tbb2_lm <- lm(`tbb-2_exp` ~ ben1_var_cat_meta, data = exp_summary)

# # summary of the linear model
# summary(ben1_var_cat_tbb2_lm)

# ## Test if ben-1 variant category is associated with tbb-1 expression ##
# ben1_var_cat_tbb1_lm <- lm(`tbb-1_exp` ~ ben1_var_cat_meta, data = exp_summary)

# # summary of the linear model
# summary(ben1_var_cat_tbb1_lm)
