library(dplyr)
library(data.table)
library(glue)
library(ggplot2)
library(broom)

#### Functions ####
# Function to create scatter plots with R2 and p-value
create_scatter_plot <- function(df, x_var, y_var, x_label, y_label) {
  plot <- ggplot2::ggplot(df, aes_string(x = x_var, y = y_var, color = "tubulin_id")) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "grey",
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "tbb1" = "blue",
        "tbb2" = "green",
        "ben1" = "purple"
      ),
      labels = c(
        "tbb1" = expression(italic("tbb-1")),
        "tbb2" = expression(italic("tbb-2")),
        "ben1" = expression(italic("ben-1"))
      )) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      color = "Beta-tubulin gene"
    ) +
    ggplot2::theme_bw()
  
  # Set x-axis limits based on x_var
  if (x_var == "BLOSUM") {
    plot <- plot + ggplot2::xlim(-4, 4)
  } else if (x_var == "Grantham") {
    plot <- plot + ggplot2::xlim(0, 215)
  } else if (x_var == "Percent_Protein") {
    plot <- plot + ggplot2::xlim(0, 100)
  }
  
  # Calculate R2 and p-value
  lm_model <- lm(as.formula(paste(y_var, "~", x_var)), data = df)
  model_summary <- broom::glance(lm_model)
  r2 <- model_summary$r.squared
  p_value <- model_summary$p.value

  # Add R2 and p-value to the plot
  plot <- plot +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
    label = paste(
      "italic(R)^2 == ", round(r2, 2), 
      "*','~italic(p) == ", format.pval(p_value, digits = 2)
    ),      hjust = 1.1,
      vjust = 1.1,
      parse = TRUE,
      size = 5,
      color = "black",
      family = "Helvetica"
    )
  return(plot)
}


create_faceted_scatter_plot <- function(df) {
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(BLOSUM, Grantham, Percent_Protein),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(Metric, Percent_Protein = "Percent Protein")
      )
  
  plot <- ggplot2::ggplot(df_long, aes(x = Value, y = mean_median_wormlength_um_delta_reg, color = tubulin_id)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "grey",
      alpha = 0.8
    ) +
    ggplot2::facet_wrap(~ Metric, scales = "free_x") +
    ggplot2::scale_color_manual(
      values = c(
        "tbb1" = "blue",
        "tbb2" = "green",
        "ben1" = "purple"
      ),
      labels = c(
        "tbb1" = "TBB-1",
        "tbb2" = "TBB-2",
        "ben1" = "BEN-1"
      )) +
    ggplot2::labs(
      x = "Metric Value",
      y = "Normalized animal length (µm)",
      color = "Beta-tubulin"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = element_text(family = "Helvetica", size = 11),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )
  
  return(plot)
}

# Function to save the plots
save_plot <- function(tplot, fn_list, w_in, h_in) {
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
    plot = tplot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 300
  )
  # save png plot
  ggplot2::ggsave(
    filename = fn_png,
    plot = tplot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 300
  )
}

# Create faceted scatter plot for c_elegans
#ce_faceted_plot <- create_faceted_scatter_plot(ce_df)



#### Inputs ####

# blast key file to filter the annotation file
blast_key_file <- "data/blast_summary/20230324_blast/blast_key.tsv"

# load the blast key
blast_key <- data.table::fread(blast_key_file) %>%
  # remove tbb-6 from the key
  dplyr::filter(ce_common != "tbb_6")

# paths to the annotation flat file
ce_anno_file <- "data/caendr_annotation/c_elegans/20231213/WI.20231213.strain-annotation.tsv"
cb_anno_file <- "data/caendr_annotation/c_briggsae/20240129/WI.20240129.strain-annotation.tsv"
ct_anno_file <- "data/caendr_annotation/c_tropicalis/20231201/WI.20231201.strain-annotation.tsv"

# isotype variant summary data
isotype_folder_id <- "20250128"

ce_var_file <- glue::glue("data/isotype_variant_table/c_elegans/{isotype_folder_id}/isotype_variant_summary.tsv")
cb_var_file <- glue::glue("data/isotype_variant_table/c_briggsae/{isotype_folder_id}/isotype_variant_summary.tsv")
ct_var_file <- glue::glue("data/isotype_variant_table/c_tropicalis/{isotype_folder_id}/isotype_variant_summary.tsv")

## HTA data ##
ce_hta_data_file <- "data/hta_summaries/elegans_df_ABZ_summary_mean_median_wormlength.csv"
cb_hta_data_file <- "data/hta_summaries/briggsae_df_ABZ_summary_mean_median_wormlength.csv"
ct_hta_data_file <- "data/hta_summaries/tropicalis_df_ABZ_summary_mean_median_wormlength.csv"

# resistance thresholds for each species
rest_thres_file <- 'data/hta_res_threshold/res_thresholds.csv'


#### Outputs ####

## Data ##

subscores_out_name <- "tables/table_S7/table_S7.tsv"

table_out_dir <- dirname(subscores_out_name)

# create the output directory if it doesn't exist
if (!dir.exists(table_out_dir)) {
  dir.create(table_out_dir, recursive = TRUE)
}

## Figures ##
ce_missense_fn <- c(
  png = "figures/figure_S10/figure_S10.png",
  eps = "figures/figure_S10/figure_S10.eps"
)

# ce_abz_missense_fn <- c(
#   png = "figures/figure_S11/figure_S11.png",
#   eps = "figures/figure_S11/figure_S11.eps"
# )

cb_missense_fn <- c(
  png = "figures/figure_S11/figure_S11.png",
  eps = "figures/figure_S11/figure_S11.eps"
)

ct_missense_fn <- c(
  png = "figures/figure_S12/figure_S12.png",
  eps = "figures/figure_S12/figure_S12.eps"
)



#### Load the isotype variant data ####

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

#### Load the annotation files ####
ce_anno_df <- data.table::fread(
  ce_anno_file,
  select = c(
    "Strains",
    "AMINO_ACID_CHANGE",
    "BLOSUM",
    "Grantham",
    "Percent_Protein",
    "GENE",
    "CONSEQUENCE",
    "TRANSCRIPT"
  )
)
cb_anno_df <- data.table::fread(
  cb_anno_file,
  select = c(
    "Strains",
    "AMINO_ACID_CHANGE",
    "BLOSUM",
    "Grantham",
    "Percent_Protein",
    "GENE",
    "CONSEQUENCE",
    "TRANSCRIPT"
  )
)
ct_anno_df <- data.table::fread(
  ct_anno_file,
  select = c(
    "Strains",
    "AMINO_ACID_CHANGE",
    "BLOSUM",
    "Grantham",
    "Percent_Protein",
    "GENE",
    "CONSEQUENCE",
    "TRANSCRIPT"
  )
)

#### Load the HTA data ####

ce_hta_data <- data.table::fread(ce_hta_data_file)
cb_hta_data <- data.table::fread(cb_hta_data_file)
ct_hta_data <- data.table::fread(ct_hta_data_file)

# Read in the resistance thresholds
res_thres_df <- data.table::fread(rest_thres_file)

ce_res_thres <- res_thres_df %>%
  dplyr::filter(species == "c_elegans") %>%
  dplyr::pull(hta_threshold)

cb_res_thres <- res_thres_df %>%
  dplyr::filter(species == "c_briggsae") %>%
  dplyr::pull(hta_threshold)

ct_res_thres <- res_thres_df %>%
  dplyr::filter(species == "c_tropicalis") %>%
  dplyr::pull(hta_threshold)

# flag the resistant strains in the HTA data for each species

ce_hta_data <- ce_hta_data %>%
  dplyr::mutate(
    resistant = mean_median_wormlength_um_delta_reg > ce_res_thres
  )

cb_hta_data <- cb_hta_data %>%
  dplyr::mutate(
    resistant = mean_median_wormlength_um_delta_reg > cb_res_thres
  )

ct_hta_data <- ct_hta_data %>%
  dplyr::mutate(
    resistant = mean_median_wormlength_um_reg_delta > ct_res_thres
  )

#### Filter the annotation file to beta-tubulins ####

# create vector of beta-tubulin transcript ids
ce_beta_tubs <- blast_key$N2_transcript_id
cb_beta_tubs <- blast_key$QX1410_transcript_id
ct_beta_tubs <- blast_key$NIC58_transcript_id

# Function to filter isotype variant summary to missense variants
filter_isotype_missense <- function(iso_sum_df, gene_id) {
  clean_call_col <- paste0(gene_id, "_clean_call")
  snv_call_col <- paste0(gene_id, "_high_impact_SNV")
  tubulin_id <- gsub("-", "", gene_id)
  
  iso_df_missense <- iso_sum_df %>%
    dplyr::filter(!!rlang::sym(clean_call_col) == "Missense") %>%
    dplyr::select(
      strain,
      #abz_hta_norm_pheno,
      #abz_hta_norm_res,
      clean_var = !!rlang::sym(clean_call_col),
      snv_call = !!rlang::sym(snv_call_col)
    ) %>%
    dplyr::mutate(tubulin_id = tubulin_id)
  
  return(iso_df_missense)
}

# Update the function to filter the annotation file
filter_blossum_grantham <- function(anno_df, transcript_ids, gene_id) {
  anno_df <- anno_df %>%
    dplyr::filter(TRANSCRIPT %in% transcript_ids) %>%
    dplyr::filter(!is.na(BLOSUM) & !is.na(Grantham)) %>%
    dplyr::filter(GENE == gene_id, CONSEQUENCE == "missense") %>%
    dplyr::select(GENE, AMINO_ACID_CHANGE, BLOSUM, Grantham, Percent_Protein, Strains) %>%
    tidyr::separate_rows(Strains, sep = ",")
  
  return(anno_df)
}

# Function to pull BLOSUM and Grantham scores
pull_blossum_grantham_scores <- function(iso_sum_df, anno_df, transcript_ids, gene_id, species) {
  iso_missense_df <- filter_isotype_missense(iso_sum_df, gene_id)
  anno_missense_df <- filter_blossum_grantham(anno_df, transcript_ids, gene_id)
  
  strain_scores_df <- iso_missense_df %>%
    dplyr::left_join(anno_missense_df, by = c("strain" = "Strains")) %>%
    dplyr::mutate(Species = species)
  
  return(strain_scores_df)
}

# Ensure all dataframes have the same structure
ensure_structure <- function(df) {
  if (nrow(df) == 0) {
    df <- tibble::tibble(
      strain = character(),
      clean_var = character(),
      snv_call = character(),
      tubulin_id = character(),
      GENE = character(),
      AMINO_ACID_CHANGE = character(),
      BLOSUM = numeric(),
      Grantham = numeric(),
      Percent_Protein = numeric(),
      Species = character()
    )
  }
  return(df)
}

# Pull scores for all species and genes
all_beta_tub_blossum_grantham <- dplyr::bind_rows(
  ensure_structure(pull_blossum_grantham_scores(ce_iso_var_sum_df, ce_anno_df, ce_beta_tubs, "ben-1", "c_elegans")),
  ensure_structure(pull_blossum_grantham_scores(cb_iso_var_sum_df, cb_anno_df, cb_beta_tubs, "ben-1", "c_briggsae")),
  ensure_structure(pull_blossum_grantham_scores(ct_iso_var_sum_df, ct_anno_df, ct_beta_tubs, "ben-1", "c_tropicalis")),
  ensure_structure(pull_blossum_grantham_scores(ce_iso_var_sum_df, ce_anno_df, ce_beta_tubs, "tbb-1", "c_elegans")),
  ensure_structure(pull_blossum_grantham_scores(cb_iso_var_sum_df, cb_anno_df, cb_beta_tubs, "tbb-1", "c_briggsae")),
  ensure_structure(pull_blossum_grantham_scores(ct_iso_var_sum_df, ct_anno_df, ct_beta_tubs, "tbb-1", "c_tropicalis")),
  ensure_structure(pull_blossum_grantham_scores(ce_iso_var_sum_df, ce_anno_df, ce_beta_tubs, "tbb-2", "c_elegans")),
  ensure_structure(pull_blossum_grantham_scores(cb_iso_var_sum_df, cb_anno_df, cb_beta_tubs, "tbb-2", "c_briggsae")),
  ensure_structure(pull_blossum_grantham_scores(ct_iso_var_sum_df, ct_anno_df, ct_beta_tubs, "tbb-2", "c_tropicalis"))
)

all_beta_tub_blossum_grantham

data.table::fwrite(all_beta_tub_blossum_grantham, subscores_out_name, sep = "\t", quote = FALSE, col.names = TRUE)

# Create a combined dataframe for c_elegans
ce_df <- all_beta_tub_blossum_grantham %>%
  dplyr::filter(Species == "c_elegans") %>%
  dplyr::left_join(ce_hta_data, by = "strain") %>%
  dplyr::left_join(ce_iso_var_sum_df %>% dplyr::select(strain, abz_hta_norm_pheno, abz_hta_norm_res), by = "strain")

# Create a combined dataframe for c_briggsae
cb_df <- all_beta_tub_blossum_grantham %>%
  dplyr::filter(Species == "c_briggsae") %>%
  dplyr::left_join(cb_hta_data, by = "strain")

ct_df <- all_beta_tub_blossum_grantham %>%
  dplyr::filter(Species == "c_tropicalis") %>%
  dplyr::left_join(ct_hta_data, by = "strain")

# Create scatter plots for each species and metric
ce_plot_list <- list(
  create_scatter_plot(
    ce_df,
    "BLOSUM",
    "mean_median_wormlength_um_delta_reg",
    "BLOSUM",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    ce_df,
    "Grantham",
    "mean_median_wormlength_um_delta_reg",
    "Grantham",
    "Mean Median Wormlength Delta"
  ), #+ ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    ce_df,
    "Percent_Protein",
    "mean_median_wormlength_um_delta_reg",
    "Percent Protein",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank())
)

ce_abz_plot_list <- list(
  create_scatter_plot(
    ce_df,
    "BLOSUM",
    "abz_hta_norm_pheno",
    "BLOSUM",
    "Normalized Animal Length (µm)"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    ce_df,
    "Grantham",
    "abz_hta_norm_pheno",
    "Grantham",
    "Normalized Animal Length (µm)"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    ce_df,
    "Percent_Protein",
    "abz_hta_norm_pheno",
    "Percent Protein",
    "Normalized Animal Length (µm)"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank())
)

cb_plot_list <- list(
  create_scatter_plot(
    cb_df,
    "BLOSUM",
    "mean_median_wormlength_um_delta_reg",
    "BLOSUM",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    cb_df,
    "Grantham",
    "mean_median_wormlength_um_delta_reg",
    "Grantham",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    cb_df,
    "Percent_Protein",
    "mean_median_wormlength_um_delta_reg",
    "Percent Protein",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank())
)

ct_plot_list <- list(
  create_scatter_plot(
    ct_df,
    "BLOSUM",
    "mean_median_wormlength_um_reg_delta",
    "BLOSUM",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    ct_df,
    "Grantham",
    "mean_median_wormlength_um_reg_delta",
    "Grantham",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank()),
  create_scatter_plot(
    ct_df,
    "Percent_Protein",
    "mean_median_wormlength_um_reg_delta",
    "Percent Protein",
    "Mean Median Wormlength Delta"
  ) + ggplot2::theme(axis.title.y = ggplot2::element_blank())
)

# Combine the plots into a single figure with shared y-axis title and a single legend

## Define theme for pub. combined plot ##

theme_pub <- function() {
  theme_bw() +
    theme(
      # remove gridlines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(family = "Helvetica", size = 11),
      ## text formatting ##
      # bold x-axis title 
      axis.title.x = element_text(face = "bold"),
      # bold y-axis title
      axis.title.y = element_blank(),
      # bold legend title
      legend.title = element_text(face = "bold"),
      # make x and y axis text black
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}

# Function to create combined plot with legend

ce_combined_p <- create_faceted_scatter_plot(ce_df)

cb_combined_p <- create_faceted_scatter_plot(cb_df)

ct_combined_p <- create_faceted_scatter_plot(
  ct_df  %>% 
    dplyr::rename(mean_median_wormlength_um_delta_reg = mean_median_wormlength_um_reg_delta)
  )

# Save the combined plots using the save_plot function
save_plot(
  tplot = ce_combined_p,
  fn_list = ce_missense_fn,
  w_in = 7.5,
  h_in = 5
)
# save_plot(
#   tplot = ce_abz_combined_plot_with_legend,
#   fn_list = ce_abz_missense_fn,
#   w_in = 7.5,
#   h_in = 7
# )
save_plot(
  tplot = cb_combined_p,
  fn_list = cb_missense_fn,
  w_in = 7.5,
  h_in = 5
)
save_plot(
  tplot = ct_combined_p,
  fn_list = ct_missense_fn,
  w_in = 7.5,
  h_in = 5
)

# # Create an extra plot where ce_combined and ce_abz_combined are presented in two columns on the same plot
# combined_ce_abz_plot <- cowplot::plot_grid(
#   ce_combined_plot_with_legend,
#   ce_abz_combined_plot_with_legend,
#   ncol = 2,
#   labels = c("A", "B"),
#   label_size = 12,
#   label_fontfamily = "Helvetica",
#   label_fontface = "bold"
# )

# save_plot(
#   tplot = combined_ce_abz_plot,
#   fn_list = c("figures/figure_S14/figure_S14.png", "figures/figure_S14/figure_S14.eps"),
#   w_in = 15,
#   h_in = 7
# )

