library(data.table)
library(tidyverse)

## Inputs ##
# URL to the gene_data.txt file from the C.elegans_C.briggsae_Embryo_Single_Cell repo
# Using GitHub's media URL to download the actual LFS file content
gene_summary_url <- "https://media.githubusercontent.com/media/livinlrg/C.elegans_C.briggsae_Embryo_Single_Cell/main/Summary_Tables/gene_data.txt"


gene_summary <- data.table::fread(gene_summary_url)

print(colnames(gene_summary))

# filter for beta tubulin genes
core_tub_gene_data <- gene_summary %>%
  dplyr::filter(
    gene %in% c(
      "ben-1", "tbb-1", "tbb-2", "mec-7", "tbb-4"
    )
  )

## JSD values ##

# Terminal cell types
jsd_term <- core_tub_gene_data %>%
  dplyr::select(
    gene,
    jsd_median = jsd_median_term,
    jsd_lower = jsd_lower_term,
    jsd_upper = jsd_upper_term
  ) %>%
  dplyr::mutate(cell_population = "Terminal")

# Progenitor cell types
jsd_pro <- core_tub_gene_data %>%
  dplyr::select(
    gene,
    jsd_median = jsd_median_pro,
    jsd_lower = jsd_lower_pro,
    jsd_upper = jsd_upper_pro
  ) %>%
  dplyr::mutate(cell_population = "Progenitor")

# joint - all cell types
jsd_joint <- core_tub_gene_data %>%
  dplyr::select(
    gene,
    jsd_median = jsd_median_joint,
    jsd_lower = jsd_lower_joint,
    jsd_upper = jsd_upper_joint
  ) %>%
  dplyr::mutate(cell_population = "Joint")


jsd_values <- rbind(jsd_term, jsd_pro, jsd_joint)

# Set the gene order (tbb-1, tbb-2, ben-1, mec-7, tbb-4)
jsd_values$gene <- factor(jsd_values$gene, levels = c("tbb-1", "tbb-2", "ben-1", "mec-7", "tbb-4"))

## Tau values ##

# terminal cell types
tau_term <- core_tub_gene_data %>%
  dplyr::select(
    gene,
    cel_tau_median = cel_tau_median_term,
    cel_tau_lower = cel_tau_lower_term,
    cel_tau_upper = cel_tau_upper_term,
    cbr_tau_median = cbr_tau_median_term,
    cbr_tau_lower = cbr_tau_lower_term,
    cbr_tau_upper = cbr_tau_upper_term
  ) %>%
  dplyr::mutate(cell_population = "Terminal")

# progenitor cell types
tau_pro <- core_tub_gene_data %>%
  dplyr::select(
    gene,
    cel_tau_median = cel_tau_median_pro,
    cel_tau_lower = cel_tau_lower_pro,
    cel_tau_upper = cel_tau_upper_pro,
    cbr_tau_median = cbr_tau_median_pro,
    cbr_tau_lower = cbr_tau_lower_pro,
    cbr_tau_upper = cbr_tau_upper_pro
  ) %>%
  dplyr::mutate(cell_population = "Progenitor")

# joint - all cell types
tau_joint <- core_tub_gene_data %>%
  dplyr::select(
    gene,
    cel_tau_median = cel_tau_median_joint,
    cel_tau_lower = cel_tau_lower_joint,
    cel_tau_upper = cel_tau_upper_joint,
    cbr_tau_median = cbr_tau_median_joint,
    cbr_tau_lower = cbr_tau_lower_joint,
    cbr_tau_upper = cbr_tau_upper_joint
  ) %>%
  dplyr::mutate(cell_population = "Joint")


tau_values <- rbind(tau_term, tau_pro, tau_joint)

# Set the gene order (tbb-1, tbb-2, ben-1, mec-7, tbb-4)
tau_values$gene <- factor(tau_values$gene, levels = c("tbb-1", "tbb-2", "ben-1", "mec-7", "tbb-4"))


tau_values <- tau_values %>%
  dplyr::mutate(tau_diff = abs(cel_tau_median - cbr_tau_median))

# Create supplementary tables of JSD and Tau values of all cells

## JSD


# JSD supplementary table
jsd_supp_table <- jsd_values %>%
  dplyr::select(gene, cell_population, jsd_median, jsd_lower, jsd_upper)


# save as supplementary table 14
data.table::fwrite(
  jsd_supp_table,
  file = "tables/table_s14_jsd_values_beta_tubulin_genes.csv"
)

## Tau

# Tau supplementary table
tau_supp_table <- tau_values %>%
  dplyr::select(
    gene, cell_population, cel_tau_median, cel_tau_lower, cel_tau_upper,
    cbr_tau_median, cbr_tau_lower, cbr_tau_upper, tau_diff
  )

# save as supplementary table 14
data.table::fwrite(
  tau_supp_table,
  file = "tables/table_s15_tau_values_beta_tubulin_genes.csv"
)


# Mean and SD of JSD values
jsd_stats <- jsd_values %>%
  dplyr::group_by(cell_population) %>%
  dplyr::summarise(
    jsd_mean = mean(jsd_median),
    jsd_sd = sd(jsd_median)
  )

print("JSD Statistics:")
print(jsd_stats)

# Mean and SD of Tau values and difference
tau_stats <- tau_values %>%
  dplyr::group_by(cell_population) %>%
  dplyr::summarise(
    cel_tau_mean = mean(cel_tau_median),
    cel_tau_sd = sd(cel_tau_median),
    cbr_tau_mean = mean(cbr_tau_median),
    cbr_tau_sd = sd(cbr_tau_median),
    tau_diff_mean = mean(tau_diff),
    tau_diff_sd = sd(tau_diff)
  )

print("Tau Statistics:")
print(tau_stats)
