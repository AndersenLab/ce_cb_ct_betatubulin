library(data.table)
library(tidyverse)

## Inputs ##
# URL to the gene_data.txt file from the C.elegans_C.briggsae_Embryo_Single_Cell repo
# Using GitHub's media URL to download the actual LFS file content
gene_summary_url <- "https://media.githubusercontent.com/media/livinlrg/C.elegans_C.briggsae_Embryo_Single_Cell/main/Summary_Tables/gene_data.txt"


gene_summary <- data.table::fread(gene_summary_url)

# print(head(gene_summary))

# filter for beta tubulin genes
core_tub_gene_data <- gene_summary %>%
  dplyr::filter(
    gene %in% c(
      "ben-1", "tbb-1", "tbb-2", "mec-7", "tbb-4"
    )
  ) %>%
  dplyr::select(
    gene, jsd_median_term, jsd_lower_term, jsd_upper_term
  )

print(core_tub_gene_data)
