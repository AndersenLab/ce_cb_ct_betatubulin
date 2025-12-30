library(tidyverse)
library(readxl)


d1_fn <- "data/expression_conservation/table_s1.xlsx"
beta_tubulin_genes <- c("tbb-1", "tbb-2", "ben-1", "tbb-4", "mec-7")


# read in all sheets
sheets <- excel_sheets(d1_fn)

data_list <- lapply(sheets, function(sheet) {
  read_excel(d1_fn, sheet = sheet)
})

names(data_list) <- sheets

jaccard_df <- data_list[["I - Jaccard distances"]]

# filter to beta-tubulin genes
beta_jaccard_df <- jaccard_df %>%
  filter(Gene_triplet %in% beta_tubulin_genes)

# save output
output_fn <- "tables/table_s9_neuronal_beta_tubulin_jaccard_distances.tsv"
write_tsv(beta_jaccard_df, output_fn)
message("Wrote output to ", output_fn)
