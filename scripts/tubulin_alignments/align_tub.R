# Install ggmsa
# devtools::install_github("YuLab-SMU/ggmsa")
library(ggmsa)
library(glue)
library(Biostrings) # Added for reading FASTA files
library(ggtext) # Added for Markdown in axis text

#### Inputs #### ----------
alignments_dir <- "data/tubulin_alignments/20240828_tubulin_protein_seqs_20240828"


#### Outputs #### ----------
out_dir <- glue::glue("figures/figure_S1")

# check if the directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#### Helper function to transform sequence names ####
transform_seq_name <- function(original_name) {
  parts <- strsplit(original_name, "_")[[1]]
  if (length(parts) != 2) {
    return(original_name) # Return original if format is unexpected
  }

  prefix_code <- parts[1]
  suffix <- parts[2]

  new_prefix <- ""
  if (prefix_code == "CE") {
    new_prefix <- "Cel"
  } else if (prefix_code == "CB") {
    new_prefix <- "Cbr"
  } else if (prefix_code == "CT") {
    new_prefix <- "Ctr"
  } else {
    stop(
      "Unexpected prefix code"
    ) # Keep original prefix if not matched
  }

  uppercase_suffix <- toupper(
    suffix
  )

  return(
    paste0(
      "*",
      new_prefix,
      "*-",
      uppercase_suffix
    )
  )
}


#### tbb-1 ####
# set path to the MSA file
tbb1_aln_file <- glue::glue("{alignments_dir}/tbb_1.fa")

# Read and transform sequence names
tbb1_seqs <- readAAStringSet(tbb1_aln_file)
names(tbb1_seqs) <- sapply(names(tbb1_seqs), transform_seq_name)

# create the alignment plot
tbb1_align_zoom <- ggmsa(
  tbb1_seqs, # Use modified sequences
  start = 175,
  end = 225,
  # position_highlight = 200, - will highlight the position 200 - no other colors
  color = "Chemistry_AA",
  seq_name = TRUE
)


#### tbb-2 ####
# set path to the MSA file
tbb2_aln_file <- glue::glue("{alignments_dir}/tbb_2.fa")

# Read and transform sequence names
tbb2_seqs <- readAAStringSet(tbb2_aln_file)
names(tbb2_seqs) <- sapply(names(tbb2_seqs), transform_seq_name)

tbb2_align_zoom <- ggmsa(
  tbb2_seqs, # Use modified sequences
  start = 175,
  end = 225,
  # position_highlight = 200, - will highlight the position 200 - no other colors
  color = "Chemistry_AA",
  seq_name = TRUE
)

#### ben-1 ####
# set path to the MSA file
ben1_aln_file <- glue::glue("{alignments_dir}/ben_1.fa")

# Read and transform sequence names
ben1_seqs <- readAAStringSet(ben1_aln_file)
names(ben1_seqs) <- sapply(names(ben1_seqs), transform_seq_name)

ben1_align_zoom <- ggmsa(
  ben1_seqs, # Use modified sequences
  start = 175,
  end = 225,
  # position_highlight = 200, - will highlight the position 200 - no other colors
  color = "Chemistry_AA",
  seq_name = TRUE
)

#### combine all three plots ####

# create a list of the three plots
plot_list <- list(tbb1_align_zoom, tbb2_align_zoom, ben1_align_zoom)

# add the same theme to all plots
plot_list <- lapply(plot_list, function(p) {
  p +
    ggplot2::theme_bw() +

    ggplot2::theme(
      plot.margin = ggplot2::margin(
        t = 2,
        r = 5,
        b = 2,
        l = 5,
        unit = "pt"
      ),
      axis.text.y = ggtext::element_markdown(
        family = "Arial", # Set font family for y-axis labels
        size = 11, # Set font size for y-axis labels
        color = "black" # Set color for y-axis labels
      ),
      axis.text.x = ggtext::element_markdown(
        family = "Arial", # Set font family for x-axis labels
        size = 11, # Set font size for x-axis labels
        color = "black" # Set color for x-axis labels
      ),
    )
})

# combine the plots
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, align = "v")

# save the plot
# Consider reducing the height if the plots are now more compact
ggplot2::ggsave(
  glue::glue("{out_dir}/figure_S1.jpg"),
  combined_plot,
  width = 7.5,
  height = 3, # Adjust height as needed. You might need more height if names are longer.
  dpi = 300
)
