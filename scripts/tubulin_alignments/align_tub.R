# Install ggmsa
# devtools::install_github("YuLab-SMU/ggmsa")
library(ggmsa)
library(glue)
library(Biostrings) # Added for reading FASTA files
library(ggtext) # Added for Markdown in axis text

# Load standardized output functions and figure mappings
source("bin/outs.R")
source("bin/materials_key.R")

#### Inputs #### ----------
alignments_dir <- "data/proc/tubulin_alignments/20251002_tubulin_protein_seqs_20251002"

#### Outputs #### ----------
# Output file paths are now managed by materials_key.R: sup_figure_fns$combined_plot

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

#### Function to load and plot tubulin alignment ####
load_and_plot_alignment <- function(gene_name, alignments_dir, start = 175, end = 225) {
  # Construct path to the MSA file
  aln_file <- glue::glue("{alignments_dir}/{gene_name}.fa")

  # Read and transform sequence names
  seqs <- readAAStringSet(aln_file)
  names(seqs) <- sapply(names(seqs), transform_seq_name)

  # Create the alignment plot
  align_plot <- ggmsa(
    seqs,
    start = start,
    end = end,
    color = "Chemistry_AA",
    seq_name = TRUE
  )

  return(align_plot)
}


#### Generate alignment plots for each gene ####
tbb1_align_zoom <- load_and_plot_alignment("tbb_1", alignments_dir)
tbb2_align_zoom <- load_and_plot_alignment("tbb_2", alignments_dir)
ben1_align_zoom <- load_and_plot_alignment("ben_1", alignments_dir)
mec7_align_zoom <- load_and_plot_alignment("mec_7", alignments_dir)
tbb4_align_zoom <- load_and_plot_alignment("tbb_4", alignments_dir)


#### combine all three plots ####

# create a list of the three plots
plot_list <- list(tbb1_align_zoom, tbb2_align_zoom, mec7_align_zoom, tbb4_align_zoom, ben1_align_zoom)

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

# save the plot using standardized save_plot function
save_plot(
  tplot = combined_plot,
  fn_list = sup_figure_fns$combined_plot,
  w_in = 7.5,
  h_in = 5
)
