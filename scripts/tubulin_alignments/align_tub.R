# Install ggmsa
# devtools::install_github("YuLab-SMU/ggmsa")
library(ggmsa)
library(glue)

#### Inputs #### ----------
alignments_dir <- "data/tubulin_alignments/20240828_tubulin_protein_seqs_20240828"


#### Outputs #### ----------
out_dir <- glue::glue("figures/tubulin_alignments/tub_alignments")

# check if the directory exists
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

#### tbb-1 ####
# set path to the MSA file
tbb1_aln_file <- glue::glue("{alignments_dir}/tbb_1.fa")

# create the alignment plot

tbb1_align_zoom <- ggmsa(
    tbb1_aln_file,
    start = 175,
    end = 225,
    # position_highlight = 200, - will highlight the position 200 - no other colors
    color = "Clustal",
    seq_name = TRUE
) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 5, b = 2, l = 5, unit = "pt")) # Adjust margins: top, right, bottom, left

# save the plot
ggplot2::ggsave(
    glue::glue("{out_dir}/tbb1_aln_zoom.jpg"),
    tbb1_align_zoom,
    width = 7.5,
    height = 7.5,
    dpi = 300
)

#### tbb-2 ####
# set path to the MSA file
tbb2_aln_file <- glue::glue("{alignments_dir}/tbb_2.fa")

tbb2_align_zoom <- ggmsa(
    tbb2_aln_file,
    start = 175,
    end = 225,
    # position_highlight = 200, - will highlight the position 200 - no other colors
    color = "Clustal",
    seq_name = TRUE
) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 5, b = 2, l = 5, unit = "pt")) # Adjust margins

# save the plot
ggplot2::ggsave(
    glue::glue("{out_dir}/tbb2_aln_zoom.jpg"),
    tbb2_align_zoom,
    width = 7.5,
    height = 7.5,
    dpi = 300
)

#### ben-1 ####
# set path to the MSA file
ben1_aln_file <- glue::glue("{alignments_dir}/ben_1.fa")

ben1_align_zoom <- ggmsa(
    ben1_aln_file,
    start = 175,
    end = 225,
    # position_highlight = 200, - will highlight the position 200 - no other colors
    color = "Clustal",
    seq_name = TRUE
) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 5, b = 2, l = 5, unit = "pt")) # Adjust margins

# save the plot
ggplot2::ggsave(
    glue::glue("{out_dir}/ben1_aln_zoom.jpg"),
    ben1_align_zoom,
    width = 7.5,
    height = 7.5,
    dpi = 300
)

#### combine all three plots ####

# create a list of the three plots
plot_list <- list(tbb1_align_zoom, tbb2_align_zoom, ben1_align_zoom)

# combine the plots
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, align = "v")

# save the plot
# Consider reducing the height if the plots are now more compact
ggplot2::ggsave(
    glue::glue("{out_dir}/tub_alignments_combined.jpg"),
    combined_plot,
    width = 7.5,
    height = 3, # Adjust height as needed
    dpi = 300
)
