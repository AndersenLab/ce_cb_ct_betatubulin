#!/usr/bin/env Rscript

# Combine missense figures into a 2x2 panel figure for manuscript
# Panel layout:
#   A (elegans)    B (briggsae)
#   C (tropicalis) D (contrortus)

library(png)
library(grid)
library(gridExtra)

# Load standardized output functions and figure mappings
source("bin/outs.R")
source("bin/materials_key.R")

# Define input files
panel_files <- c(
  a = "data/structure/raw_images/elegans_missense.png",
  b = "data/structure/raw_images/briggsae_missense.png",
  c = "data/structure/raw_images/tropicalis_missense.png",
  d = "data/structure/raw_images/contrortus_missense.png"
)

# Check that all files exist
missing <- panel_files[!file.exists(panel_files)]
if (length(missing) > 0) {
  stop("Missing files: ", paste(names(missing), missing, sep = ": ", collapse = ", "))
}

# Read PNG images
images <- lapply(panel_files, readPNG)

# Convert to grobs with panel labels
create_labeled_panel <- function(img, label) {
  img_grob <- rasterGrob(img, interpolate = TRUE)

  # Create label grob (top-left corner)
  label_grob <- textGrob(
    label,
    x = unit(0.02, "npc"),
    y = unit(0.98, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )

  # Combine image and label
  gTree(children = gList(img_grob, label_grob))
}

# Create labeled panels
panels <- mapply(create_labeled_panel, images, names(panel_files), SIMPLIFY = FALSE)

# Arrange in 2x2 grid
missense_structure_combined <- arrangeGrob(
  grobs = panels,
  nrow = 2,
  ncol = 2,
  padding = unit(0.5, "line")
)

# Figure dimensions in inches (width: 2.63-7.5 in, height: <8.75 in)
fig_width <- 7 # inches
fig_height <- 7 # inches

# Save combined figure using standardized save_grid_plot function
save_grid_plot(
  grid_plot = missense_structure_combined,
  fn_list = sup_figure_fns$missense_structure_combined,
  w_in = fig_width,
  h_in = fig_height
)

cat("Combined figure saved to:", sup_figure_fns$missense_combined$png, "\n")
