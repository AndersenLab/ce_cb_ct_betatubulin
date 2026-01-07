#!/usr/bin/env Rscript

# Combine missense figures into a 2x2 panel figure for manuscript
# Panel layout:
#   A (elegans)    B (briggsae)
#   C (tropicalis) D (contrortus)

library(png)
library(grid)
library(gridExtra)

# Define input files
panel_files <- c(
  a = "figures/elegans_missense.png",
  b = "figures/briggsae_missense.png",
  c = "figures/tropicalis_missense.png",
  d = "figures/contrortus_missense.png"
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
combined <- arrangeGrob(
  grobs = panels,
  nrow = 2,
  ncol = 2,
  padding = unit(0.5, "line")
)

# Output file
output_file <- "figures/missense_combined.png"

# Figure dimensions in inches (adjust as needed for manuscript)
fig_width <- 7    # inches
fig_height <- 7   # inches
fig_res <- 600    # DPI

# Save combined figure
png(output_file, width = fig_width, height = fig_height, res = fig_res, units = "in")
grid.draw(combined)
dev.off()

cat("Combined figure saved to:", output_file, "\n")
