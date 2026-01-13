# Standard output functions for saving figures and tables
# All figures are saved as .tiff, .png, and .jpg at 600 dpi
# Width should be between 2.63 and 7.5 inches
# Height should be less than 8.75 inches

#' Save a ggplot2 plot in multiple formats
#'
#' Saves a plot as .tiff, .png, and .jpg files at 600 dpi with white background.
#'
#' @param tplot A ggplot2 plot object
#' @param fn_list List of filenames with png, jpg, and tiff paths
#' @param w_in Width of the plot in inches (should be between 2.63 and 7.5)
#' @param h_in Height of the plot in inches (should be less than 8.75)
#'
#' @return No return value, called for side effects
#'
#' @examples
#' fn_list <- list(
#'   png = "figures/S1_FIG/S1_FIG.png",
#'   jpg = "figures/S1_FIG/S1_FIG.jpg",
#'   tiff = "figures/S1_FIG/S1_FIG.tiff"
#' )
#' save_plot(my_plot, fn_list, 7.5, 6)
save_plot <- function(tplot, fn_list, w_in, h_in) {
  # Create the output directories if they don't exist
  for (fn in fn_list) {
    folder <- dirname(fn)
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
  }

  fn_png <- fn_list$png
  fn_jpg <- fn_list$jpg
  fn_tiff <- fn_list$tiff

  # save png plot
  ggplot2::ggsave(
    filename = fn_png,
    plot = tplot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  # save jpg plot
  ggplot2::ggsave(
    filename = fn_jpg,
    plot = tplot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  # save tiff plot
  ggplot2::ggsave(
    filename = fn_tiff,
    plot = tplot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 600,
    bg = "white",
    compression = "lzw"
  )
}


#' Save a tree plot in multiple formats
#'
#' Saves a tree plot (from ggtree) as .tiff, .png, and .jpg files at 600 dpi.
#'
#' @param tree_plot A ggtree plot object
#' @param fn_list List of filenames with png, jpg, and tiff paths
#' @param w_in Width of the plot in inches (should be between 2.63 and 7.5)
#' @param h_in Height of the plot in inches (should be less than 8.75)
#'
#' @return No return value, called for side effects
save_tree <- function(tree_plot, fn_list, w_in, h_in) {
  # Use the same save_plot function for trees
  save_plot(tree_plot, fn_list, w_in, h_in)
}


#' Save a grid-based plot in multiple formats
#'
#' Saves a grid arrangement (from grid/gridExtra) as .tiff, .png, and .jpg files at 600 dpi.
#' This function is for plots created with grid.draw() rather than ggplot2.
#'
#' @param grid_plot A grid grob object (from arrangeGrob, gTree, etc.)
#' @param fn_list List of filenames with png, jpg, and tiff paths
#' @param w_in Width of the plot in inches (should be between 2.63 and 7.5)
#' @param h_in Height of the plot in inches (should be less than 8.75)
#'
#' @return No return value, called for side effects
save_grid_plot <- function(grid_plot, fn_list, w_in, h_in) {
  # Create the output directories if they don't exist
  for (fn in fn_list) {
    folder <- dirname(fn)
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
  }

  fn_png <- fn_list$png
  fn_jpg <- fn_list$jpg
  fn_tiff <- fn_list$tiff

  # Save as PNG
  png(fn_png, width = w_in, height = h_in, units = "in", res = 600, bg = "white")
  grid::grid.draw(grid_plot)
  dev.off()

  # Save as JPG
  jpeg(fn_jpg, width = w_in, height = h_in, units = "in", res = 600, bg = "white")
  grid::grid.draw(grid_plot)
  dev.off()

  # Save as TIFF
  tiff(fn_tiff, width = w_in, height = h_in, units = "in", res = 600, compression = "lzw", bg = "white")
  grid::grid.draw(grid_plot)
  dev.off()
}


#' Save data as a compressed TSV file
#'
#' @param data A data frame or data.table to save
#' @param output_file Path to the output file (should end in .tsv or .tsv.gz)
#'
#' @return No return value, called for side effects
save_compressed_tsv <- function(data, output_file) {
  # Create the output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save the data as a compressed .tsv file
  data.table::fwrite(data, output_file, sep = "\t", compress = "gzip")
}


#' Save data as a TSV file
#'
#' @param data A data frame or data.table to save
#' @param output_file Path to the output file
#'
#' @return No return value, called for side effects
save_tsv <- function(data, output_file) {
  # Create the output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save the data as a .tsv file
  data.table::fwrite(data, output_file, sep = "\t")
}


#' Save data as a CSV file
#'
#' @param data A data frame or data.table to save
#' @param output_file Path to the output file
#'
#' @return No return value, called for side effects
save_csv <- function(data, output_file) {
  # Create the output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save the data as a .csv file
  data.table::fwrite(data, output_file, sep = ",")
}
