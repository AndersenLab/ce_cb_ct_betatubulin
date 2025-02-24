#### libs ####
library(glue)
library(data.table)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)


#### Load data ####


#### Plot panel A ####

# # Filter to only strains with high-impact variants in
# # either ben-1, tbb-1, or tbb-2.
# # add a tag for the species
# # and remove any strains that dont have lat long
# ce_ben1_var_df <- ce_df %>%
#   dplyr::filter(
#     `ben-1_var` != "NA | NA"
#   ) %>%
#   dplyr::mutate(species = "C. elegans") %>%
#   dplyr::filter(!is.na(latitude) | !is.na(longitude))

# n_ce_ben1_var <- nrow(ce_ben1_var_df)

# print(
#   glue::glue(
#     "Number of C. elegans strains with high-impact ben-1 variants and coordinates: {n_ce_ben1_var}"
#   )
# )

# cb_ben1_var_df <- cb_df %>%
#   dplyr::filter(
#     `ben-1_var` != "NA | NA"
#   ) %>%
#   dplyr::mutate(species = "C. briggsae") %>%
#   dplyr::filter(!is.na(latitude) | !is.na(longitude))

# n_cb_ben1_var <- nrow(cb_ben1_var_df)

# print(
#   glue::glue(
#     "Number of C. briggsae strains with high-impact ben-1 variants and coordinates: {n_cb_ben1_var}"
#   )
# )

# ct_ben1_var_df <- ct_df %>%
#   dplyr::filter(
#     `ben-1_var` != "NA | NA"
#   ) %>%
#   dplyr::mutate(species = "C. tropicalis") %>%
#   dplyr::filter(!is.na(latitude) | !is.na(longitude))

# n_ct_ben1_var <- nrow(ct_ben1_var_df)

# print(
#   glue::glue(
#     "Number of C. tropicalis strains with high-impact ben-1 variants and coordinates: {n_ct_ben1_var}"
#   )
# )
#' Filter isotype variant summary data frame for high-impact variants and coordinates
#'
#' @param df A data frame containing isotype variant summary data.
#' @param variant_col The name of the column containing variant information.
#' @param species_id A string denoting the species ("C. elegans", "C. briggsae", or "C. tropicalis").
#' @return A filtered data frame with only strains that have high-impact variants and coordinates, with an added species column.
#' @import dplyr
filter_high_impact_variants <- function(df, variant_col, species_id) {
  filtered_df <- df %>%
    dplyr::filter(
      !!sym(variant_col) != "No variant"
    ) %>%
    dplyr::mutate(
      species = species_id,
      !!variant_col := as.character(!!sym(variant_col))
    ) %>%
    dplyr::filter(!is.na(latitude) & !is.na(longitude)) %>%
    dplyr::select(strain, isotype, latitude, longitude, !!sym(variant_col), species)

  n_filtered <- nrow(filtered_df)
  print(glue::glue("Number of {species_id} strains with high-impact {variant_col} variants and coordinates: {n_filtered}"))

  return(filtered_df)
}

# # Apply the function to each species data frame
# ce_ben1_var_df <- filter_high_impact_variants(ce_df, "ben-1_var", "C. elegans")
# cb_ben1_var_df <- filter_high_impact_variants(cb_df, "ben-1_var", "C. briggsae")
# ct_ben1_var_df <- filter_high_impact_variants(ct_df, "ben-1_var", "C. tropicalis")

#' Combine filtered data frames
#'
#' @param ... Data frames output by the `filter_high_impact_variants` function.
#' @return A combined data frame from all input data frames.
#' @import dplyr
combine_filtered_data_frames <- function(...) {
  combined_df <- dplyr::bind_rows(...) %>%
    # Set the levels to the species for plotting
    dplyr::mutate(
      species = factor(species, levels = c("C. elegans", "C. briggsae", "C. tropicalis"))
    )
  return(combined_df)
}

#' Convert data frame to sf object
#'
#' @param df A data frame containing combined data from all species.
#' @return A list containing an sf object and the world basemap.
#' @import sf
#' @import rnaturalearth
#' @import dplyr
convert_to_sf <- function(df) {
  # Load basemap
  world <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf"
  ) %>%
    dplyr::filter(geounit != "Antarctica")

  # Convert the combined data frame to sf object
  all_var_sf <- sf::st_as_sf(
    df,
    coords = c("longitude", "latitude"),
    crs = sf::st_crs(world),
    remove = TRUE
  )

  return(list(all_var_sf = all_var_sf, world = world))
}

# # Combine the filtered data frames
# combined_df <- combine_filtered_data_frames(ce_ben1_var_df, cb_ben1_var_df, ct_ben1_var_df)

# # Convert the combined data frame to sf object
# sf_objects <- convert_to_sf(combined_df)
# all_var_sf <- sf_objects$all_var_sf
# world <- sf_objects$world

# # Apply the function to the filtered data frames
# sf_objects <- combine_and_convert_to_sf(ce_ben1_var_df, cb_ben1_var_df, ct_ben1_var_df)
# all_var_sf <- sf_objects$all_var_sf
# world <- sf_objects$world

#' Plot high-impact variants on a world map
#'
#' @param world An sf object containing the world map.
#' @param all_var_sf An sf object containing the combined data of all species with high-impact variants.
#' @return A ggplot object representing the map with high-impact variants.
#' @import ggplot2
plot_high_impact_variants_map <- function(world, all_var_sf) {
  all_var_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = "white") +
    ggplot2::geom_sf(
      data = all_var_sf,
      ggplot2::aes(fill = species),
      shape = 21,
      size = 1.5,
      alpha = 0.7
      ) +
    ggplot2::scale_fill_manual(
      values = c(
        "C. elegans" = "#FFA500",
        "C. briggsae" = "#53886c",
        "C. tropicalis" = "#0619BC"
      ),
      labels = c(
        "C. elegans" = expression(italic("C. elegans")),
        "C. briggsae" = expression(italic("C. briggsae")),
        "C. tropicalis" = expression(italic("C. tropicalis"))
      ),
      name = "Species"
    ) +
    # ggplot2::scale_color_manual(
    #   values = c(
    #     "C. elegans" = "#FFA500",
    #     "C. briggsae" = "#53886c",
    #     "C. tropicalis" = "#0619BC"
    #   ),
    #   labels = c(
    #     "C. elegans" = expression(italic("C. elegans")),
    #     "C. briggsae" = expression(italic("C. briggsae")),
    #     "C. tropicalis" = expression(italic("C. tropicalis"))
    #   ),
    #   name = "Species"
    # ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = c(.1, .2),
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      text = ggplot2::element_text(family = "Arial", size = 11)
    )

  return(all_var_map)
}

# Generate the plot using the function
# all_var_map <- plot_high_impact_variants_map(world, all_var_sf)

#### Save plot ####

#' Save a plot of high-impact variants on a world map
#'

save_plot <- function(fn_list, plot) {
  
  # get the folder name from the first file name
  folder <- dirname(fn_list$eps)
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  } else {

  }

  fn_png <- fn_list$png
  fn_eps <- fn_list$eps

  # Save as eps
  ggplot2::ggsave(
    filename = glue::glue("{fn_eps}"),
    plot = plot,
    width = 7,
    height = 4,
    units = "in",
    dpi = 300
  )

  # Save as png
  ggplot2::ggsave(
    filename = glue::glue("{fn_png}"),
    plot = plot,
    width = 7,
    height = 4,
    units = "in",
    dpi = 300
  )
}

# # Example usage
# save_plot("ben_1", figure_out_dir, all_var_map)
