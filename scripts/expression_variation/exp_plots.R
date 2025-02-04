# Scipts to plot expression variation

#' Create a Scatter Plot of Beta-Tubulin Expression vs. Normalized Response
#'
#' Generates a scatter plot of beta-tubulin expression levels against normalized response.
#' Points are colored based on a specified categorical variable.
#' Includes a linear regression line and a horizontal threshold line.
#'
#' @param exp_data A data frame containing the expression data.
#' @param x_column_id A string specifying the name of the column for the x-axis variable (expression levels).
#' @param y_column_id A string specifying the name of the column for the y-axis variable (normalized response).
#' @param fill_column_id A string specifying the name of the column for the categorical variable used to fill colors.
#' @param fill_scale A named vector specifying colors for each category in the `fill_column_id`.
#' @param x_label A string or expression for the x-axis label.
#' @param y_label A string for the y-axis label.
#' @param fill_label A string or expression for the legend title.
#' @param res_threshold A numeric value specifying the y-axis threshold for drawing a horizontal line.
#'
#' @return A ggplot object representing the scatter plot.
#' @import ggplot2
#' @examples
#' create_expression_scatter_plot(
#'   exp_data = exp_summary,
#'   x_column_id = "ben-1_exp",
#'   y_column_id = "strain_norm_abz_response",
#'   fill_column_id = "ben1_var_cat_meta",
#'   fill_scale = meta_cat_cols,
#'   x_label = expression(italic("ben-1") * " expression (TPM)"),
#'   y_label = "Normalized ABZ Response",
#'   fill_label = expression(italic("ben-1") * " consequence"),
#'   res_threshold = norm_res_threshold
#' )
create_expression_scatter_plot <- function(
    exp_data,
    x_column_id,
    y_column_id,
    fill_column_id,
    fill_scale,
    x_label,
    y_label,
    fill_label,
    res_threshold) {
  # Perform linear regression
  lm_formula <- as.formula(glue::glue("`{y_column_id}` ~ `{x_column_id}`"))
  lm_model <- lm(lm_formula, data = exp_data)
  summary_lm <- summary(lm_model)
  p_value <- summary_lm$coefficients[2, 4]
  r_squared <- summary_lm$r.squared

  plot <- ggplot2::ggplot(
    exp_data,
    ggplot2::aes(
      x = !!sym(x_column_id),
      y = !!sym(y_column_id)
    )
  ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      color = "grey"
    ) +
    ggplot2::geom_hline(
      yintercept = res_threshold,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = !!sym(fill_column_id)
      ),
      shape = 21,
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = fill_scale
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      fill = fill_label,
      caption = glue::glue("p-value: {format(p_value, digits = 3)}, R²: {format(r_squared, digits = 3)}")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )

  out = list(plot = plot, p_value = p_value, r_squared = r_squared)

  return(out)
}

#' Create a Boxplot of Expression Data
#'
#' Generates a boxplot of expression data with jittered points.
#' Points are colored based on a specified categorical variable.
#' Includes statistical comparisons between groups.
#'
#' @param data A data frame containing the expression data.
#' @param x_col A string specifying the name of the column for the x-axis variable.
#' @param y_col A string specifying the name of the column for the y-axis variable.
#' @param x_label A string or expression for the x-axis label.
#' @param y_label A string or expression for the y-axis label.
#' @param comparisons_list A list of character vectors specifying the groups to compare.
#' @param fill_column A string specifying the name of the column for the fill variable.
#' @param fill_scale A named vector specifying colors for each category in the `fill_column`.
#'
#' @return A ggplot object representing the boxplot.
create_expression_boxplot <- function(data, x_col, y_col, x_label, y_label, comparisons_list, fill_column, fill_scale) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      fill = .data[[fill_column]]
    )
  ) +
    ggplot2::geom_boxplot(
      outliers = FALSE
    ) +
    ggplot2::geom_jitter(
      width = 0.2,
      height = 0,
      size = 2
    ) +
    ggplot2::scale_fill_manual(values = fill_scale) +
    ggplot2::labs(
      x = x_label,
      y = y_label
    ) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      paired = FALSE,
      label = "p.signif",
      comparisons = comparisons_list
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#### Plot ben-1 expression x tbb-1/tbb-2 expression ####

plot_exp_x_exp <- function(
    data,
    x_col,
    y_col,
    fill_col,
    shape_col,
    x_label,
    y_label,
    fill_legend_label,
    shape_legend_label) {
  # Check if the specified columns exist in the data frame
  if (!(x_col %in% colnames(data))) {
    stop(glue::glue("Column '{x_col}' not found in the data frame"))
  }
  if (!(y_col %in% colnames(data))) {
    stop(glue::glue("Column '{y_col}' not found in the data frame"))
  }
  if (!(fill_col %in% colnames(data))) {
    stop(glue::glue("Column '{fill_col}' not found in the data frame"))
  }
  if (!(shape_col %in% colnames(data))) {
    stop(glue::glue("Column '{shape_col}' not found in the data frame"))
  }

  # Perform linear regression
  lm_formula <- as.formula(glue::glue("`{y_col}` ~ `{x_col}`"))
  lm_model <- lm(lm_formula, data = data)
  summary_lm <- summary(lm_model)
  p_value <- summary_lm$coefficients[2, 4]
  r_squared <- summary_lm$r.squared

  # Create the scatter plot with linear regression line
  scatter_plot <- data %>%
    ggplot(aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(fill = !!sym(fill_col), shape = !!sym(shape_col)), size = 3, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "grey") +
    labs(
      x = x_label,
      y = y_label,
      fill = fill_legend_label,
      shape = shape_legend_label,
      caption = glue::glue("p-value: {format(p_value, digits = 3)}, R²: {format(r_squared, digits = 3)}")
    ) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21)) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "white")) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) + # Force the fill legend
    theme_minimal()

  return(scatter_plot)
}

# add the ben-1 metadata to the expression summary w/ normalized phenos
#' Add Ben-1 Metadata
add_ben1_metadata <- function(df) {
  ben1_var_cats_meta <- c(
    "Deletion" = "SV",
    "Duplication" = "SV",
    "Inversion" = "SV",
    "Transposon insertion" = "SV",
    "Frameshift" = "Frame altering",
    "Inframe_deletion" = "Frame altering",
    "Missense" = "Missense",
    "Start/stop altering" = "Start/stop altering",
    "No variant" = "No variant",
    "Low ben-1 expression" = "Low ben-1 expression"
  )

  # Check if all unique values in `ben-1_clean_call` have an assigned value in `ben1_var_cats_meta`
  unique_calls <- unique(df$`ben-1_clean_call`)
  missing_calls <- setdiff(unique_calls, names(ben1_var_cats_meta))

  if (length(missing_calls) > 0) {
    stop(glue::glue("The following `ben-1_clean_call` values are not present in `ben1_var_cats_meta`: {paste(missing_calls, collapse = ', ')}"))
  }

  df <- df %>%
    mutate(
      ben1_var_cat_meta = factor(`ben-1_clean_call`, levels = names(ben1_var_cats_meta), labels = ben1_var_cats_meta)
    )

  return(df)
}
