# `beta_tub_expression_plots.R`
Plotting script for beta-tubulin expression and previous ABZ phenotypes
### inputs:

Hardcoded paths to:
- CE isotype variant summary file

- dataframe with previous ABZ phenotypes (sourced from `normalize_phenos.R`)
- outputs of `get_res_strains()` for normalized phenotype and both previous assays
### outputs:
#### plots (not for manuscript)
- Scatter plot of `ben-1` expression vs. normalized ABZ response, colored by `ben-1` variant category: `ben1_expression_x_bz_response_variant_cat.jpg`
- Scatter plot of `tbb-1` expression vs. normalized ABZ response, colored by `ben-1` variant category: `tbb1_expression_x_bz_response_variant_cat_meta.jpg`
- Scatter plot of `tbb-2` expression vs. normalized ABZ response, colored by `ben-1` variant category: `tbb2_expression_x_bz_response_variant_cat_meta.jpg`
- Faceted scatter plot of `ben-1` expression vs. ABZ response by assay, colored by `ben-1` variant category: `ben1_expression_x_bz_response_variant_cat_meta_assay.jpg`
- Faceted scatter plot of `tbb-1` expression vs. ABZ response by assay, colored by `ben-1` variant category: `tbb1_expression_x_bz_response_variant_cat_meta_assay.jpg`
- Faceted scatter plot of `tbb-2` expression vs. ABZ response by assay, colored by `ben-1` variant category: `tbb2_expression_x_bz_response_variant_cat_meta_assay.jpg`
- Boxplot of `ben-1` expression by `ben-1` variant category: `ben1_var_cat_x_ben1_exp.jpg`
- Boxplot of `tbb-1` expression by `ben-1` variant category: `ben1_var_cat_x_tbb1_exp.jpg`
- Boxplot of `tbb-2` expression by `ben-1` variant category: `ben1_var_cat_x_tbb2_exp.jpg`
- Boxplot of normalized ABZ response by `ben-1` variant category: `ben1_var_cat_x_bz_response.jpg`
- Scatter plots of `ben-1` expression vs. `tbb-1` and `tbb-2` expression: `beta_tub_expression_x_expression.jpg`
#### figures (for manuscript)
- Combined plot of `ben-1` expression vs. normalized ABZ response scatter plot and `ben-1` expression boxplot by `ben-1` variant category: `main_figure.jpg`
- Combined plot of `tbb-1` and `tbb-2` expression scatter plots and boxplots by `ben-1` variant category: `supplementary_figure_tbb1_tbb2.jpg`
- Combined faceted scatter plots of `ben-1`, `tbb-1`, and `tbb-2` expression vs. ABZ response by assay: `supplementary_figure_scatter_x_assay.jpg`
# Archive
## `ce_tubulin_expression_variation.R`
Script to look at natural variation in beta-tubulin expression
## `expression_thresholds.R`
Script to look at expression levels of beta-tubulin in resistant and sensitive strains using normalized phenotypes
### inputs:
- path to normalized phenotypes for common strains (from `normalize_phenos.R`)
### outputs:
- "{figure_out_dir}/expression_distribution_w_thresholds.jpg"