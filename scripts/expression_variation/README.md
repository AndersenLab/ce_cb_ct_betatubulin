# Plot beta-tubulin expression level x ABZ response

## `exp_plots.R`
Contains functions used to plot 

## `beta_tub_expression_plots.R`
Code to generate manuscript plots

Plots:
- Plot of ben-1 expression and ABZ response (Figure S2)
- Plot of other beta-tubulin expression (tbb-1, tbb-2, mec-7, and tbb-4) and ABZ response (Figure S3)
  - `other_beta_tubulin_abz` prefix


## 20250723 Updates
- Adding mec-7 and tbb-4 to the plot with tbb-1 and tbb-2
  - Added columns containing mec-7 and tbb-4 expression levels to the select function to generate `exp_summary_df`
  - Added columns in NA filtering
  - Updated to have function to apply theme adjustments
  - Updated output filenames to match current manuscript
  - Switch panel labels to lower case