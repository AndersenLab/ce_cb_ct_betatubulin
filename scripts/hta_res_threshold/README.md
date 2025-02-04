Scripts to set the res threshold for HTA data collected in paper.

## HTA Resistance Threshold Script

The `hta_res_threshold.R` script is designed to calculate resistance thresholds for HTA data included in the paper. It identifies strains that are resistant based on a specified phenotype column and a reference strain.


### Usage

1. **Input Data**: The script reads HTA data from CSV files for different species (`C. elegans`, `C. briggsae`, `C. tropicalis`).
2. **Calculate Resistance Thresholds**: It calculates resistance thresholds for each species using the `get_res_strains_ref` function.
3. **Output Data**: The results are saved to a CSV file in the `tables/hta_res_threshold` directory.


