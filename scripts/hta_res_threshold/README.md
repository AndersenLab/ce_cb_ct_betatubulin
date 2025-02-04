Scripts to set the res threshold for HTA data collected in paper.

## HTA Resistance Threshold Script

The `hta_res_threshold.R` script is designed to calculate resistance thresholds for HTA (High Throughput Assay) data. It identifies strains that are resistant based on a specified phenotype column and a reference strain.

### Key Functions

- **get_res_strains_ref**: Calculates the resistance threshold and identifies resistant strains.
  - **Parameters**:
    - `exp_summary_df`: Data frame containing experimental summary with columns `strain` and a dynamic phenotype column.
    - `pheno_col`: Name of the phenotype column to be used for resistance calculation. The column should contain normalized animal length vales. 
    - `ref_strain`: Reference strain to calculate the resistance threshold.
    - `threshold_per`: Threshold percentage to classify strains as resistant (default is 0.50).
  - **Returns**: A list containing the minimum and maximum values of the phenotype column, resistant strains, and the calculated threshold.

### Usage

1. **Input Data**: The script reads HTA data from CSV files for different species (`C. elegans`, `C. briggsae`, `C. tropicalis`).
2. **Calculate Resistance Thresholds**: It calculates resistance thresholds for each species using the `get_res_strains_ref` function.
3. **Output Data**: The results are saved to a CSV file in the `data/proc/hta_res_threshold` directory.


