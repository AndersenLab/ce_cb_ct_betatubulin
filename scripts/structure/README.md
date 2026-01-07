# Overview 
Scripts to analyze structure of ABZ bound to beta-tubulins

# 1. Create Pymol scripts to visualize variants on AF structures
```{zsh}
python draw.py data/raw/variants.csv > draw.pml
```

# 2. Run Pymol to generate images
```{zsh}
pymol draw.pml
```

# 3. Combine plots into a single figure
```{zsh}
Rscript combine_plots.R
```