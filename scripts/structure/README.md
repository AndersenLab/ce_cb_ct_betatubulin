# Overview
Scripts to analyze structure of ABZ bound to beta-tubulins

# 1. Create PyMOL scripts to visualize variants on AF structures
```{zsh}
python scripts/structure/draw.py data/structure/variants.csv
```
This generates `scripts/structure/draw_structures.pml`

# 2. Run PyMOL to generate images
```{zsh}
pymol scripts/structure/draw_structures.pml
```
Images are saved to `data/structure/raw_images/`

# 3. Combine plots into a single figure
```{zsh}
Rscript scripts/structure/combine_missense_panels.R
```
Output: `data/structure/raw_images/missense_combined.png`