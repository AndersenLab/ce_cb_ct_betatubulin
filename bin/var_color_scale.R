# Colors of high-impact variants for plots

# variants in WI
strain_var_colors <- c(
  "Deletion" = "#065275",
  "Duplication" = "cornflowerblue",
  "Frameshift" = "#88C3D0",
  "Inframe_deletion" = "#39B185",
  "Inversion" = "#9BE3BD",
  "Missense" = "#FDDE9C",
  "Splice Donor" = "#FF928E",
  "Start/stop altering" = "darkred",
  "Low ben-1 expression" = "#FFDCDE",
  "Transposon insertion" = "#620273"
)

# set meta cataegory colors
meta_cat_cols <- c(
  "SV" = "#065275",
  "Frame altering" = "#88C3D0",
  "Missense" = "#FDDE9C",
  "Start/stop altering" = "darkred",
  "No variant" = "white",
  "Low ben-1 expression" = "#FFDCDE"
)

# Define custom labels with one label italicized
custom_meta_cat_labels <- c(
  "SV" = "SV",
  "Frame altering" = "Frame altering",
  "Missense" = "Missense",
  "Start/stop altering" = "Start/stop altering",
  "No variant" = "No variant",
  "Low ben-1 expression" = expression("Low " * italic("ben-1") * " expression")
)