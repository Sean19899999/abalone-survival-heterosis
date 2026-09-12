# 00_common_theme.R — Shared color palette and theme for Figure 1-4
# Unified with p2.fulllength (Figure 5-10) decisions.md conventions

library(ggplot2)
library(RColorBrewer)

# --- Genotype colors (Set1, consistent with p2.fulllength PCA) ---
col_genotype <- c(
  "DD" = "#0072B2",   # steel blue  (H. discus hannai)
  "GG" = "#D55E00",   # vermillion  (H. gigantea)
  "GD" = "#CC79A7"    # reddish purple  (Hybrid)
)

# --- Age colors ---
col_age <- c(
  "1-year" = "#FDB462",  # orange-yellow
  "2-year" = "#B3B3B3"   # gray
)

# --- Condition colors ---
col_condition <- c(
  "Control"  = "#FFFFFF",
  "Infected" = "#333333"
)

# --- Genotype labels for plots (italic Latin names) ---
genotype_labels <- c(
  "DD" = expression(italic("H. discus hannai")),
  "GG" = expression(italic("H. gigantea")),
  "GD" = "Hybrid"
)

# Plain text version - use code labels directly (DD/GG/GD)
# Latin names are documented in Figure Documentation files
genotype_labels_plain <- c(
  "DD" = "DD",
  "GG" = "GG",
  "GD" = "GD"
)

# Parseable version for scale_x_discrete with parse=TRUE
# Use: scale_x_discrete(labels = function(x) parse(text = genotype_labels_parse[x]))
genotype_labels_parse <- c(
  "DD" = "italic('H. discus hannai')",
  "GG" = "italic('H. gigantea')",
  "GD" = "'Hybrid'"
)

# --- Unified theme ---
theme_phenotype <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.title    = element_text(size = 14),
      axis.text     = element_text(size = 13),
      legend.title  = element_text(size = 13, face = "bold"),
      legend.text   = element_text(size = 12),
      strip.background = element_blank(),
      strip.text    = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin   = margin(10, 10, 10, 10)
    )
}

# --- Significance annotation helper ---
sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "ns")))
}

# --- Genotype ordering (parents first, hybrid last) ---
genotype_order <- c("DD", "GG", "GD")

# --- Map old naming conventions ---
# SS/ss → GG, SD/sd → GD, DD/dd → DD
remap_genotype <- function(x) {
  x <- toupper(x)
  x[x == "SS"] <- "GG"
  x[x == "SD"] <- "GD"
  factor(x, levels = genotype_order)
}
