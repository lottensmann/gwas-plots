#!/usr/bin/env Rscript

# =============================================================================
# Manhattan Plot from GWAS summary statistics
#
# Description:
#   This script reads a GWAS summary statistics file for a given trait,
#   and generates a Manhattan plot.
#
# Usage:
#   ./manhattan_plot.R <trait>
#
# Arguments:
#   <trait>    Name of the trait; used to locate the input file and name the output.
#              Input file should be named "<trait>_gwas.txt.gz"
#
# Input:
#   A gzipped text file with the following required columns:
#     - CHROM   (integer): Chromosome number
#     - Pos     (integer): Basepair position
#     - logP    (numeric): -log10(p-value)
#     - lead    (boolean): TRUE if the variant is a lead variant
#     - gene    (string):  Gene name annotated to the variant (required for lead or novel variants)
#     - status  (string): "known" or "novel"
#     - variant (string):  Variant identifier
#
# Output:
#   A Manhattan plot saved as a PDF named "<trait>_manhattan_plot.pdf"
#
# Dependencies:
#   data.table, ggplot2, ggrepel, dplyr
#
# Notes:
#   - Significance threshold set at logP > 7.30103 (P < 5e-8)
#   - logP values capped at 100 for visualization
#   - Known lead variants are colored orange, novel variants are colored red
#   - Lead and novel variants are annotated by variant identifier and gene name
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

#extract trait name from argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide a trait name. The GWAS summary statistics should be located in a file named trait_gwas.txt.gz")
trait <- args[1]

# Read in GWAS summary statistics file with required columns CHROM, Pos, logP, lead (boolean), gene, status (known or novel)
gwas <- fread(
  file = paste0(trait, "_gwas.txt.gz"),
  sep = " ",
  header = TRUE,
  data.table = FALSE
)

# Subset significant and sampled non-significant variants to avoid large image size
sig.dat <- filter(gwas, logP > 7.30103)
notsig.dat <- gwas %>%
  filter(logP <= 7.30103) %>%
  slice_sample(prop = 0.2)
gwas.dat <- bind_rows(sig.dat, notsig.dat)


# Obtain cumulative position for x-axis of manhattan plot

# Reorder file by chromosome and position
gwas.dat <- gwas.dat %>% 
  arrange(CHROM, Pos)

# Compute max position per chromosome
chr_sizes <- gwas.dat %>% 
  group_by(CHROM) %>%
  summarize(chr_len = max(Pos), .groups = "drop")

# Add cumulative offset
chr_sizes <- chr_sizes %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

# Add offsets to data and compute cumulative position
gwas.dat <- gwas.dat %>%
  left_join(chr_sizes, by = "CHROM") %>%
  mutate(BPcum = Pos + offset)

# Obtain position of chromosome numbers on x-axis
axis.set <- gwas.dat %>%
  group_by(CHROM) %>%
  summarise(center = mean(range(BPcum)), .groups = "drop")

# Cap y-axis showing -log10(P-value) at 100
gwas.dat <- gwas.dat %>%
  mutate(logP = pmin(logP, 100))

# Set color and shape of points, alternating chromosomes colored in different colors, mark lead and novel variants
gwas.dat <- gwas.dat %>%
  mutate(
    CHRwA = case_when(
      lead & status == "novel" ~ 4,
      lead ~ 3,
      CHROM %% 2 == 0 ~ 2,
      TRUE ~ 1
    ),
    CHRwA = factor(CHRwA, levels = 1:4)
  )


# Create and save the plot as pdf
options(bitmapType = 'cairo')

plot <- ggplot(gwas.dat, aes(x = BPcum, y = logP,
                             col = CHRwA, shape = CHRwA, size = CHRwA)) +
  geom_point() +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "black") + #horizontal line at genome-wide significance threshold
  geom_label_repel(
    data = filter(gwas.dat, lead),
    aes(label = paste(variant, gene)), #lead or novel variants are annotated with gene name
    size = 2, color = "black", nudge_y = 0.5, segment.alpha = 0.5
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = axis.set$center,
    labels = axis.set$CHROM
  ) +
  scale_y_sqrt(
    name = "-log10(P-value)",
    breaks = c(0, 5, 10, 15, 30, 60, 80, 100),
    labels = c(0, 5, 10, 15, 30, 60, 80, 100)
  ) +
  scale_color_manual(values = c("blue4", "blue", "orange", "red")) + #color of alternating chromosomes, lead and novel variants 
  scale_shape_manual(values = c(16, 16, 18, 18)) + #variants depicted as circles, lead and novel variants shown as diamonds
  scale_size_manual(values = c(2, 2, 4, 4)) + #larger size for lead and novel variants
  labs(title = paste("Manhattan plot", trait)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Save to file
ggsave(paste0(trait, "_manhattan_plot.pdf"), plot, width = 30, height = 12, units = "cm", dpi = 300)
