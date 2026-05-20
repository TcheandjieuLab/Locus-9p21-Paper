#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)

# ---- Command line arguments ----
option_list <- list(
  make_option(c("--sumstats"), type = "character", help = "Path to summary statistics file"),
  make_option(c("--ld_matrix_file"), type = "character", help = "Path to LD matrix file (.ld)"),
  make_option(c("--snplist_file"), type = "character", help = "Path to SNP list file (.snplist)"),
  make_option(c("--subset_snps_file"), type = "character", default = NULL,
              help = "Optional file of SNP IDs to subset (one per line) [default: all SNPs in region]"),
  make_option(c("--output"), type = "character", default = "output_plot.pdf",
              help = "Output PDF filename [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Load summary statistics ----
sumstats <- read.table(opt$sumstats, header = TRUE, comment.char = "=")
sumstats <- sumstats %>%
  mutate(
    significant = ifelse(P < 5e-8, "Significant", "Not Significant"),
    effect_size = abs(1 - OR),
    MAF = ifelse(A1_FREQ < 0.5, A1_FREQ, 1 - A1_FREQ),
    log_pvalue = -log10(P),
    MB = POS / 1e6
  ) %>%
  filter(MAF >= 0.01)

# Optional SNP subset
if (!is.null(opt$subset_snps_file)) {
  keep_snps <- readLines(opt$subset_snps_file)
  sumstats <- sumstats %>% filter(ID %in% keep_snps)
} else {
  # Default to filter by region. Subset to the region that we know we have LD info for BBJ. THis will change once Satoshi gives us the expanded set
  sumstats <- sumstats %>% filter(POS >= 21600000 & POS <= 22200000)
  #Do this for Afr Het
  #sumstats <- sumstats %>% filter(POS >= 21300000 & POS <= 22200000)
}

# Identify SNPs by ascending p-value
sorted_snps <- sumstats %>% arrange(P) %>% pull(ID)

# ---- Load LD matrix and SNP list ----
ld <- read.table(opt$ld_matrix_file, header = FALSE)
ld_snps <- read.table(opt$snplist_file, header = FALSE)

rownames(ld) <- ld_snps$V1
colnames(ld) <- ld_snps$V1

# Square LD if it contains negatives
if (any(ld < 0, na.rm = TRUE)) {
  ld <- ld^2
}

# Get lead SNP (lowest p) that is present in LD
lead_snp <- sorted_snps[which(sorted_snps %in% rownames(ld))[1]]
if (is.na(lead_snp)) stop("No lead SNP found in LD matrix.")

# Add LD scores to sumstats
sumstats <- sumstats %>%
  mutate(LD_score = map_dbl(ID, ~ ifelse(.x %in% colnames(ld), ld[lead_snp, .x], 0)))

# ---- Plot ----
muted_colors <- scale_color_gradientn(
  colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
  na.value = "grey"
)


# Set y-axis max to at least 9
y_max <- max(10, ceiling(max(sumstats$log_pvalue, na.rm = TRUE)))

p <- ggplot(sumstats, aes(x = MB, y = log_pvalue, color = LD_score)) +
  geom_point(size = 3) +
  geom_point(data = sumstats %>% filter(ID == lead_snp),
             aes(x = MB, y = log_pvalue),
             shape = 23, size = 5, fill = "black", color = "black") +
  muted_colors +
  scale_x_continuous(labels = scales::comma_format(scale = 1)) +
  scale_y_continuous(limits = c(0, y_max)) +
  theme_minimal() +
  labs(
    x = "Genomic Position (Mb)",
    y = expression(-log[10](italic(P))),
    color = "LD"
  ) +
  theme(text = element_text(size = 10, face = "bold"),
        legend.position = "right") +
  geom_vline(xintercept = c(21.94, 22.13), color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 22.07, color = "black", linetype = "dashed", size = 0.7) + 
theme(
  axis.title.x = element_blank(),         # removes x-axis label
  axis.text.x = element_text(size = 14),  # enlarges x-axis tick labels
  axis.title.y = element_text(size = 14), # keeps y-axis label large
  axis.text.y = element_text(size = 14),  # optional: enlarge y-axis tick labels
  text = element_text(face = "bold"),
  legend.position = "right"
)

# ---- Save ----
ggsave(opt$output, plot = p, width = 11, height = 3)


