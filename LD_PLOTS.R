#!/usr/bin/env Rscript
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(reshape2)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))
suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(svglite)))



# ----- Argument Parser -----
option_list <- list(
  make_option(c("--ld_matrix_file"), type = "character", help = "Path to LD matrix file"),
  make_option(c("--snplist_file"), type = "character", help = "Path to SNP list file"),
  make_option(c("--subset_snps_file"), type = "character", help = "File containing SNPs of interest (one per line)"),
  make_option(c("--output"), type = "character", default = NULL, help = "Optional output PNG file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----- Load Input Files -----
# Read SNP list with full IDs (chr:pos:ref:alt)
snp_df <- read.table(opt$snplist_file, header = FALSE, stringsAsFactors = FALSE)
colnames(snp_df) <- "ID"
snp_df <- snp_df %>% separate(ID, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE)

# Read LD matrix
ld_mat <- read.table(opt$ld_matrix_file, header = FALSE, row.names = snp_df$ID)
colnames(ld_mat) <- snp_df$ID
ld_mat <- as.matrix(ld_mat)

if (any(ld_mat < 0, na.rm = TRUE)) {
  message("Detected negative values in LD matrix. Squaring it to compute R²...")
  ld_mat <- ld_mat^2
}

# Read list of SNPs of interest from file
snps_of_interest <- scan(opt$subset_snps_file, what = character())

# Subset matrix

# Subset SNPs that are in both the LD matrix and the input list
subset_snps <- intersect(snps_of_interest, colnames(ld_mat))
cat(sprintf(
  "Found %d out of %d SNPs in the LD matrix (%s).\n",
  length(subset_snps),
  length(snps_of_interest),
  opt$ld_matrix_file
))
if (length(subset_snps) < 2) stop("Need at least 2 SNPs of interest present in the LD matrix.")

# Extract position and sort subset_snps by position. Do not assume an ordered snp list was past
subset_snps_ordered <- subset_snps[order(as.integer(sub("\\d+:(\\d+):.*", "\\1", subset_snps)))]


# Subset and ensure matrix is symmetric
ld_mat_sub <- ld_mat[subset_snps, subset_snps, drop = FALSE]
ld_mat_sub[upper.tri(ld_mat_sub)] <- t(ld_mat_sub)[upper.tri(ld_mat_sub)]  # enforce symmetry


# Melt for ggplot and ensure proper ordering
LD_long <- melt(ld_mat_sub, varnames = c("SNP1", "SNP2"), value.name = "LD_value") %>%
  mutate(SNP1 = factor(SNP1, levels = subset_snps_ordered),
         SNP2 = factor(SNP2, levels = subset_snps_ordered)) %>%
  filter(as.numeric(SNP1) >= as.numeric(SNP2)) %>%
  filter(!is.na(LD_value))  # ensure no NA values sneak in


# Warn if LD values were missing
missing_count <- sum(is.na(ld_mat_sub))
if (missing_count > 0) {
  cat(sprintf("Warning: %d missing LD values in submatrix (after subsetting).\n", missing_count))
}

# Label only first and last SNPs
first_last_snps <- c(head(subset_snps_ordered, 1), tail(subset_snps_ordered, 1))
diagonal_labels <- LD_long %>%
  filter(SNP1 == SNP2 & SNP1 %in% first_last_snps) %>%
  mutate(label_x = SNP2, label_y = SNP1, snp_label = as.character(SNP1))

# ----- Plot -----
p <- ggplot(LD_long, aes(x = SNP2, y = SNP1, fill = LD_value)) +
  geom_tile(color = "white", linewidth = 0) +  # force tile edges to be seamless
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, limits = c(0, 1), name = expression(R^2)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 40, 30, 10)
  ) +
  labs(fill = expression(R^2), x = "SNPs", y = "") +
  geom_text(data = diagonal_labels,
            aes(x = label_x, y = label_y, label = snp_label),
            angle = -45, hjust = 0, vjust = 1, size = 2.5,
            nudge_x = 0.8, nudge_y = -0.8) +
  coord_fixed(clip = "off")


# Output
# Handle output file format and path (supports vector formats for Illustrator)
if (!is.null(opt$output)) {
  if (grepl("\\.pdf$", opt$output)) {
    ggsave(opt$output, plot = p, width = 6, height = 5, device = cairo_pdf)
  } else if (grepl("\\.svg$", opt$output)) {
    ggsave(opt$output, plot = p, width = 6, height = 5, device = "svg")
  } else if (grepl("\\.eps$", opt$output)) {
    ggsave(opt$output, plot = p, width = 6, height = 5, device = "eps")
  } else if (grepl("\\.png$", opt$output)) {
    ggsave(opt$output, plot = p, width = 6, height = 5, dpi = 600)
  } else {
    stop("Output file must end with .pdf, .svg, .eps, or .png")
  }
} else {
  print(p)
}


