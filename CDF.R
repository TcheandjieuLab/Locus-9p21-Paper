#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
Usage:
Rscript plot_maf_cdf.R <sig_snps.txt> <input_config.tsv> <out_prefix>

input_config.tsv columns:
group   file    id_col  freq_col sep

Example:
AFR_Homozygous  afr_homo.tsv ID A1_FREQ \\t
South_Asian     sas.tsv.gz   SNP.chrps EAF \\t
")
}

sig_file <- args[1]
config_file <- args[2]
out_prefix <- args[3]

# Optional region filter
region_start <- ifelse(length(args) >= 4, as.numeric(args[4]), 21940000)
region_end   <- ifelse(length(args) >= 5, as.numeric(args[5]), Inf)

# -----------------------------
# Read significant SNP list
# -----------------------------
sig_snps <- fread(sig_file, header = FALSE)$V1

sig_df <- data.frame(
  id = sig_snps,
  pos = as.numeric(sub("^\\d+:(\\d+):.+", "\\1", sig_snps))
) %>%
  filter(!is.na(pos), pos >= region_start, pos <= region_end)

message("Number of significant SNPs in region: ", nrow(sig_df))


# Helper functions
read_assoc <- function(file, sep = "\t") {
  fread(
    file,
    sep = sep,
    header = TRUE,
    comment.char = "=",
    data.table = FALSE
  )
}

make_maf_df <- function(file, group, id_col, freq_col, sep = "\t") {
  df <- read_assoc(file, sep = sep)

  if (!id_col %in% colnames(df)) {
    stop("ID column not found in ", file, ": ", id_col)
  }

  if (!freq_col %in% colnames(df)) {
    stop("Frequency column not found in ", file, ": ", freq_col)
  }

  df <- df[df[[id_col]] %in% sig_df$id, ]

  df$MAF <- ifelse(
    df[[freq_col]] <= 0.5,
    df[[freq_col]],
    1 - df[[freq_col]]
  )

  data.frame(
    group = group,
    SNP = df[[id_col]],
    MAF = df$MAF
  )
}

cumulative_maf <- function(df, thresholds) {
  df <- df %>% filter(!is.na(MAF))

  data.frame(
    threshold = thresholds,
    cumulative_proportion = sapply(
      thresholds,
      function(t) mean(df$MAF <= t)
    )
  )
}

# Read input config
config <- fread(config_file, data.table = FALSE)

required_cols <- c("group", "file", "id_col", "freq_col", "sep")
missing_cols <- setdiff(required_cols, colnames(config))

if (length(missing_cols) > 0) {
  stop("Missing required columns in config: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Process each population/group
# -----------------------------
maf_all <- bind_rows(
  lapply(seq_len(nrow(config)), function(i) {
    make_maf_df(
      file = config$file[i],
      group = config$group[i],
      id_col = config$id_col[i],
      freq_col = config$freq_col[i],
      sep = config$sep[i]
    )
  })
)

maf_all$group <- gsub("_", " ", maf_all$group)

fwrite(maf_all, paste0(out_prefix, ".maf_values.tsv"), sep = "\t")

# -----------------------------
# Compute CDF
# -----------------------------
thresholds <- seq(0, 0.5, by = 0.001)

cdf_all <- maf_all %>%
  group_by(group) %>%
  group_modify(~ cumulative_maf(.x, thresholds)) %>%
  ungroup()

fwrite(cdf_all, paste0(out_prefix, ".maf_cdf.tsv"), sep = "\t")

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(
  cdf_all,
  aes(
    x = threshold * 100,
    y = cumulative_proportion * 100,
    color = group
  )
) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "MAF Threshold (%)",
    y = "Cumulative % of Significant SNPs",
    title = paste0(
      "MAF Distribution of Significant SNPs",
      " (", region_start / 1e6, " - ",
      ifelse(is.infinite(region_end), "end", region_end / 1e6),
      " Mb)"
    ),
    color = "Population"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  paste0(out_prefix, ".maf_cdf.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  paste0(out_prefix, ".maf_cdf.pdf"),
  p,
  width = 8,
  height = 6
)