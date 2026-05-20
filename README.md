# Locus-9p21-Paper

Scripts and analyses used for the 9p21 locus paper.

---

## Haplotype Statistics Analysis (Tcheandjieu et al. 2022)

### 1. `haplo_analysis_template.R`

Implements haplotype-based association analyses following the framework described in Tcheandjieu et al. (2022).

---

## Hamming Distance Analysis

### 1. `construct_haplotype.R`

Constructs haplotypes from phased VCF files.

**Workflow:**
- Parses phased VCF input
- Filters variants to SNPs listed in `hamming_distance_snps_updated.txt`
- Outputs haplotype table for downstream analyses

**Note:** Ensure all SNPs in `hamming_distance_snps_updated.txt` are present in the input VCF.

### 2. `hamming_distance.R`

Calculates Hamming distance between observed haplotypes and a theoretical reference haplotype.

**Workflow:**
- Uses the theoretical all–non-effect allele haplotype from `all_noneffect_allele.txt`
- Computes Hamming distance relative to haplotypes generated in Step 1
- Performs association analyses using calculated Hamming distances

---

## Generate LocusZoom Plots with LD Information

### 1. `generate_lz.R`

Generates locus zoom plots using summary statistics and an LD matrix.

### Example Run

```bash
Rscript generate_lz.R \
    --sumstats association.PHENO1.glm.logistic.hybrid \
    --ld_matrix_file 9p21.ld \
    --snplist_file 9p21.snplist \
    --output lz_plot.png
```

---

## LD Heatmap Plotting

### 1. `LD_PLOTS.R`

Generates publication-quality LD heatmaps from a linkage disequilibrium (LD) matrix and a user-defined subset of variants.

**Features:**
- Loads and validates LD matrices and SNP identifiers (`chr:pos:ref:alt` format)
- Automatically converts signed correlation matrices (*r*) to \(R^2\) when negative values are detected
- Subsets the LD matrix to variants of interest and orders SNPs by genomic position
- Enforces matrix symmetry and handles missing values
- Creates lower-triangular LD heatmaps using **ggplot2**
- Supports publication-ready raster (`.png`) and vector (`.pdf`, `.svg`, `.eps`) outputs suitable for Illustrator editing

---

## Cumulative MAF Distribution Analysis

### 1. `CDF.R`

Generates cumulative minor allele frequency (MAF) distribution plots for significant variants across multiple populations or ancestry groups.

**Workflow:**
- Loads a list of significant SNPs
- Reads association summary statistic files using a configurable input table
- Filters variants to a specified genomic region
- Computes minor allele frequencies (MAF)
- Calculates cumulative proportions across MAF thresholds
- Produces publication-ready plots and tabular outputs

### Example Run

```bash
Rscript plot_maf_cdf.R \
    all_sig_snps.txt \
    input_config.tsv \
    maf_distribution \
    21940000 \
    22130000
```

---

## Iterative Conditional Analysis

### 1. `conditional_analysis.sh`

Performs iterative PLINK2 conditional association analysis within a user-specified genomic region.

**Workflow:**
- Runs PLINK2 association analysis using configurable genotype, phenotype, and covariate inputs
- Supports standard association filters (`MAF`, `MAC`, `MACH R²`, etc.)
- Automatically identifies the top associated SNP from each run (lowest p-value)
- Appends the top SNP to a running condition list
- Repeats conditional analysis for a user-defined number of iterations

**Features:**
- Region-based analysis (`CHR`, `REGION_START`, `REGION_END`)
- Automatic iterative conditioning
- Handles PLINK scientific-notation p-values correctly
