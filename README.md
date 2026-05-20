# Locus-9p21-Paper

Haplo stats analysis (Tcheandjieu et al 2022): haplo_analysis_template.R

Hamming distance analysis: 
1. construct_haplotype.R: construct haplotypes from phased VCF. Filters to snps in hamming_distance_snps_updated.txt. Make sure all snps are included in VCF
2. hamming_distance.R: calculate hamming distance between theoretical all noneffect allele (all_noneffect_allele.txt) and hap table from step1. Conduct associations

Generate locus zoom plots with LD file: 
1. generate_lz.R
2. Example run: Rscript generate_lz.R --sumstats association.PHENO1.glm.logistic.hybrid --ld_matrix_file 9p21.ld --snplist_file 9p21.snplist --output lz_plot.png

Calculate CDF: 
1. CDF.R
2. This script generates cumulative minor allele frequency (MAF) distribution plots for a set of significant variants across multiple populations or ancestry groups. Given a list of significant SNPs and a configurable set of association summary statistic files, the script
3. Example run: Rscript plot_maf_cdf.R all_sig_snps.txt input_config.tsv maf_distribution 21940000 22130000
