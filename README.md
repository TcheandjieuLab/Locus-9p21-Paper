# Locus-9p21-Paper

Haplo stats analysis (Tcheandjieu et al 2022): haplo_analysis_template.R

Hamming distance analysis: 
1. construct_haplotype.R: construct haplotypes from phased VCF. Filters to snps in hamming_distance_snps_updated.txt. Make sure all snps are included in VCF
2. hamming_distance.R: calculate hamming distance between theoretical all noneffect allele (all_noneffect_allele.txt) and hap table from step1. Conduct associations

