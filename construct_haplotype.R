library(data.table)
library(stringr)
library(dplyr)


#usage Rscript construct_haplotypes.R input.vcf output.txt
  #Rscript construct_haplotypes.R hamming.vcf bbj_haptable.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).vcf", call.=FALSE)
} 

#load in 9p21 snps
CAD = "hamming_distance_snps.txt"
cad <- fread(CAD,header=FALSE)
sigsnps <- cad$V1


#Load VCF file
vcf <- vcfR::read.vcfR(args[1])

#Subset vcf to significant SNPs only
vcf <- vcf[vcf@fix[,"ID"] %in% sigsnps,] #subset to significant snps only


#Make genotype matrix
gt <- vcfR::extract_gt_tidy(vcf)
gt[c("hap1","hap2")] <- str_split_fixed(gt$gt_GT_alleles,"\\|",2)

hap1_gt <- gt[,c('Key','Indiv','hap1')]
hap2_gt <- gt[,c('Key','Indiv','hap2')]

hap1 <- tidyr::pivot_wider(
  data = hap1_gt,
  id_cols = .data$Key,
  names_from = .data$Indiv,
  values_from = .data$hap1
)

hap2 <- tidyr::pivot_wider(
  data = hap2_gt,
  id_cols = .data$Key,
  names_from = .data$Indiv,
  values_from = .data$hap2
)

extract_hap <- function(hap,col_name){
  temp_hap = hap #hap2
  hap_result = data.frame(sample=character(),hap=character())
  for (i in 2:dim(temp_hap)[2]){
    seq = paste(unlist(temp_hap[,i]),collapse="")
    sample = colnames(temp_hap)[i]
    hap_result[nrow(hap_result) + 1,] = c(sample,seq)
  }
  colnames(hap_result) <- c("sample",col_name)
  return(hap_result)
}

hap1_result <- extract_hap(hap1,"hap1")
hap2_result <- extract_hap(hap2,"hap2")

final_haps <- merge(hap1_result,hap2_result,by="sample")

write.table(final_haps,file=args[2],sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


