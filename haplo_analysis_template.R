library(tidyverse)
library(dplyr)
library(data.table)
library(haplo.stats)
library(ggplot2)
library(data.table)


setwd("")
meta1 <- read.csv("LOAD IN META FILE") #CAD cases, PCs, sex

#If your CAD case/control is not coded as 1/0, you will need to do that
#CASE = 1
#Control = 0
meta <- meta %>% 
  mutate(CAD_coded = ifelse("FILL IN"))


#Load in PED file for 6 lead SNPs. Create genotype matrix
df <- fread("9p21_6snps.ped") 
geno <- df[,7:18] #Grab all genotype columns. Should be the last 12 columns (6 snps * 2 alleles)
geno <- data.frame(geno)
rownames(geno) <- df$V1


#Subset meta to only include individuals that we have haplotypes/genotypes for...
#Then order meta based on geno
meta <- meta[meta$eid %in% rownames(geno),]
geno <- geno[rownames(geno) %in% meta$eid,]
#Order meta to match geno matrix order
meta2 <- meta[match(rownames(geno), as.character(meta$eid)),]
#check they are in same order
identical(rownames(geno), as.character(meta2$eid)) #TRUE

geno.desc <- summaryGeno(geno, miss.val=c(0,NA)) #Returns info how many loci/genotypes are missing for each individiual
full.cases <- geno.desc$missing0 == 6 #All sites are full. This is stringent. Haplo.stat can handle individuals with missing SNPs but we choose to ignore these for now

#Remove individuals with missing genotypes
geno2 <- geno[full.cases,]
meta3 <- meta2[full.cases,]

#IGNORE 45 to 47 if you are NOT dealing with multiple populations
#Subset to ______ individuals
african_index <- meta3$ancestry == "FILL IN ANCESTRY"
geno3 <- geno2[african_index,] #subset geno
meta3 <- meta3[african_index,] #subst meta

#Pass genotype matrix to setupGeno. 
geno.glm <- setupGeno(geno3)

#Set up data frame for haplo.glm. Pass geno.glm, and any meta data for GLM
glm_df <- data.frame(geno.glm,
                     CAD=unlist(meta3$CAD_coded),
                     sex=unlist(meta3$sex),
                     PC1=unlist(meta3$PC1),
                     PC2=unlist(meta3$PC2),
                     PC3=unlist(meta3$PC3),
                     PC4=unlist(meta3$PC4),
                     PC5=unlist(meta3$PC5),
                     PC6=unlist(meta3$PC6),
                     PC7=unlist(meta3$PC7),
                     PC8=unlist(meta3$PC8),
                     PC9=unlist(meta3$PC9),
                     PC10=unlist(meta3$PC10)
)



###################################
#We will run haplo.glm to extract haplotype frequencies as well as conduct an CAD association 
###################################

# Haplo.glm has a few quirks that make it a little difficult to run:
# In order to run haplo.glm, a reference haplotype must be used to calculate Beta/ORs.
# We want to run haplo.glm twice, using a different reference/base each time. AGTTCA and AACATT
# We can choose the base by using the "haplo.base" parameter in haplo.glm.control, however
# this parameter only takes an index and not a literal string (e.g "AGTTCA").
# In order to know the index of the base that we want, we need to first run haplo.glm
# WITHOUT passing haplo.base and then look at the output to identify index associated with our desired reference.
#Once we identify the index number, we can rerun haplo.glm and pass that index.
#The index does not change across runs so this works just fine

# Code below to run without passing specific haplo.base:


locus.label = c("rs10811650","rs10757269","rs1537370","rs4007642","rs1537375","rs1333046")

#haplo.base not passed
hap_assoc_1 <- haplo.glm(CAD ~ geno.glm + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                       family = binomial,
                       data=glm_df,
                       locus.label = locus.label,
                       control = haplo.glm.control(em.c=haplo.em.control(min.posterior=1e-04,max.iter=20000,max.haps.limit=1e9,n.try=30,insert.batch.size=4),haplo.effect = "add",haplo.freq.min = 0.001))

#Now we look at output. Specifically look at the table labeled "Haplotypes:"
hap_assoc_1

#Each row is a haplotype labeled as geno.glm.#. The last row is labeled "haplo.base".
#The haplotype associated with "haplo.base" was the base chosen by haplo.glm because we did not pass anyting 
#If the the current base is one of the two haplotypes that we want (AGTTCA or AACATT),
#then great. These are the results that you want for one of our references. 
#While we are here, we can look at the other haplotypes in the table,
#and identify the other haplotype (AGTTCA or AACATT) and write down the index associated with it. 
#Then we will rerun hapo.glm and pass the index that you recorded.
#if neither of our desired haplotypes were haplo.base then we will have to recode the index for both of them
#and run haplo.glm two more times 

#Example: AGTTCA was labeled as haplo.base in hap_assoc_1, so that analysis is done
#Furthermore, AACATT was labeled as geno.glm.11 in hap_assoc_1, so we will now run the below code with haplo.base=11

hap_assoc_2 <- haplo.glm(CAD ~ geno.glm + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                       family = binomial,
                       data=glm_df,
                       locus.label = locus.label,
                       control = haplo.glm.control(em.c=haplo.em.control(min.posterior=1e-04,max.iter=20000,max.haps.limit=1e9,n.try=30,insert.batch.size=4),haplo.effect = "add",haplo.freq.min = 0.001,haplo.base=11))

save(hap_assoc_1,file="hap_assoc_FM_FILL.IN.BASE.RData")
save(hap_assoc_2,file="hap_assoc_FM_FILL.IN.BASE.RData")


#PLOT HAPLOTYPE FREQUENCY PLOT
#We can use the output from hap_assoc_1 or hap_assoc_2 to produce frequency plot


temp <- data.frame(freq = hap_assoc_1['haplo.freq'],hap_assoc_1['haplo.unique'])
temp$haplotype <- paste(temp[,2],temp[,3],temp[,4],temp[,5],temp[,6],temp[,7],sep="")
temp2 <- temp[,c("haplo.freq","haplotype")]
hap_table <- temp2[temp2$haplo.freq > 0.001,] #Remove very low frequency haplotypes for plotting purposes 


ggplot(hap_table, aes(y=haplo.freq, x=haplotype)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("Haplotype Frequency") + 
  xlab("Haplotype")


#MAKE GLM results TABLE
make_haplo_table <- function(df){
  df_haplo <- data.frame(df$haplo.freq,df$haplo.unique)
  df_haplo$haplotype <- paste(df_haplo$rs10811650,
                              df_haplo$rs10757269,
                              df_haplo$rs1537370,
                              df_haplo$rs4007642,
                              df_haplo$rs1537375,
                              df_haplo$rs1333046,sep="")
  
  df_haplo <- df_haplo[,c(8,1)]
  return(df_haplo)
}

make_summary<- function(df){
  summ <- summary(df)
  
  summ_hap <- summ$haplotypes
  summ_hap$id<- rownames(summ_hap)
  
  summ_coe <- data.frame(summ$coefficients)
  summ_coe$id <- rownames(summ_coe)
  
  hap.coef<- merge(summ_hap,summ_coe,all=T)
  hap.coef<- hap.coef[!is.na(hap.coef$hap.freq),]
  
  hap.coef$haplotype <- paste(hap.coef$rs10811650,
                              hap.coef$rs10757269,
                              hap.coef$rs1537370,
                              hap.coef$rs4007642,
                              hap.coef$rs1537375,
                              hap.coef$rs1333046,sep="")
  return(hap.coef)
  
}

hap_summary <- make_summary(hap_assoc_1)
hap_haplo <- make_haplo_table(hap_assoc_1)
hap_final <- merge(hap_haplo,hap_summary,by="haplotype",all.x = TRUE)


write.csv(hap_final,file = "glm_results_FILL.IN.BASE.csv")




