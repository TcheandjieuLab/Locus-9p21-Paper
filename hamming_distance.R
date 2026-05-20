
library(dplyr)
library(ggplot2)
library(data.table)
library("stringr") 
library(cowplot)
library(tidyr)



# Define the path to the directory with your files
path <- ""
# Assuming your haplotype files are split up in chunks. List all files with a specific pattern, e.g., ".txt"
file_list <- list.files(path = path, pattern = "*_haptable.txt", full.names = TRUE)

# Load each file, add a filename column, and combine them
combined_df <- lapply(file_list, function(file) {
  # Read in the data
  data <- read.table(file,header=TRUE)
  # Add a new column with the file name (use `basename` to keep just the file name)
  data$group <- ""
  return(data)
}) %>%
  bind_rows()  # Combine all data frames into one


#Load in all effect allele haplotype.
effect_allele <- read.table("all_noneffect_allele.txt",header = FALSE)
EAH <- effect_allele$V1[1]

# Function to calculate Hamming distance
hamming_distance <- function(hap1, hap2) {
  sum(strsplit(hap1, "")[[1]] != strsplit(hap2, "")[[1]])
}

#For each individual, calculate the hamming distance to EAH for both hap1 and hap2
combined_df <- combined_df %>%
  mutate(hamming_distance_hap1 = sapply(hap1, hamming_distance, hap2 = EAH)) %>%
  mutate(hamming_distance_hap2 = sapply(hap2, hamming_distance, hap2 = EAH))

#Calculate dosage, sum up hamming distance values for hap1 and hap2
combined_df$dosage <- combined_df$hamming_distance_hap1 + combined_df$hamming_distance_hap2

#Write hamming distance values.
outdf <- combined_df[,c("hamming_distance_hap1","hamming_distance_hap2","dosage")]
write.table(outdf,file="hamming_table.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


##### PLOTTING ##########

create_plot <- function(data, ancestry) {
  temp <- data[data$group==ancestry,]
  
  freq_df <- temp %>%
    count(hamming_distance_hap1) %>%
    rename(frequency = n) %>%
    mutate(proportion = frequency/dim(temp)[1])
  
  ggplot(freq_df, aes(x = dosage, y = proportion)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    labs(x = "Hamming Distance ", y = "Proportion of Individuals") +
    ylim(0,.12) + 
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10),  # Adjust text size as needed
          axis.text.y = element_text(size = 10)) + 
    ggtitle(ancestry)

}

plot1 <- create_plot(combined_df, "")

ggsave("hamming_plot_.pdf", plot = final_plot, width = 12, height = 10)

################# Haplotype association ##############
#Load in meta data that includes CAD case/control, PCs and sex
meta <- read.table("")
#Load in related individuals to remove from analysis 
related <- read.table("")
#remove related
meta <- meta[!meta$IID %in% related$V1,]
#Merge meta with hamming distance df
df <- merge(meta,combined_df,by="IID")


########
## Association with hap1 and hap2 seperately using CONTINUOUS hamming distance

#Using hamming distance for hap1. Run second time with age
model_hap1 <- glm(CAD ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + hamming_distance_hap1,
             data = df,
             family = binomial())

#Using hamming distance for hap2. Run second time with age
model_hap2 <- glm(CAD ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + hamming_distance_hap2,
             data = df,
             family = binomial())


#Process results and save 
coef_hap1 <- summary(model_hap1)$coefficients
hap1_hamming <- coef_hap1[grep("hamming_distance_hap1", rownames(coef_hap1)), , drop = FALSE]
hap1_hamming <- as.data.frame(hap1_hamming)
hap1_hamming$hap <- "hap1"

coef_hap2 <- summary(model_hap2)$coefficients
hap2_hamming <- coef_hap2[grep("hamming_distance_hap2", rownames(coef_hap2)), , drop = FALSE]
hap2_hamming <- as.data.frame(hap2_hamming)
hap2_hamming$hap <- "hap2"

#combine hap1 and hap2 results
hamminggroup <- rbind(hap1_hamming,hap2_hamming)
write.table(hamminggroup,file="hamming_seperated_continuous.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")


########
## Association with hap1 and hap2 seperately using CATEGORICAL hamming distance

#Discretize hamming distance 
df$hamming_group_hap1 <- cut(df$hamming_distance_hap1, 
                        breaks = 4, 
                        labels = c("Low", "Medium", "High", "Very High"), 
                        include.lowest = TRUE)

df$hamming_group_hap2 <- cut(df$hamming_distance_hap2, 
                        breaks = 4, 
                        labels = c("Low", "Medium", "High", "Very High"), 
                        include.lowest = TRUE)

#Using hamming group for hap1  
model_hap1group <- glm(CAD ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + hamming_group_hap1,
             data = df,
             family = binomial())

#Using hamming group for hap2 
model_hap2group <- glm(CAD ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + hamming_group_hap2,
             data = df,
             family = binomial())

cat_hap1 <- summary(model_hap1group)$coefficients
hap1_cat <- cat_hap1[grep("hamming_group_hap1", rownames(cat_hap1)), , drop = FALSE]
hap1_cat <- as.data.frame(hap1_cat)
hap1_cat$hap <- "hap1"

cat_hap2 <- summary(model_hap2group)$coefficients
hap2_cat <- cat_hap2[grep("hamming_group_hap2", rownames(cat_hap2)), , drop = FALSE]
hap2_cat <- as.data.frame(hap2_cat)
hap2_cat$hap <- "hap2"

catgroup <- rbind(hap1_cat,hap2_cat)
write.table(catgroup,file="hamming_seperated_categorical.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")



########
## Association using dosage CONTINUOUS

model_dosage <- glm(CAD ~ PC1 + PC2 + PC3 + PC4 + PC5 + Sex + dosage,
                       data = df,
                       family = binomial())

coef_dosage <- summary(model_dosage)$coefficients
dosage_hamming <- coef_dosage[grep("dosage", rownames(coef_dosage)), , drop = FALSE]
dosage_hamming <- as.data.frame(dosage_hamming)
write.table(dosage_hamming,file="hamming_dosage_continuous.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")

########
## Association using dosage CATEGORICAL

df$dosage_categorical <- cut(df$dosage, 
                                      breaks = 4, 
                                      labels = c("Low", "Medium", "High", "Very High"), 
                                      include.lowest = TRUE)

df$dosage_categorical <- relevel(factor(df$dosage_categorical), ref = "Low")

model_dosagecat <- glm(CAD ~ PC1 + PC2 + PC3 + PC4 + PC5 + Sex + dosage_categorical,
                    data = df,
                    family = binomial())

coef_dosagecat <- summary(model_dosagecat)$coefficients
dosagecat_hamming <- coef_dosagecat[grep("dosage_categorical", rownames(coef_dosagecat)), , drop = FALSE]
dosagecat_hamming < as.data.frame(dosagecat_hamming)
write.table(dosagecat_hamming,file="hamming_dosage_categorical.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")





