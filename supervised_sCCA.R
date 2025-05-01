#Supervised sparse CCA

#SPARSE SUPERVISED METHOD
library(PMA)
library(tidyverse)

#load all necessary data
metadata <- read.csv("/Users/elise/Downloads/STAT_437/analysis_all/metadata.csv")
exp_data <- read.csv("/Users/elise/Downloads/STAT_437/analysis_all/collapsed_expression.csv", row.names = 1)
copy_num_data <- read.csv("/Users/elise/Downloads/STAT_437/analysis_all/collapsed_CGH.csv", row.names = 1)
copy_num_chromosomes <- read.csv("/Users/elise/Downloads/STAT_437/collapsed_CGH_chr.csv")
exp_chromosomes <- read.csv("/Users/elise/Downloads/STAT_437/analysis_all/exp_chromosomes.csv")

# map gene expression and copy number to the samples in metadata
sample_map <- metadata %>%
  select(geo_accession_exp, geo_accession_copynum, Final_Microarray_Diagnosis, Follow_Up_Status) %>%
  filter(geo_accession_exp %in% colnames(exp_data) &
           geo_accession_copynum %in% colnames(copy_num_data))

# match to samples
exp_data_matched <- exp_data[, sample_map$geo_accession_exp]
copy_num_data_matched <- copy_num_data[, sample_map$geo_accession_copynum]
exp_data_matched <- apply(exp_data_matched, 2, as.numeric)
copy_num_data_matched <- apply(copy_num_data_matched, 2, as.numeric)

# need to have samples as rows
x <- t(exp_data_matched)
z <- t(copy_num_data_matched)

# add back in the gene and copy number location names
colnames(x) <- rownames(exp_data)
colnames(z) <- rownames(copy_num_data)

#remove columns containing NAs
bad_cols <- which(colSums(!is.finite(z)) > 0)
if(length(bad_cols) > 0) {
  z <- z[, -bad_cols]
}

#create vector of disease subtype associated with each sample
subtype_by_sample <- sample_map$Final_Microarray_Diagnosis

#find the genes most associated with subtype using ANOVA
y <- as.numeric(as.factor(subtype_by_sample))
f_exp <- apply(x, 2, function(feature) summary(aov(feature ~ y))[[1]][["Pr(>F)"]][1])
sort(f_exp, decreasing = FALSE)[1:30]
sig_genes <- f_exp[f_exp < 0.05]

#find the copy number locations most associated with subtype using ANOVA
f_copynum <- apply(z, 2, function(feature) summary(aov(feature ~ y))[[1]][["Pr(>F)"]][1])
sig_copy <- f_copynum[f_copynum < 0.05] 

#filter datasets down to just these features
x_filt <- x[, names(sig_genes)]
z_filt <- z[, names(sig_copy)]

#scale the data
x_scale <- scale(x_filt, scale = TRUE, center = TRUE)
z_scale <- scale(z_filt, scale = TRUE, center = TRUE)

#create a dataframe to store results
results_sup_df <- data.frame(Chromosome = c("Chromosome"), Subtype = 0)

# Subset copy number down to one chromosome at a time and perform supervised sCCA
for (i in 1:22){
  chrom_i <- subset(copy_num_chromosomes, CHROMOSOME == paste0("chr", i))$group_assign
  
  z_i <- z_scale[, colnames(z_scale) %in% chrom_i]
  x_i <- x_scale
  
  z_i <- z_i[, colSums(z_i) != 0]
  x_i <- x_i[, colSums(x_i) != 0]
  
  perm.out.sup.subtype <- CCA.permute(x_i, z_i,
                                      typex = "standard",
                                      typez = "ordered",
                                      nperms = 10,
                                      outcome = "quantitative",
                                      y = y)
  
  #store results in dataframe
  results_sup_df <- rbind(results_sup_df,
                          data.frame(Chromosome = as.character(i),
                                     Subtype = round(mean(perm.out.sup.subtype$pvals), 4)))
}

#save the results
write.csv(format(results_sup_df, digits = 4, nsmall = 4), "/Users/elise/Downloads/CCA_results_supervised_sub.csv")
