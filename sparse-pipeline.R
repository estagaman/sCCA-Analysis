library(PMA)
library(tidyverse)

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

bad_cols <- which(colSums(!is.finite(z)) > 0)
if(length(bad_cols) > 0) {
  z <- z[, -bad_cols]
}

subtype_by_sample <- sample_map$Final_Microarray_Diagnosis #fill this in with the subtypes from each sample, in the order of the samples
survival_by_sample <- sample_map$Follow_Up_Status

#### run it ######
results_df <- data.frame(Chromosome = c("Chromosome"), Percent_Genes_on_Chr = 0, P_val = 0, Subtype = 0, Survival = 0)

# Subset copy number down to one chromosome at a time
for (i in 1:22){
  # keep only the copy number information from that chromosome
  chrom_i <- subset(copy_num_chromosomes, CHROMOSOME == paste0("chr", i))$group_assign
  z_i <- z[, colnames(z) %in% chrom_i] #some data manipulation has to happen for this to work
  
  # scale the data
  x <- scale(x)
  z_i <- scale(z_i)
  
  x <- x[, colSums(x) != 0]
  z_i <- z_i[, colSums(z_i) != 0]
  
  # This finds the optimal sparsity level using permutations which basically randomizes our data to remove any real correlations and then tests 10 diff penalty values to see which gives the best separation between real and random results 
  perm.out <- CCA.permute(x, z_i,
                          typex = "standard", 
                          typez = "ordered", 
                          nperms = 10)  
  # basic structure for CCA 
  cca.out <- CCA(x, z_i,
                 typex = "standard",
                 typez = "ordered",
                 K = 1, 
                 penaltyx = perm.out$bestpenaltyx,
                 penaltyz = perm.out$bestpenaltyz,
                 v = perm.out$v.init)   
  print(cca.out)
  
  non_zero_indices <- which(cca.out$u != 0)
  non_zero_genes <- colnames(x)[non_zero_indices]
  non_zero_coefficients <- cca.out$u[non_zero_indices]
  genes_chr_i <- subset(exp_chromosomes, exp_chromosomes$chromosome_name == i)
  results_by_gene <- data.frame(Gene = non_zero_genes, 
                                Coefficient = non_zero_coefficients, 
                                Chromosome = non_zero_genes %in% genes_chr_i$Gene.Symbol)
  
  percent_on_chrom <- sum(results_by_gene$Chromosome)/nrow(results_by_gene)
  
  perm.out.sup.subtype <- CCA.permute(x, z = z_i,
                                      typex = "standard",   # there is also ordered type as well for which adds fused lasso penalty not sure if we need this 
                                      typez = "ordered", ##changed this to fused lasso
                                      nperms = 10,
                                      outcome = "quantitative",
                                      y = as.numeric(as.factor(subtype_by_sample))
  )
  
  #filter out samples where survival is NA
  NA_x <- subset(sample_map, is.na(sample_map$Follow_Up_Status))$geo_accession_exp
  NA_z <- subset(sample_map, is.na(sample_map$Follow_Up_Status))$geo_accession_copynum
  
  z_i <- z_i[, !colnames(z_i) %in% NA_z]
  x <- x[, !colnames(x) %in% NA_x]
  
  perm.out.sup.survival <- CCA.permute(x, z = z_i,
                                       typex = "standard",   # there is also ordered type as well for which adds fused lasso penalty not sure if we need this 
                                       typez = "ordered", ##changed this to fused lasso
                                       nperms = 10,
                                       outcome = "quantitative",
                                       y = as.numeric(as.factor(survival_by_sample))
  )
  
  results_df <- rbind(results_df,
                      data.frame(
                        Chromosome = as.character(i),
                        Percent_Genes_on_Chr = percent_on_chrom,
                        P_val = round(mean(perm.out$pvals), 4),
                        Subtype = round(mean(perm.out.sup.subtype$pvals), 4),
                        Survival = round(mean(perm.out.sup.survival$pvals), 4)
                        ))
}

#for each chromosome, we would need to somehow: 
#1) save the important results
#2) know what percentage of genes that had non-zero weights were from chromosome i 

#save the results
write.csv(format(results_df, digits = 4, nsmall = 4), "/Users/elise/Downloads/CCA_results_collapsed.csv")
