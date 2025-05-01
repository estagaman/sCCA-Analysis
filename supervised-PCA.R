# Load libraries
library(mixOmics)
library(tidyverse)
library(dplyr)
# Load data
metadata <- read.csv("/Users/usid/Downloads/metadata.csv")
exp_data <- read.csv("/Users/usid/Downloads/collapsed_expression.csv", row.names = 1)
col_CGH <- read.csv("/Users/usid/Downloads/collapsed_CGH.csv")

# Match to samples in metadata
sample_map <- metadata %>%
  dplyr::select(geo_accession_exp, geo_accession_copynum, Final_Microarray_Diagnosis) %>%
  filter(geo_accession_exp %in% colnames(exp_data) & 
           geo_accession_copynum %in% colnames(col_CGH))

exp_data_matched <- exp_data[, sample_map$geo_accession_exp]
col_CGH_matched <- col_CGH[, sample_map$geo_accession_copynum]

# Convert to numeric
exp_data_matched <- apply(exp_data_matched, 2, as.numeric)
col_CGH_matched <- apply(col_CGH_matched, 2, as.numeric)

# Transpose and scale
x_exp <- scale(t(exp_data_matched))
x_cgh <- scale(t(col_CGH_matched))

# Set matching rownames
rownames(x_exp) <- sample_map$geo_accession_exp
rownames(x_cgh) <- sample_map$geo_accession_exp  

# Define outcome labels
y_mixomics <- factor(sample_map$Final_Microarray_Diagnosis)

# Combine datasets into a list
X_list <- list(expression = x_exp, copy_number = x_cgh)

# Run supervised PCA using block plsda
block_plsda_model <- block.splsda(X = X_list,
                                  Y = y_mixomics,
                                  ncomp = 3)

# Plot results
plotIndiv(block_plsda_model,
          comp = c(1, 2),
          group = y_mixomics,
          ind.names = FALSE,
          pch = 16,
          legend = TRUE,
          title = "Supervised PCA: Expression and Copy Number")

