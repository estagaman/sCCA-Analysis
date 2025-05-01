###### Importing data to work with from GEO
library("GEOquery")
library("biomaRt")
library("WGCNA")

query_data <- getGEO("GSE11318")
    #2 files found: 
        #1) GSE11318-GPL570_series_matrix.txt.gz - [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
        #2) GSE11318-GPL6400_series_matrix.txt.gz - NCI NimbleGen Homo sapiens HG17 Whole Genome 385K Tiling Set version 1

gpl570 <- getGEO("GPL570") #get just one platform
gpl6400 <- getGEO("GPL6400") #get the other platform

head(Table(gpl570)) #viewing the metadata
head(Table(gpl6400)) #has CHROMOSOME category

#There are 17350 gene expression measurements and 386165 copy number measurements. 
    #(In the raw data set, more gene expression measurements are available. 
    #However, we limited the analysis to genes for which we knew the chromosomal location, 
        #and we averaged expression measurements for genes for which multiple measurements were available.)

#subset the genes 
#keep only the genes we know the chromosomal location for
#average together the measurement for genes for which multiple measurements are available

by_gene <- Table(gpl570)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #use the Ensembl database to assign chromosomal locations to genes

# Fetch the chromosomal locations for each gene in our data
gene_chr_map <- getBM(attributes = c("hgnc_symbol", "chromosome_name"),
                      filters = "hgnc_symbol",
                      values = by_gene$`Gene Symbol`,
                      mart = mart)

#merge the detected chromosomal information together with the gene information 
by_gene_with_chr <- merge(by_gene, gene_chr_map, by.x = "Gene Symbol", by.y = "hgnc_symbol")

unique(by_gene_with_chr$chromosome_name) #inspect the chromosome assignments

#create a list of possible chromosomes
chrom_options <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

#subset data to only genes we know the chromosomal location for
only_chrom <- subset(by_gene_with_chr, by_gene_with_chr$chromosome_name %in% chrom_options) #this cuts down to 82500 instead of 91000

#filter the data down to the genes contained in only_chrom and chromosomes in chrom_options
gene_chr_map <- subset(gene_chr_map, gene_chr_map$chromosome_name %in% chrom_options)
gene_chr_map <- subset(gene_chr_map, gene_chr_map$hgnc_symbol %in% unique(only_chrom$"Gene Symbol"))

#take out just the gene and chrom info - save as csv for analysis
gene_chrom <- only_chrom[, c("Gene Symbol", "chromosome_name")]
gene_chrom <- unique(gene_chrom)
write.csv(gene_chrom, "/home/2025/estagaman/STAT_437_code/for analysis/exp_chromosome.csv", row.names = FALSE)

###### Now, we need to average counts for probes that are from the same gene#######
#pull expression data
exp <- exprs(query_data[[1]])

#filter rows down to the IDs in gene_by_chrom
exp <- exp[only_chrom$'ID', ]

write.csv(only_chrom, "/home/2025/estagaman/STAT_437_code/intermediate_matrices/exp_metadata.csv") #save the metadata

#collapse rows from the same gene and rename
rownames(exp) <- seq(1, length(rownames(exp)))
rowGroup <- only_chrom$'Gene Symbol'
rowID <- rownames(exp)
method = 'Average'

exp_by_gene <- collapseRows(exp, rowGroup = rowGroup, rowID = rowID, method = method) #leaves 17312 genes

collapsed_expr <- exp_by_gene$datETcollapsed
collapsed_expr_df <- data.frame(Gene = rownames(collapsed_expr), collapsed_expr)
write.csv(collapsed_expr_df, "/home/2025/estagaman/STAT_437_code/intermediate_matrices/collapsed_expression.csv", row.names = FALSE)


###### Lastly, we retrieve the metadata for each dataset
meta_570 <- pData(query_data[[1]])
meta_6400 <- pData(query_data[[2]])

write.csv(meta_570, "/home/2025/estagaman/STAT_437_code/intermediate_matrices/metadata_expression.csv", row.names = TRUE)
write.csv(meta_6400, "/home/2025/estagaman/STAT_437_code/intermediate_matrices/metadata_copynumber.csv", row.names = TRUE)

copy_num <- exprs(query_data[[2]])
write.csv(copy_num, "/home/2025/estagaman/STAT_437_code/intermediate_matrices/copy_num_counts.csv", row.names = TRUE) #save the copy number counts

####importing the data back in so I can subset to one chromosome 
copy_info <- fData(query_data[[2]])
copy_counts <- exprs(query_data[[2]])

#save the copy number info for subsetting my chromosome
copy_info <- copy_info[, c("ID", "CHROMOSOME")]
write.csv(copy_info, "/home/2025/estagaman/STAT_437_code/intermediate_matrices/copy_num_chromosomes.csv", row.names = FALSE)

#subset the metadata down a little 
meta_570 <- pData(query_data[[1]])
meta_6400 <- pData(query_data[[2]])

#clean the data up
remove <- c("Tissue: ", "Disease state: ", "Individual: ", "Clinical info: Submitting diagnosis: ", "Clinical info: Final microarray diagnosis: ", "Clinical info: Follow up status: ", "Clinical info: Follow up years: ", "Clinical info: Chemotherapy: ", "Clinical info: ECOG performance status: ", "Clinical info: Stage: ", "Clinical info: LDH ratio: ")
col_names <- c("Tissue", "Disease_State", "Individual", "Diagnosis", "Final_Microarray_Diagnosis", "Follow_Up_Status", "Follow_Up_Years", "Chemotherapy", "ECOG_Performance_Status", "Cancer_Stage", "LDH Ratio")
old_col_names <- c("characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4", "characteristics_ch1.5", "characteristics_ch1.6", "characteristics_ch1.7", "characteristics_ch1.8", "characteristics_ch1.9", "characteristics_ch1.10", "characteristics_ch1.11", "characteristics_ch1.12")
for (i in 1:length(remove)) {
  # Create the new cleaned column using the new name
  meta_570[[col_names[i]]] <- gsub(remove[i], "", meta_570[[old_col_names[i]]])
  
  # Remove the old column
  meta_570[[old_col_names[i]]] <- NULL
}

exp_cols_to_keep <- c("geo_accession", col_names)
meta_570 <- meta_570[, exp_cols_to_keep]

meta_6400$geo_accession_exp <- meta_570$geo_accession

meta_6400 <- meta_6400[, c("geo_accession", "Individual:ch2")]
individual_order <- meta_570$'Individual'
meta_6400 <- meta_6400[match(individual_order, meta_6400$`Individual:ch2`), ]

meta_570$geo_accession_copynum <- meta_6400$geo_accession
meta_570$geo_accession_exp <- meta_570$geo_accession
meta_570$geo_accession <- NULL

write.csv(meta_570, "/home/2025/estagaman/STAT_437_code/chr1/metadata.csv")
