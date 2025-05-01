#collapsing the copy number data by adjacent locations

#the paper performed a collapse for computational reasons, as the original dataset contains ~300,000 copy number locations
library("data.table")
library("WGCNA")

#import copy number raw data and the chromosomes that match each location
copy_num_data <- as.data.frame(fread("/Users/elise/Downloads/STAT_437/analysis_all/copy_num_counts.csv"))
copy_num_chromosomes <- as.data.frame(fread("/Users/elise/Downloads/STAT_437/analysis_all/copy_num_chromosomes.csv"))

rownames(copy_num_data) <- copy_num_data$V1
copy_num_data$V1 <- NULL

#remove X and Y chromosomes from the data
copy_num_chromosomes <- subset(copy_num_chromosomes, !copy_num_chromosomes$CHROMOSOME %in% c("chrX", "chrY"))
copy_num_data <- copy_num_data[copy_num_chromosomes$ID, ]

#start assign each location into a group of adjacent locations
all_group_assign <- c()

for (i in 1:22){ #for each chromosome
  copy_chr <- subset(copy_num_chromosomes, copy_num_chromosomes$CHROMOSOME == paste0("chr", i)) #find the copy number locations on that chromosome
  num_groups <- ceiling(nrow(copy_chr)/10) #identify the number of groups needed if each group has a size of 10
  
  group_assign <- c() #group assignments for just that chromosome
  for (g in 1:num_groups){ #for each group
    group_num <- paste0(i, "_", g) #create a group name using chromosome # and group #
    groups <- c(rep(group_num, 10)) #repeat it 10 times
    group_assign <- c(group_assign, groups) #add it to our group assignments
  }
  group_assign <- group_assign[1:nrow(copy_chr)] #make sure our group assignment vector is only as long as our copy number information 
  
  all_group_assign <- c(all_group_assign, group_assign) #add to our group assignments vector for all chromosomes
}

copy_num_chromosomes$group_assign <- all_group_assign #assign groups for the entire copy number dataset

rownames(copy_num_data) <- seq(1, length(rownames(copy_num_data))) #change rownames 
rowGroup <- copy_num_chromosomes$group_assign #want to group by group assignment

copy_num_dt <- as.data.table(copy_num_data) #turn our data into a data table
copy_num_dt[, group := rowGroup] #group by the assignments I made
collapsed_avg <- copy_num_dt[, lapply(.SD, mean), by = group] #average the counts for each group
collapsed_avg <- as.data.frame(collapsed_avg) #convert back to a dataframe, because that is easier for me to work with
rownames(collapsed_avg) <- collapsed_avg$group #reassign rownames
collapsed_avg$group <- NULL #remove the group column

#save results as a csv
write.csv(collapsed_avg, "/Users/elise/Downloads/STAT_437/analysis_all/collapsed_CGH.csv")

copy_num_chromosomes$ID <- NULL
copynum_chromosomes <- unique(copy_num_chromosomes)

#save the chromosomes matching each group name
write.csv(copynum_chromosomes, "/Users/elise/Downloads/STAT_437/collapsed_CGH_chr.csv")
