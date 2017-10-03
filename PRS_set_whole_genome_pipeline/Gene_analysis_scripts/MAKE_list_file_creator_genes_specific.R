# Making a MAKE-FILE #
### Start Timer
ptm <- proc.time()

## New environment
e <- new.env()

## Load packages ## 
library(stringr)
library(data.table)
library(plyr)

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

# specify the different input tables #
Validation_datatable_bed_file <- args[3]
Training_name <- args[4]
Validation_name <- args[5]
Gene_directory <- args[6]
Gene_regions <- args[7] # whether to include/exclude the regulatory regions of a gene
chromosomes_to_parse <- as.numeric(args[c(8:length(args))])
print(chromosomes_to_parse)

if(Gene_regions == "normal"){ 
new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))

for (i in chromosomes_to_parse){
  new_data_frame[i,1] <- paste0(Gene_directory,"/chromosome_",i,"_",Validation_name,"_",Training_name,"_normal_gene_regions")
}

filename <- paste0(Gene_directory,"/make_plink_file_", Validation_name, "_", Training_name,"_normal_gene_regions.txt")
write.table (new_data_frame, file = filename, col.names = F, row.names = F,quote = F)
}
if(Gene_regions == "extended"){
  new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))
  
  for (i in chromosomes_to_parse){
    new_data_frame[i,1] <- paste0(Gene_directory,"/chromosome_",i,"_",Validation_name,"_",Training_name,"_extended_gene_regions")
  }
  
  filename <- paste0(Gene_directory,"/make_plink_file_", Validation_name, "_", Training_name,"_extended_gene_regions.txt")
  write.table (new_data_frame, file = filename, col.names = F, row.names = F,quote = F)
}