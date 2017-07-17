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
Pathway_directory <- args[6]
current_pathway <- args[7]
chromosomes_to_parse <- as.numeric(args[c(8:length(args))])
print(chromosomes_to_parse)

output_directory <- paste0("./", Training_name, "_", Validation_name, "_output/")

unread_pathways_one <- read.delim(paste0(Pathway_directory,"/MAGMA_empty_files_after_analysis.txt"), header = F, stringsAsFactors = F)
unread_pathways_two <- fread(paste0(Pathway_directory,"/Pathways_analysis_empty_pathways_info_file.txt"), header = F, stringsAsFactors = F)
unread_pathways_one <- as.data.table(unread_pathways_one)
unread_pathways_two <- as.data.table(unread_pathways_two)

setnames(unread_pathways_one, c("Pathway","chromosome"))
setnames(unread_pathways_two, c("Pathway","chromosome"))
unread_pathways_total <- rbind(unread_pathways_one,unread_pathways_two)

items_to_remove <- unread_pathways_total[which(unread_pathways_total$Pathway == current_pathway)]

if (nrow(items_to_remove > 0)){
  chromosomes_not_involved <- as.vector(items_to_remove$chromosome)
  chromosomes_to_parse <- chromosomes_to_parse[! chromosomes_to_parse %in% chromosomes_not_involved]
}

new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))

for (i in 1:length(chromosomes_to_parse)){
  new_data_frame[i,1] <- paste0(Pathway_directory,"/",current_pathway,"/chromosome_",chromosomes_to_parse[i],"_",Validation_name,"_",Training_name,"_",current_pathway)
}

filename <- paste0(Pathway_directory,"/",current_pathway,"/make_full_plink_",current_pathway,".txt")
write.table (new_data_frame, file = filename, col.names = F, row.names = F, quote = F)
