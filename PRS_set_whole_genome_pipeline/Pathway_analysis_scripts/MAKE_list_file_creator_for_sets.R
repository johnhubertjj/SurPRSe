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
Gene_regions <- args[8]
chromosomes_to_parse <- as.numeric(args[c(9:length(args))])
print(chromosomes_to_parse)

output_directory <- paste0("./", Training_name, "_", Validation_name, "_output/")

#Validation_datatable_bed_file <- "ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr"
#Training_name <- "CLOZUK_PGC2noclo"
#Validation_name <- "ALSPAC"
#Pathway_directory <- paste0("./", Training_name, "_", Validation_name, "_output/Pathways/")
#current_pathway <- "Morphology_schizophrenia_related"
#chromosomes_to_parse <- seq(1,22)


# Read in both of the error files to check which pathways you want to analyse
unread_pathways_one <- read.delim(paste0(Pathway_directory,"/MAGMA_empty_files_after_analysis_",Gene_regions,".txt"), header = F, stringsAsFactors = F)
unread_pathways_two <- fread(paste0(Pathway_directory,"/Pathways_analysis_empty_pathways_info_file.txt"), header = F, stringsAsFactors = F)

# For the moment I create Dummy files so that R will always recognise that file instead of searching to check if the file exits.
# That way there is less chance of problems for different OS and I have more control within R
if (nrow(unread_pathways_one) == 1 & nrow(unread_pathways_two) == 1){
  new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))
  
  for (i in 1:length(chromosomes_to_parse)){
    new_data_frame[i,1] <- paste0(Pathway_directory,"/",current_pathway,"/chromosome_",chromosomes_to_parse[i],"_",Validation_name,"_",Training_name,"_",current_pathway,"_",Gene_regions,"_gene_regions")
  }
  
  filename <- paste0(Pathway_directory,"/",current_pathway,"/make_full_plink_",current_pathway,"_",Gene_regions,".txt")
  write.table (new_data_frame, file = filename, col.names = F, row.names = F, quote = F)
  quit()
}

unread_pathways_one <- as.data.table(unread_pathways_one)
unread_pathways_two <- as.data.table(unread_pathways_two)
cat("Some chromosomes within the pathway do not have any SNPs contained within them for ", current_pathway," these will be ignored")

if (nrow(unread_pathways_one) == 1 & nrow(unread_pathways_two) != 1){
  unread_pathways_two <- unread_pathways_two[-1]
  rm(unread_pathways_one)
  cat("No empty chromosomes found in R analysis for ", current_pathway,"but some found in MAGMA analysis")
  
  # Print out a new make file with the relevant chromosomes removed
  setnames(unread_pathways_two, c("Pathway","chromosome"))
  unread_pathways_total <- unread_pathways_two
  
} else if (nrow(unread_pathways_one) != 1 & nrow(unread_pathways_two) == 1){
    unread_pathways_one <- unread_pathways_one[-1]
    rm(unread_pathways_two)
    cat("No empty chromosomes found in MAGMA analysis for ", current_pathway,"but some found in R analysis")
    
    # Print out a new make file with the relevant chromosomes removed
    setnames(unread_pathways_one, c("Pathway","chromosome"))
    unread_pathways_total <- unread_pathways_one
    
} else if (nrow(unread_pathways_one) != 1 & nrow(unread_pathways_two) != 1){
    # Remove the top line stating what the file is now that you know that the file exists
    unread_pathways_one <- unread_pathways_one[-1]
    unread_pathways_two <- unread_pathways_two[-1]

    # Print out a new make file with the relevant chromosomes removed
    setnames(unread_pathways_one, c("Pathway","chromosome"))
    setnames(unread_pathways_two, c("Pathway","chromosome"))
    unread_pathways_total <- rbind(unread_pathways_one,unread_pathways_two)
}

items_to_remove <- unread_pathways_total[which(unread_pathways_total$Pathway == current_pathway)]
dimensions_of_items_to_remove <- dim(items_to_remove)

if (dimensions_of_items_to_remove[1] != 0){
  chromosomes_not_involved <- as.vector(items_to_remove$chromosome)
  chromosomes_to_parse <- chromosomes_to_parse[! chromosomes_to_parse %in% chromosomes_not_involved]
}

if(length(chromosomes_to_parse) == 0){
  cat("No SNPs found for pathway ",current_pathway," analysis will stop")
  write(x = current_pathway, file = "Pathways_with_no_SNPs.txt",append = T)
  stop()
}

new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))

for (i in 1:length(chromosomes_to_parse)){
  new_data_frame[i,1] <- paste0(Pathway_directory,"/",current_pathway,"/chromosome_",chromosomes_to_parse[i],"_",Validation_name,"_",Training_name,"_",current_pathway,"_", Gene_regions,"_gene_regions")
}

filename <- paste0(Pathway_directory,"/",current_pathway,"/make_full_plink_",current_pathway,"_", Gene_regions, ".txt")
write.table (new_data_frame, file = filename, col.names = F, row.names = F, quote = F)
