#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

### Start Timer
ptm <- proc.time()

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

# specify the different input tables #
validation_set_full_name_without_chromosome <- args[3]
training_set_name <- args[4]
validation_set_name <- args[5]
Summary_analysis <- args[6]
Chromosomes_to_analyse <- as.numeric(args[c(7:length(args))])

# insert library
library(data.table)

if(Summary_analysis == FALSE){
## Read in the latest Summary stats tables after converting to one table
for (i in Chromosomes_to_analyse){
  assign(paste0("Training_table", i), fread(paste0("./",training_set_name,"_",validation_set_name,"_extrainfo/Summary_SNP_info_for_QC_of_",validation_set_name,"_",training_set_name,"_chr",i,".txt"), header = T),envir = .GlobalEnv)
}
}else{
  for (i in Chromosomes_to_analyse){
    assign(paste0("Training_table", i), fread(paste0("./",training_set_name,"_",validation_set_name,"_extrainfo/SUMMARY_USING_MAF_genotype_Summary_SNP_info_for_QC_of_",validation_set_name,"_",training_set_name,"_chr",i,".txt"), header = T),envir = .GlobalEnv)
  }
}

l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_analyse) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}
combined_final_table <- rbindlist(l)


for (i in Chromosomes_to_analyse){
  assign(paste0("Training_table", i), fread(paste0("./",training_set_name,"_",validation_set_name,"_output/",validation_set_full_name_without_chromosome,"_chr",i,"_consensus_with_",training_set_name,"_flipped_alleles_no_duplicates.bim"), header = T),envir = .GlobalEnv)
}

l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_analyse) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}
last_row_of_summary <- rbindlist(l)

summary_row <- data.frame(stage = 9, x = paste0("Number of SNPs after MAF, HWE and genotype missingness in ",validation_set_name), y = nrow(last_row_of_summary),stringsAsFactors = F)


# Final Summary information
combined_final_table <- data.table(combined_final_table, key="stage")

final_table <- combined_final_table[, list(x=unique(x), 
          y=sum(y)), by=stage]

final_table <- rbind(final_table, summary_row)

write.table(final_table, file = paste0("SNPs_lost_to_allele_flipping_and_identifiers_",training_set_name,"_",validation_set_name,".txt"), col.names = T, row.names = F, quote = F, sep = "\t")

