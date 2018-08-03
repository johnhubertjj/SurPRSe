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
Training_name <- args[3]
Validation_name <- args[4]
Validation_set_name <- args[5]
chromosomes_to_parse <- as.numeric(args[c(6:length(args))])
print(chromosomes_to_parse)

new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))


for (i in chromosomes_to_parse){
  new_data_frame[i,1] <- paste0("./",Training_name,"_",Validation_name,"_output/",Validation_set_name,i,"_consensus_with_", Training_name,"_flipped_alleles_no_duplicates")
}

filename <- paste0("./",Training_name,"_",Validation_name,"_output/Make_file_",Validation_name,"_consensus_with_", Training_name,"_flipped_alleles_no_duplicates_full_genome_before_clumping.txt")
write.table (new_data_frame, file = filename, col.names = F, row.names = F,quote = F)
