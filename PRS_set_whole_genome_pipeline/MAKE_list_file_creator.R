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
chromosomes_to_parse <- as.numeric(args[c(6:length(args))])
print(chromosomes_to_parse)

new_data_frame <- data.frame(filenames=rep(0,length(chromosomes_to_parse)))

for (i in chromosomes_to_parse){
  new_data_frame[i,1] <- paste0(Validation_datatable_bed_file,i,"_",Training_name,"_CLUMPED")
}

filename <- paste0("./",Training_name,"_",Validation_name,"_output/",Validation_name,"_",Training_name,"_FULL_GENOME_CLUMPED_MAKE_LIST_INPUT.txt")
write.table (new_data_frame, file = filename, col.names = F, row.names = F,quote = F)
