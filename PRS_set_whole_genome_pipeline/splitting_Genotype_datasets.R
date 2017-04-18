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


Validation_datatable_bim_file <- paste0(args[4],".bim")
n_of_chrms <- args[5]

# Set working directory
setwd(".")

# Read in bim file
Validation_dataset <- fread(Validation_datatable_bim_file)

for (i in n_of_chrms){
  current_matrix <- Validation_dataset[V1 == i]
  current_SNPs <- current_matrix$V2 
  write(current_SNPs, file = paste0("convert_to_chromosomes_GeM", i, ".txt"))
}