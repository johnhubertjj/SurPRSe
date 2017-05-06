### Start Timer
ptm <- proc.time()

## Load packages ## 
library(data.table)

## setting working directory
setwd(".")

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

Training_name_full_unseparated <- args[3]
Training_set_name <- args[4]
Training_name_full_unseparated <- "/Users/johnhubert/Documents/testing_cross_disorder/BPSCZ.bp_v_scz.results.txt"

## Reading in PGC data
## Select for CHR 22
Training_data <- fread(Training_name_full_unseparated)
Parsing_colnames <- colnames(Training_data)

## what type of files do we have?
if(length(grep("CHR", Parsing_colnames)) == 1){
Chromosome_name <- "CHR_name=TRUE"
}else{
Chromosome_name <- "CHR_name=FALSE"  
}

if(length(grep("SNP", Parsing_colnames)) == 1){
  SNP_name <- "SNP_name=TRUE"
}else{
  SNP_name <- "SNP_name=FALSE"  
}

if(length(grep("BP", Parsing_colnames)) == 1){
  BP_name <- "BP_name=TRUE"
}else{
  BP_name <- "BP_name=FALSE"  
}

if(length(grep("A1", Parsing_colnames)) == 1){
  A1_name <- "A1_name=TRUE"
}else{
  A1_name <- "A1_name=FALSE"  
}

if(length(grep("A2", Parsing_colnames)) == 1){
  A2_name <- "A2_name=TRUE"
}else{
  A2_name <- "A2_name=FALSE"  
}

if(length(grep("FRQ", Parsing_colnames)) == 1 || length(grep("FRQ", Parsing_colnames)) == 2){
  MAF_calculation_summary <- "MAF_calculation_summary=TRUE"
  number_of_frequency_columns <- paste0("Number_of_frequency_columns=",grep("FRQ", Parsing_colnames))
}else{
  MAF_calculation_summary <- "MAF_calculation_summary=FALSE"  
  number_of_frequency_columns <- "Number_of_frequency_columns=NA"
}

if(length(grep("INFO", Parsing_colnames)) == 1){
  INFO_name <- "INFO_name=TRUE"
}else{
  INFO_name <- "INFO_name=FALSE"  
}

if(length(grep("OR", Parsing_colnames)) == 1){
  OR_name <- "OR_name=TRUE"
}else{
  OR_name <- "OR_name=FALSE"  
}

if(length(grep("BETA", Parsing_colnames)) == 1){
  BETA_name <- "BETA_name=TRUE"
}else{
  BETA_name <- "BETA_name=FALSE"  
}

if(length(grep("SE", Parsing_colnames)) == 1){
  SE_name <- "SE_name=TRUE"
}else{
  SE_name <- "SE_name=FALSE"  
}

if(length(grep("\\bP\\b", Parsing_colnames)) == 1){
  P_name <- "P_value_name=TRUE"
}else{
  P_name <- "P_value_name=FALSE"  
}

fileConn<-file(paste0("new_PRS_set_arguments_for_", Training_set_name, ".txt")
               
writeLines(c(Chromosome_name, SNP_name, BP_name, A1_name, A2_name, MAF_calculation_summary, number_of_frequency_columns,INFO_name, OR_name, BETA_name, SE_name, P_name), fileConn)
close(fileConn)


setkey(PGC_data,CHR)

# separating out chromosome into certain sections
# ideally for future, one of the shell script arguments will be the number of segmentations you wish to separate the chromosome into (maybe even taking from the PBS options)

# below is aka:
  # after setting the key to the chromosome number,
  # take a subset of the data matching the chromsome number
  # and write the subsetted table to a new table for each chromosome. 
for (chromosome.number in 1:22) { 
write.table(eval(parse(text = paste0("PGC_data[J(",chromosome.number,")]"))), file = paste0("PGC_table",chromosome.number,".txt"), quote = F, row.names = F)
}

