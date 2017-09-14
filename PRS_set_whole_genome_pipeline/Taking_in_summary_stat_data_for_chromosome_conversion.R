### Start Timer
ptm <- proc.time()

## Load packages ## 
library(data.table)

## setting working directory
setwd(".")

## new_env
e <- new.env()

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

Training_name_full_unseparated <- args[3]
Training_set_name <- args[4]
Validation_set_name <- args[5]
MAF_summary <- args[6]
INFO_summary <- args[7]
SE_summary <- args[8]
SE_threshold <- as.numeric(args[9])

Chromosomes_to_split <- as.numeric(args[c(10:length(args))])
print(Chromosomes_to_split)

# Chromosomes_to_split <- args[9]


# Training_name_full_unseparated <- "/Users/johnhubert/Documents/testing_cross_disorder/BPSCZ.bp_v_scz.results.txt"
# Training_set_name <- "BPvsSCZ"

# Function to calculate what type of Training dataset we have
What_summary_columns_do_we_have <- function(Parsing_colnames){

## what type of files do we have?
if(length(grep("\\bCHR\\b", Parsing_colnames)) == 1 | length(grep("\\bHG19CHR\\b", Parsing_colnames)) == 1 | length(grep("\\bCHROMOSOME\\b", Parsing_colnames)) == 1 ){
  Chromosome_name <- "CHR_name=TRUE"
  if(length(grep("\\bHG19CHR\\b", Parsing_colnames)) == 1){
      place.of.change.snp <- grep("\\bHG19CHR\\b", Parsing_colnames)
      Parsing_colnames[place.of.change.snp] <- "CHR"
  }
  if(length(grep("\\bCHROMOSOME\\b", Parsing_colnames)) == 1){
    place.of.change.snp <- grep("\\bCHROMOSOME\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.snp] <- "CHR"
  }
  
}else{
  Chromosome_name <- "CHR_name=FALSE"  
}

if(length(grep("\\bSNP\\b", Parsing_colnames)) == 1 | length(grep("\\bMARKERNAME\\b", Parsing_colnames)) == 1 | length(grep("\\bSNPID\\b", Parsing_colnames)) == 1 | length(grep("\\bRSID\\b", Parsing_colnames)) == 1 ) {
  SNP_name <- "SNP_name=TRUE"
  if(length(grep("\\bMARKERNAME\\b", Parsing_colnames)) == 1){
      place.of.change.snp <- grep("\\bMARKERNAME\\b", Parsing_colnames)
      Parsing_colnames[place.of.change.snp] <- "SNP"
  }
  if(length(grep("\\bSNPID\\b", Parsing_colnames)) == 1){
      place.of.change.snp <- grep("\\bSNPID\\b", Parsing_colnames)
      Parsing_colnames[place.of.change.snp] <- "SNP"
  }
  if(length(grep("\\bRSID\\b", Parsing_colnames)) == 1){
    place.of.change.snp <- grep("\\bRSID\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.snp] <- "SNP"
  }

}else{
  SNP_name <- "SNP_name=FALSE"  
  stop("Training dataset does NOT have a 'SNP' column essential for the pipeline to work, please add a 'SNP' column")
  
}

if(length(grep("\\bBP\\b", Parsing_colnames)) == 1 | length(grep("\\bPOS\\b", Parsing_colnames)) == 1 | length(grep("\\bPOSITION\\b", Parsing_colnames)) == 1){
  BP_name <- "BP_name=TRUE"
  if(length(grep("\\bPOS\\b", Parsing_colnames)) == 1){
    place.of.change.bp <- grep("\\bPOS\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.bp] <- "BP"
  }
  
  if(length(grep("\\bPOSITION\\b", Parsing_colnames)) == 1){
    place.of.change.bp <- grep("\\bPOSITION\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.bp] <- "BP"
  }
  
}else{
  BP_name <- "BP_name=FALSE"  
  stop("Training dataset does NOT have a 'BP' column essential for the pipeline to work, please add a 'BP' column")
}
  
if(length(grep("\\bFRQ", Parsing_colnames)) == 1 || length(grep("\\bFRQ", Parsing_colnames)) == 2){
  MAF_calculation_summary <- paste0("MAF_Summary=", MAF_summary)
  number_of_frequency_columns <- paste0("Number_of_frequency_columns=", grep("FRQ", Parsing_colnames))
}else{
  MAF_calculation_summary <- "MAF_summary=FALSE"  
  number_of_frequency_columns <- "Number_of_frequency_columns=NA"
  warning("MAF frequency calculation is impossible without an allele frequency estimation, please calculate MAF using Validation dataset or other genotyped data")
  if(MAF_summary == "TRUE"){
    warning("MAF_summary input argument has been changed from TRUE to FALSE due to constraints warned above")
  }
}
  
if(length(grep("\\bA1\\b", Parsing_colnames)) == 1 | length(grep("\\bREF\\b", Parsing_colnames)) == 1 ){
  A1_name <- "A1_name=TRUE"
  if(length(grep("\\bREF\\b", Parsing_colnames)) == 1){
    place.of.change.bp <- grep("\\bREF\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.bp] <- "A1"
    warning("REF allele heading found, converting to A1 heading and inferring that it is the reference allele to which the OR and/or BETA values are based on. IF NOT reference allele, please change file to appropriate heading, for safety, MAF will not be calculated using training set")
    MAF_calculation_summary <- "MAF_summary=FALSE"  
  }
}else{
  A1_name <- "A1_name=FALSE"
  stop("Alternative Allele is not present in the Training dataset or is not named \"A1\", please change column headers or add a column with Alternative allele")
}

if(length(grep("\\bA2\\b", Parsing_colnames)) == 1| length(grep("\\bALT\\b", Parsing_colnames)) == 1 ){
  A2_name <- "A2_name=TRUE"
  if(length(grep("\\bALT\\b", Parsing_colnames)) == 1){
    place.of.change.bp <- grep("\\bALT\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.bp] <- "A2"
    warning("ALT allele heading found, converting to A2 heading and inferring that it is the alternative allele which WILL NOT match the BETA/OR column for this allele. IF NOT alternative allele, please change file to appropriate heading, for safety, MAF will not be calculated using this training set")
    MAF_calculation_summary <- "MAF_summary=FALSE"  
  }
}else{
  A2_name <- "A2_name=FALSE"  
  warning("Reference allele is not present in the Training dataset or is not named \"A2\", please change column headers or add a column with Reference allele")
}


if(length(grep("\\bINFO\\b", Parsing_colnames)) == 1){
  INFO_name <- "INFO_summary=TRUE"
  if (INFO_summary == "FALSE"){
    warning("INFO score has been supplied but will NOT be used for post-GWAS-QC, if human error, please change the argument to: INFO_Summary=TRUE")
    INFO_name <- "INFO_summary=FALSE"
  }
}else{
  INFO_name <- "INFO_summary=FALSE" 
  warning("INFO score not avaliable in summary stats dataset, will attempt to use Standard error instead (removing SNPS with SE >= ", SE_threshold, " )")
  SE_summary <- "TRUE"
  if (INFO_summary == "TRUE"){
    warning("INFO_summary input argument has been changed from TRUE to FALSE due to constraints warned above")
  }
}

if(length(grep("\\bOR\\b", Parsing_colnames)) == 1){
  OR_name <- "OR_name=TRUE"
}else{
  OR_name <- "OR_name=FALSE"  
}

if(length(grep("\\bBETA\\b", Parsing_colnames)) == 1){
  BETA_name <- "BETA_name=TRUE"
}else{
  BETA_name <- "BETA_name=FALSE"  
}

if(BETA_name == "BETA_name=FALSE" & OR_name == "OR_name=FALSE"){
  stop("No Odds ratios or BETA scores supplied in Training set, please calculate at least one before attempting PRSset")
}
  
if(length(grep("\\bSE\\b", Parsing_colnames)) == 1 & INFO_name == "INFO_summary=FALSE"){
  SE_name <- "SE_summary=TRUE"
  warning("SE column found, will use SE as a Post-GWAS-QC measure")
  
} else if (length(grep("\\bSE\\b", Parsing_colnames)) == 1 & INFO_name == "INFO_summary=TRUE") {
  SE_name <- "SE_summary=FALSE"
  
  if(SE_summary == "TRUE"){
    warning("SE_summary input argument has been changed from TRUE to FALSE because no Standard Error column is avaliable in the dataset")
  }
  if (INFO_name == "INFO_summary=FALSE"){
    warning("No INFO score or SE estimation for the Training set has been provided, No post-QC with INFO score or Standard error will be performed")
  }
}

if(length(grep("\\bP\\b", Parsing_colnames)) == 1 | length(grep("\\bPVAL\\b", Parsing_colnames)) == 1 | length(grep("\\bP_VALUE\\b", Parsing_colnames)) == 1){
  P_name <- "P_value_name=TRUE"
  if(length(grep("\\bPVAL\\b", Parsing_colnames)) == 1){
    place.of.change.P <- grep("\\bPVAL\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.P] <- "P"
  }
  if(length(grep("\\bP_VALUE\\b", Parsing_colnames)) == 1){
    place.of.change.P <- grep("\\bP_VALUE\\b", Parsing_colnames)
    Parsing_colnames[place.of.change.P] <- "P"
  }
}else{
  P_name <- "P_value_name=FALSE"
  stop("No column matching \"P\" are supplied in this dataset, please re-name headers or calculate P_values")
}
assign("Parsing_colnames",Parsing_colnames, envir = e)
  
fileConn<-file(paste0("./",Training_set_name,"_",Validation_set_name,"_extrainfo/new_PRS_set_arguments_for_", Training_set_name, ".txt"))
               
writeLines(c(Chromosome_name, SNP_name, BP_name, A1_name, A2_name, MAF_calculation_summary, number_of_frequency_columns, INFO_name, OR_name, BETA_name, SE_name, P_name), fileConn)
close(fileConn)
}

# Get wd for checking
getwd()

## Reading in PGC data
## Select for CHR 22
Training_data <- fread(Training_name_full_unseparated)
Parsing_colnames <- colnames(Training_data)
Parsing_colnames <- toupper(Parsing_colnames)

## Run function to determine what type of data we have and change to standardised Column names
What_summary_columns_do_we_have(Parsing_colnames)
colnames(Training_data) <- e$Parsing_colnames


if(length(grep("CHR", Parsing_colnames)) == 1){
setkey(Training_data, CHR)
}else{
  stop("No chromosome column found in summary stats dataset, re-organise dataset so a chromosome number identifier is a column")
}

# separating out chromosome into certain sections

# below is aka:
  # after setting the key to the chromosome number,
  # take a subset of the data matching the chromsome number
  # and write the subsetted table to a new table for each chromosome. 
for (chromosome.number in Chromosomes_to_split) { 
write.table(eval(parse(text = paste0("Training_data[J(",chromosome.number,")]"))), file = paste0(Training_set_name,"_table", chromosome.number,".txt"), quote = F, row.names = F)
}

