## This script would inherit another script (probably python this time) that decides whether or not the summary stats /
## data has enough information to run an analysis the way you wish /
## and sets arguments accordingly.

### Start Timer
ptm <- proc.time()

### Library
library(data.table)

## set wd
setwd(".")

##################################################
# Checking location for serial or batch analysis #
##################################################

System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- args[11]
  
} else if (whereami == "c1020109") {
  
  # Preparing to run in Job array
  AI <- Sys.getenv("PBS_ARRAY_INDEX")
  chromosome.number <- as.numeric(AI)
  
} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}

########################################
# adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
getwd()
print(args)

# specify the different input tables #
Training_datatable <- paste0(args[3],".txt")
Training_datatable_output <- paste0("./output/", args[3], "_new.txt")
Training_name <- args[4]
MAF_decision <- args[5]
MAF_threshold <- args[6]
INFO_decision <- args[7]
Info_threshold <- as.numeric(args[8])
SE_decision <- args[9]
SE_threshold <- as.numeric(args[10])

# Read in Training set
PGC_table <- fread(Training_datatable)
changed_PGC_table <- copy(PGC_table)

if (SE_decision == "TRUE"){
  changed_PGC_table <- changed_PGC_table[SE > SE_threshold]
}

if (INFO_decision == "TRUE") {
  # remove all SNPs with INFO > 0.9
  Info_threshold <- as.numeric(args[7])
  changed_PGC_table <- changed_PGC_table[INFO > Info_threshold]
}

# Do I want to remove SNPs based on MAF
if (MAF_decision == "TRUE") {
 
  # This is assuming that the FRQ_A columns exist 
  Allele_cases_location <- grep("FRQ_A_\\d+",colnames(changed_PGC_table), perl = T)
  Allele_controls_location <-grep("FRQ_U_\\d+",colnames(changed_PGC_table), perl = T)

  # Select out the case and control columns without prior knowledge of which is which and rename to FRQ_A and FRQ_U
  # Keep old column names for reassignment
  old_column_names <- colnames(changed_PGC_table[,c(Allele_cases_location, Allele_controls_location), with = F])
  new_column_names_cases <- gsub(pattern = "(?<=FRQ_A)_\\d+", replacement = "", x = names(changed_PGC_table[,Allele_cases_location, with = F]), perl = T)
  new_column_names_controls <- gsub(pattern = "(?<=FRQ_U)_\\d+", replacement = "", x = colnames(changed_PGC_table[,Allele_controls_location,with = F]), perl = T)

  # set new names for columns
  setnames(changed_PGC_table, c(Allele_cases_location,Allele_controls_location), c(new_column_names_cases, new_column_names_controls))

  # Change all values to MAF,
  # Frequencies are just stated so in order to find the MAF you need to find the values for which the frequency is above 0.5 for 
  # each cases and controls individually, and subtract one from each in order to find MAF
  Vector_of_major_alleles_FRQ_A <- which(changed_PGC_table$FRQ_A > 0.5)
  Vector_of_major_alleles_FRQ_U <- which(changed_PGC_table$FRQ_U > 0.5)
  changed_PGC_table <- changed_PGC_table[Vector_of_major_alleles_FRQ_A, FRQ_A_MAF := (1 - FRQ_A)]
  changed_PGC_table <- changed_PGC_table[Vector_of_major_alleles_FRQ_U, FRQ_U_MAF := (1 - FRQ_U)]
  changed_PGC_table <- changed_PGC_table[is.na(FRQ_A_MAF), FRQ_A_MAF := FRQ_A]
  changed_PGC_table <- changed_PGC_table[is.na(FRQ_U_MAF), FRQ_U_MAF := FRQ_U]
  
  # remove all SNPs lower than 0.01
  FRQ_to_keep_helpful3 <- which(changed_PGC_table$FRQ_A_MAF <= 0.01 & changed_PGC_table$FRQ_U_MAF <= 0.01)
  changed_PGC_table <- changed_PGC_table[!FRQ_to_keep_helpful3]
  
  #Find out which columns have one value below MAF 0.01
  FRQ_U_helpful2 <- which(changed_PGC_table$FRQ_A_MAF <= 0.01 & changed_PGC_table$FRQ_U_MAF > 0.01)
  FRQ_A_helpful1 <- which(changed_PGC_table$FRQ_A_MAF > 0.01 & changed_PGC_table$FRQ_U_MAF <= 0.01)

  # add extra column with MAF for each SNP according to the rule, keep the minimum allele frequency (except those below 0.01)
  changed_PGC_table[,MAF := pmin(FRQ_A_MAF,FRQ_U_MAF)]
  changed_PGC_table[FRQ_A_helpful1, MAF := changed_PGC_table[FRQ_A_helpful1,FRQ_A_MAF]]
  changed_PGC_table[FRQ_U_helpful2, MAF := changed_PGC_table[FRQ_U_helpful2,FRQ_U_MAF]]
  
  # Do a check 
  if (all(changed_PGC_table$MAF == changed_PGC_table$FRQ_A_MAF | changed_PGC_table$FRQ_U_MAF) == FALSE){
    warning("MAF is not neccessarily correct for this chromosome, check back please")
  }

  # Remove un-required columns for previous MAF thresholding
  changed_PGC_table[, c("FRQ_A_MAF", "FRQ_U_MAF") := NULL ]
  
  # Reassign the column names back to the table
  setnames(changed_PGC_table, c(Allele_cases_location, Allele_controls_location), old_column_names)
}

  
# write new table overwriting the previous PGC_new_table
write.table(changed_PGC_table, file = Training_datatable_output, row.names = F, quote = F)
cat("Number of SNPs in", Training_name, "before MAF and/or INFO correction: N=", nrow(PGC_table))
cat("Number of SNPs in", Training_name, "after MAF and/or INFO correction: N=", nrow(changed_PGC_table))

#End Timer
proc.time() - ptm
