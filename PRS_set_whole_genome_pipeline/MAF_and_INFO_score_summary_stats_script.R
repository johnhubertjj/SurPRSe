## This script would inherit another script (probably python this time) that decides whether or not the summary stats /
## data has enough information to run an analysis the way you wish /
## and sets arguments accordingly.

### Start Timer
ptm <- proc.time()

### Library
library(data.table)

## set wd
setwd(".")

########################################
# adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
getwd()
print(args)

# specify the different input tables #
Training_datatable <- paste0(args[3],".txt")
Training_name <- args[4]
Validation_name <- args [5]
Training_datatable_output <- paste0("./", Training_name, "_", Validation_name,"_output/", args[3], "_new.txt")
MAF_decision <- args[6]
MAF_threshold <- as.numeric(args[7])
INFO_decision <- args[8]
Info_threshold <- as.numeric(args[9])
SE_decision <- args[10]
SE_threshold <- as.numeric(args[11])
chromosome.number <- args[12]

# Read in Training set
PGC_table <- fread(Training_datatable)
changed_PGC_table <- copy(PGC_table)

cat("Number of SNPS in ", Training_name,"before MAF/INFO/SE script Chr:",chromosome.number, "N=" ,nrow(PGC_table))

### Remove Duplicated SNPs in Training here ###
PGC.duplicated.removed <- which(duplicated(PGC_table$SNP))
PGC.duplicated.removed.rev <- which(duplicated(PGC_table$SNP, fromLast = T))

### One-line duplication check - Training ###
length(PGC_table$SNP)
if (length(PGC.duplicated.removed) >= 1){
  PGC_duplicate_SNPS <- PGC_table$SNP[(duplicated(PGC_table$SNP) | duplicated(PGC_table$SNP, fromLast = TRUE)) ]
  PGC_table <- PGC_table[!(duplicated(PGC_table$SNP) | duplicated(PGC_table$SNP, fromLast = TRUE)) ]
}
length(PGC_table$SNP)

cat("Number of SNPS in ",Training_name," after duplication check one Chr:",chromosome.number, "N=" ,nrow(PGC_table))

if (SE_decision == "TRUE"){
  changed_PGC_table <- changed_PGC_table[SE < SE_threshold]
  cat("SNPs with SE > ", SE_threshold, " removed")
}

if (INFO_decision == "TRUE") {
  # remove all SNPs with INFO > 0.9
  changed_PGC_table <- changed_PGC_table[INFO > Info_threshold]
  cat("SNPs with INFO > ", Info_threshold, " removed")
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
  FRQ_to_keep_helpful3 <- which(changed_PGC_table$FRQ_A_MAF <= MAF_threshold & changed_PGC_table$FRQ_U_MAF <= MAF_threshold)
  changed_PGC_table <- changed_PGC_table[!FRQ_to_keep_helpful3]
  
  #Find out which columns have one value below MAF 0.01
  FRQ_U_helpful2 <- which(changed_PGC_table$FRQ_A_MAF <= MAF_threshold & changed_PGC_table$FRQ_U_MAF > MAF_threshold)
  FRQ_A_helpful1 <- which(changed_PGC_table$FRQ_A_MAF > MAF_threshold & changed_PGC_table$FRQ_U_MAF <= MAF_threshold)

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
  cat("SNPs with MAF > ", MAF_threshold, " removed")
}

  
# write new table overwriting the previous PGC_new_table
write.table(changed_PGC_table, file = Training_datatable_output, row.names = F, quote = F)
cat("Number of SNPs in", Training_name, "before MAF and/or INFO and/or SE correction: N=", nrow(PGC_table))
cat("Number of SNPs in", Training_name, "after MAF and/or INFO and/or SE correction: N=", nrow(changed_PGC_table))

#End Timer
proc.time() - ptm
