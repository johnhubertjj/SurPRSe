### Start Timer
ptm <- proc.time()

########################################
# Adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
getwd()
print(args)

# specify the different input tables #
in_serial <- args[3]
Training_name <- args[4]
Validation_name <- args[5]

# Set up library
library("data.table")

# Set working directory
setwd(".")


if (in_serial == 'TRUE') {
  
# Read in MAF data
MAF_clumped_CLOZUK_data <- fread(paste0("./output/",Validation_name,"_",Training_name,"_FULL_GENOME.frq"))
vector_for_place_of_SNP <- c(1:length(MAF_clumped_CLOZUK_data$CHR))
MAF_clumped_CLOZUK_data[,place.names := vector_for_place_of_SNP]
PGC_useful_columns <- fread(paste0("./output/",Training_name,"_table_for_python.txt"))

# use data table structure to merge tables together and add MAF while removing useless columns
PGC_CLOZUK_python_table <- merge(MAF_clumped_CLOZUK_data,PGC_useful_columns, by = c('CHR','SNP',"A1","A2"), all.x = T, all.y = F)

PGC_CLOZUK_python_table <- PGC_CLOZUK_python_table[, NCHROBS := NULL][,.(CHR,SNP,BP,A1,A2,BETA,P,MAF,place.names)]

# Check alleles
test_1 <- any(nchar(PGC_CLOZUK_python_table$A1) > 1)
test_2 <- any(nchar(PGC_CLOZUK_python_table$A2) > 1)
ones_to_remove <- which(is.na(PGC_CLOZUK_python_table$BP))

# Produce error if it has not worked
if (test_1 == T | test_2 == T) {
  cat("Check these SNPs: ", PGC_CLOZUK_python_table$SNP[ones_to_remove])
  stop("you have some alleles which are larger than 1 bp long, check before clumping")
}

# INSERT CONTROL SO THAT NO EXTRA COLUMNS ARE ADDED (as you've specified that you want to include all x columns)

# Read in .raw file
## File
file <- paste0("./output/",Validation_name,"_",Training_name,"_FULL_GENOME.raw")

## Create connection
con <- file(description=file, open="r")
tmp <- scan(file=con, nlines=1, quiet=TRUE,what = "character")

## Hopefully you know the number of lines from some other source or
#com <- paste("wc -l ", file, " | awk '{ print $1 }'", sep="")
#n <- system(command=com, intern=TRUE)

# RAW.CLOZUK.CLUMPED <- fread("./output/CLOZUK_PGC_FULL_GENOME.raw")
#testa <- PGC_CLOZUK_python_table$place.names[ones_to_remove]
#testb <- PGC_CLOZUK_python_table[ones_to_remove]
#testc <- PGC_CLOZUK_python_table$SNP[ones_to_remove]
#write(testb$SNP,file = "./output/SNPs_which_were_erroneous.txt")
#PGC_CLOZUK_python_table <- PGC_CLOZUK_python_table[-ones_to_remove]

# Flow control for assigning Row ID
setkey(PGC_CLOZUK_python_table, place.names)
testing_col_names <- paste0(PGC_CLOZUK_python_table$SNP,"_",PGC_CLOZUK_python_table$A1)
testing_col_names_raw <- tmp

# Test if the ROW.ID's match that of the converted RAW dataset and stop if they do not
if(all(testing_col_names == testing_col_names_raw[7:length(testing_col_names_raw)]) != T){
  PGC_CLOZUK_python_table[, ROW.ID := c(testing_col_names)]
  write.table(PGC_CLOZUK_python_table, file = paste0("./output/ERROR_ON_ROW_NAMES", Validation_name, "_GWAS_BGE_CLUMPED_", Training_name, "_MAF_FINAL.txt"), quote = F, row.names = T, col.names = F)
  stop("ROW.ID from input data does not match genotype data, check the column names of both .raw file and the merged file")
}else{
  PGC_CLOZUK_python_table[, ROW.ID := c(testing_col_names)]
  PGC_CLOZUK_python_table <- PGC_CLOZUK_python_table[, place.names := NULL]
}

# rearrange order to match Python script downstream
# order is ROW_ID, CHR, SNP, BETA, BP, A1, P, MAF, A2
setcolorder(PGC_CLOZUK_python_table, c(9,1,2,6,3,4,7,8,5)) 

# Write MAF file to working_directory
write.table(PGC_CLOZUK_python_table, file = paste0("./output/", Validation_name, "_GWAS_BGE_CLUMPED_", Training_name, "_MAF_FINAL.txt"), quote = F, row.names = T, col.names = F)
}

# End timer 
proc.time() - ptm