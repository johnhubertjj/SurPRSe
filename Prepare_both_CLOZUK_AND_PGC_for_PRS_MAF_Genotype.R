### Start Timer
ptm <- proc.time()

##################################################
# Checking location for serial or batch analysis #
##################################################
args <- commandArgs(trailingOnly = T)
System_information <- Sys.info()
whereami <- System_information['user']

# Set up library
library("data.table")

# Set working directory
setwd(".")

###############################
## Location Finding Function ##
###############################

Batch.assignment <- function(wherami) {
  if (whereami == "johnhubert") {
    assign("chromosome.number", 22, envir = .GlobalEnv)
    
  } else if (whereami == 'JJ') {
    assign("chromosome.number", 22, envir = .GlobalEnv)
    
  } else if (whereami == "c1020109") {
    
    # Preparing to run in Job array
    AI <- Sys.getenv("PBS_ARRAY_INDEX")
    AI <- as.numeric(AI)
    assign("chromosome.number", AI, envir = .GlobalEnv)
    
  } else {
    stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
  }
}

# Find out whether you want to run the script in a batch format or in serial
if( args[1] == 'Batch') {
  
  # Assign chromosome number as job array variable
  Batch.assignment(whereami)
  
  # Read in MAF data
  MAF_clumped_CLOZUK_data <- fread(paste0("./output/CLOZUK_GWAS_BGE_CLUMPED_chr",chromosome.number,".frq"))
  PGC_useful_columns <- fread(paste0("./output/PGC_table_for_python",chromosome.number,".txt"))
  
  # use data table structure to merge tables together and add MAF while removing useless columns
  PGC_CLOZUK_python_table <- merge(MAF_clumped_CLOZUK_data,PGC_useful_columns, by = c('CHR','SNP',"A1","A2"), all.x = T, all.y = F)
  PGC_CLOZUK_python_table[, NCHROBS := NULL][,.(CHR,SNP,BP,A1,A2,BETA,P,MAF)]
  ###INSERT CONTROL SO THAT NO EXTRA COLUMNS ARE ADDED (as you've specified that you want to include all x columns)
  
  #Read in .raw file
  RAW.CLOZUK.CLUMPED <- fread(paste0("./output/CLOZUK_GWAS_BGE_CLUMPED_chr",chromosome.number,".raw"))
  
  # Flow control for assigning Row ID
  testing_col_names <- paste0(PGC_CLOZUK_python_table$SNP,"_",PGC_CLOZUK_python_table$A1)
  testing_col_names_raw <- colnames(RAW.CLOZUK.CLUMPED)
  
  if(all(testing_col_names == testing_col_names_raw[7:length(testing_col_names_raw)]) != T){
    stop("ROW.ID from input data does not match genotype data, check the column names of both .raw file and the merged file")
  }else{
    PGC_CLOZUK_python_table[, ROW.ID := c(testing_col_names)]
  }
  # rearrange order to match Python script downstream
  setcolorder(PGC_CLOZUK_python_table, c(9,1,2,6,3,4,7,8,5))
  
  # Write MAF file to working_directory
  write.table(PGC_CLOZUK_python_table, file = paste0("./output/CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL",chromosome.number,".txt"),quote = F, row.names = T, col.names = F)
  
} else if (args[1] == 'Serial') {
  
# Read in MAF data
MAF_clumped_CLOZUK_data <- fread("./output/CLOZUK_PGC_FULL_GENOME.frq")
vector_for_place_of_SNP <- c(1:length(MAF_clumped_CLOZUK_data$CHR))
MAF_clumped_CLOZUK_data[,place.names := vector_for_place_of_SNP]
PGC_useful_columns <- fread("./output/PGC_table_for_python.txt")

# use data table structure to merge tables together and add MAF while removing useless columns
PGC_CLOZUK_python_table <- merge(MAF_clumped_CLOZUK_data,PGC_useful_columns, by = c('CHR','SNP',"A1","A2"), all.x = T, all.y = F)
PGC_CLOZUK_python_table[, NCHROBS := NULL][,.(CHR,SNP,BP,A1,A2,BETA,P,MAF,place.names)]

# Check alleles
test_1 <- any(nchar(PGC_CLOZUK_python_table$A1) > 1)
test_2 <- any(nchar(PGC_CLOZUK_python_table$A2) > 1)

if (test_1 == T | test_2 == T) {
  warning("you have some alleles which are larger than 1 bp long, check before clumping")
}

# INSERT CONTROL SO THAT NO EXTRA COLUMNS ARE ADDED (as you've specified that you want to include all x columns)

# Read in .raw file
## File
file <- "./output/CLOZUK_PGC_FULL_GENOME.raw"

## Create connection
con <- file(description=file, open="r")
tmp <- scan(file=con, nlines=1, quiet=TRUE,what = "character")

## Hopefully you know the number of lines from some other source or
#com <- paste("wc -l ", file, " | awk '{ print $1 }'", sep="")
#n <- system(command=com, intern=TRUE)

# RAW.CLOZUK.CLUMPED <- fread("./output/CLOZUK_PGC_FULL_GENOME.raw")
ones_to_remove <- which(is.na(PGC_CLOZUK_python_table$BP))
#testa <- PGC_CLOZUK_python_table$place.names[ones_to_remove]
#testb <- PGC_CLOZUK_python_table[ones_to_remove]
#testc <- PGC_CLOZUK_python_table$SNP[ones_to_remove]
#write(testb$SNP,file = "./output/SNPs_which_were_erroneous.txt")
#PGC_CLOZUK_python_table <- PGC_CLOZUK_python_table[-ones_to_remove]

# Flow control for assigning Row ID
setkey(PGC_CLOZUK_python_table, place.names)
testing_col_names <- paste0(PGC_CLOZUK_python_table$SNP[place_names_for_raw],"_",PGC_CLOZUK_python_table$A1[place_names_for_raw])
testing_col_names_raw <- tmp

if(all(testing_col_names == testing_col_names_raw[7:length(testing_col_names_raw)]) != T){
  stop("ROW.ID from input data does not match genotype data, check the column names of both .raw file and the merged file")
}else{
  PGC_CLOZUK_python_table[, ROW.ID := c(testing_col_names)]
  PGC_CLOZUK_python_table <- PGC_CLOZUK_python_table[, place.names := NULL]
}

# rearrange order to match Python script downstream
setcolorder(PGC_CLOZUK_python_table, c(9,1,2,6,3,4,7,8,5)) 

# Write MAF file to working_directory
write.table(PGC_CLOZUK_python_table, file = paste0("./output/CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL.txt"), quote = F, row.names = T, col.names = F)

}

# End timer 
proc.time() - ptm