### Start Timer
ptm <- proc.time()

##################################################
# Checking location for serial or batch analysis #
##################################################

System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- 22
  
} else if (whereami == 'JJ') {
  chromosome.number <- 22
  
} else if (whereami == "c1020109") {
  
  # Preparing to run in Job array
  AI <- Sys.getenv("PBS_ARRAY_INDEX")
  chromosome.number <- as.numeric(AI)
  
} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}
# Set up library
library("data.table")

# Set working directory
setwd(".")

# Read in MAF data
MAF_clumped_CLOZUK_data <- fread(paste0("./output/CLOZUK_GWAS_BGE_CLUMPED_chr",chromosome.number,".frq"))
PGC_useful_columns <- fread(paste0("./output/PGC_table_for_python",chromosome.number,".txt"))

# use data table structure to merge tables together and add MAF while removing useless columns
PGC_CLOZUK_python_table <- merge(MAF_clumped_CLOZUK_data,PGC_useful_columns, by = c('CHR','SNP',"A1","A2"), all.x = T, all.y = F)
PGC_CLOZUK_python_table[, NCHROBS := NULL][,.(CHR,SNP,BP,A1,A2,BETA,P,MAF)]
###INSERT CONTROL SO THAT NO EXTRA COLUMNS ARE ADDED (as you've specified that you want to include all x columns)

#Read in .raw file
RAW.CLOZUK.CLUMPED <- fread(paste0("./output/test_clumped_CLOZUK_PGC_PRS_recode_chr",chromosome.number,".raw"))

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
write.table(PGC_CLOZUK_python_table, file = paste0("./output/CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL",chromosome.number,".txt"),quote = F,row.names = T,col.names = F)

# End timer 
proc.time() - ptm