# Preparing to run in Job array
AI <- Sys.getenv("PBS_ARRAY_INDEX")
chromosome.number <- as.numeric(AI)

# Set up library
library("data.table")

# Read in MAF data
MAF_clumped_CLOZUK_data <- fread("CLOZUK_GWAS_BGE_CLUMPED_chr22.frq")
PGC_useful_columns <- fread("PGC_table_for_python22.txt")

# use data table structure to merge tables together and add MAF while removing useless columns
PGC_CLOZUK_python_table <- merge(MAF_clumped_CLOZUK_data,PGC_useful_columns, by = c('CHR','SNP',"A1","A2"), all.x = T, all.y = F)
PGC_CLOZUK_python_table[, NCHROBS := NULL][,.(CHR,SNP,BP,A1,A2,BETA,P,MAF)]
###INSERT CONTROL SO THAT NO EXTRA COLUMNS ARE ADDED (as you've specified that you want to include all x columns)

#Read in .raw file
RAW.CLOZUK.CLUMPED <- fread("test_clumped_CLOZUK_PGC_PRS_recode.raw")


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
write.table(PGC_CLOZUK_python_table, file = "CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL.txt",quote = F,row.names = T,col.names = F)

