# Preparing to run in Job array
AI <- Sys.getenv("PBS_ARRAY_INDEX")
chromosome.number <- as.numeric(AI)

### Library

library(data.table)

### Read in updated_raw data
setwd(".")
PGC <- fread(paste0("PGC_table",chromosome.number,".txt"))
CLOZUK.data <- fread(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))
Common_SNPs <- scan(paste0("chr",chromosome.number,"PGC_CLOZUK_common_SNPs.txt"),what = "numeric")

### Check for duplicated SNPs.PGC
PGCa <- which(duplicated(PGC$SNP))
PGCb <- which(duplicated(PGC$SNP,fromLast = T))
PGC_duplicated_SNPs <- c(a,b)
removed_SNPs <- PGC[PGC_duplicated_SNPs]
removed_SNPs <- removed_SNPs$SNP

### Check for duplicated SNPs.CLOZUK
### Add rownames so it is more readable

cloza <- which(duplicated(CLOZUK$V2))
clozb <- which(duplicated(PGC$SNP,fromLast = T))
removed_SNPs_cloz <- c(cloza,clozb)
removed_SNPs_cloz <- CLOZUK[removed_SNPs_cloz]
removed_SNPs_cloz <- removed_SNPs_cloz$V2

### Extract duplicate SNPs
extract_SNPs <-c(removed_SNPs,removed_SNPs_cloz)

### write SNPs back to file
write(extract_SNPs,file = "extracted_Duplicate_snps_chr",chromosome.number,".txt")
which(removed_SNPs %in% Common_SNPs)
