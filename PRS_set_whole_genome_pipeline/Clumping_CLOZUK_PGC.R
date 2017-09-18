### Start Timer
ptm <- proc.time()

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)


###############################
# Checking location for serial or batch analysis #
###############################

System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- args[7]
  
} else if (whereami == 'JJ') {
  chromosome.number <- args[7]
  
} else if (whereami == "c1020109") {
  
  # Preparing to run in Job array
  AI <- Sys.getenv("PBS_ARRAY_INDEX")
  chromosome.number <- as.numeric(AI)
  
} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}

# specify the different input tables #
Training_datatable <- paste0("./output/", args[3],"_new.txt")
Validation_datatable_bim_file <- paste0("./output/", args[4], "_2.bim")
Common_SNPs <- paste0("./output/chr", chromosome.number, args[5], "_", args[6], "_common_SNPs.txt")

### Library ###

library(data.table)

### Read in updated_raw data
setwd(".")
getwd()
PGC <- fread(Training_datatable)
CLOZUK.data <- fread(Validation_datatable_bim_file)
Common_SNPs <- scan(Common_SNPs, what = "numeric")

### Check for duplicated SNPs.PGC
PGCa <- which(duplicated(PGC$SNP))
PGCb <- which(duplicated(PGC$SNP, fromLast = T))
PGC_duplicated_SNPs <- c(PGCa,PGCb)
removed_SNPs <- PGC[PGC_duplicated_SNPs]
removed_SNPs <- removed_SNPs$SNP

### Check for duplicated SNPs.CLOZUK
### Add rownames so it is more readable

cloza <- which(duplicated(CLOZUK.data$V2))
clozb <- which(duplicated(CLOZUK.data$V2,fromLast = T))
removed_SNPs_cloz <- c(cloza,clozb)
removed_SNPs_cloz <- CLOZUK.data[removed_SNPs_cloz]
removed_SNPs_cloz <- removed_SNPs_cloz$V2

### Extract duplicate SNPs
extract_SNPs <-c(removed_SNPs,removed_SNPs_cloz)
extract_SNPs <- unique(extract_SNPs)

### write SNPs back to file
write(extract_SNPs,file = paste0("./output/extracted_Duplicate_snps_chr", chromosome.number,".txt"))
which(removed_SNPs %in% Common_SNPs)

#End Timer
proc.time() - ptm
