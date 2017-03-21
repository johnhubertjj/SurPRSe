### Start Timer
ptm <- proc.time()

########################################
# Adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
print(args)
getwd()

# specify the different input tables #
Training_name <- args[3]
Validation_name <- args[4]
Sig_thresholds <- as.numeric(args[5])
Validation_dataset <- args[6]
Magma_or_no_Magma_analysis <- args[7]
INFO_decision <- args[8]

# Set up library
library("data.table")

# Set working directory
setwd(".")

# Output a file if no output is specified
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} #else if (length(CLOZUK_data_set_args)==1) {
  # default output file
 # CLOZUK_data_set_args[3] = "./output/out.txt"
#}
CLOZUK_magma_input2 <- fread(Validation_dataset)
names(CLOZUK_magma_input2) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
PGC <- fread(paste0("./output/combined_", Training_name, "_table_with_CHR.POS_identifiers.txt"))
PGC_CLOZUK_merged <- merge(CLOZUK_magma_input2, PGC, by = "SNP", all = F)

if (INFO_decision == "YES") {
  # remove all SNPs with INFO > 0.9
  Info_threshold <- as.numeric(args[9])
  PGC_CLOZUK_merged <- PGC_CLOZUK_merged[INFO > Info_threshold]
}

Sig_thresholds2 <- c("0.0001","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")

x <- NULL
SNPs <- list(x,x,x,x,x,x,x,x,x)
names(SNPs) <- Sig_thresholds

for (i in 1:length(Sig_thresholds)) {
  SNPs_to_extract_for_magma <- which(PGC_CLOZUK_merged$P <= Sig_thresholds[i])
  reference_file <- PGC_CLOZUK_merged[SNPs_to_extract_for_magma]
  reference_file <- reference_file[,.(CHR.x,SNP,P,BETA)]
  names(reference_file) <- c("CHR","SNP","P","BETA")
  SNPs[[i]] <- reference_file
}

if (Magma_or_no_Magma_analysis == "TRUE") {
  
  for (i in 1:length(Sig_thresholds)) {
    write.table(SNPs[[i]], file = paste0("./MAGMA_set_analysis/", Sig_thresholds2[i], Validation_name, "_", Training_name, "_P_Vals_reference_file_magma_analysis.txt"), quote = F, row.names = F)
}
}else{
  for (i in 1:length(Sig_thresholds)) {
    write.table(SNPs[[i]], file = paste0(Sig_thresholds2[i], Validation_name, "_", Training_name, "_P_Vals_reference_file.txt"), quote = F, row.names = F)
  }
}





