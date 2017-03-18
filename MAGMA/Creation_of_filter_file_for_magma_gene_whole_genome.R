
# Read in arguments for script
CLOZUK_data_set_args <- commandArgs(trailingOnly = T)
print(CLOZUK_data_set_args)

Sig <- as.numeric(CLOZUK_data_set_args[1])
  
library(data.table)

# Output a file if no output is specified
if (length(CLOZUK_data_set_args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(CLOZUK_data_set_args)==1) {
  # default output file
  CLOZUK_data_set_args[3] = "./output/out.txt"
}

# CLOZUK_magma_input2 <- fread(CLOZUK_data_set_args[1])

CLOZUK_magma_input2 <- fread("./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.bim")
names(CLOZUK_magma_input2) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
PGC <- fread ("./output/PGC_table_for_python.txt")
PGC_CLOZUK_merged <- merge(CLOZUK_magma_input2,PGC,by = "SNP", all = F)

sig <- c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)
sig2 <- c("0.0001","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")
x <- NULL
SNPs <- list(x,x,x,x,x,x,x,x,x)
names(SNPs) <- sig

for (i in 1:length(sig)) {
  SNPs_to_extract_for_magma <- which(PGC_CLOZUK_merged$P <= sig[i])
  reference_file <- PGC_CLOZUK_merged[SNPs_to_extract_for_magma]
  reference_file <- reference_file[,.(CHR.x,SNP,P,BETA)]
  names(reference_file) <- c("CHR","SNP","P","BETA")
  SNPs[[i]] <- reference_file
}

for (i in 1:length(sig)) {
  write.table(SNPs[[i]], file = paste0("./output/",sig2[i],"CLOZUK_PGC_P_Vals_reference_file.txt"),quote = F,row.names = F)
}





