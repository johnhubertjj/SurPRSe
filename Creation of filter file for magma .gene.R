library(data.table)

CLOZUK_magma_input2 <- fread("CLOZUK_GWAS_BGE_chr22_magma_input_2.bim")
names(CLOZUK_magma_input2) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
PGC <- fread ("PGCtable_22")
PGC_CLOZUK_merged <- merge(CLOZUK_magma_input2,PGC,by = "SNP", all = F)

sig <- c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)

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
  write.table(SNPs[[i]], file = paste0(sig[i],"CLOZUK_PGC_P_Vals_reference_file.txt"),quote = F,row.names = F)
}

