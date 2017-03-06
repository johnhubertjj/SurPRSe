# Obtaining SNPs within Genes common to PRS and MAGMA
sig <- c("0.0001","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")

# loop through the significance thresholds and get the SNPs required to run the scripts in MAGMA
for (i in 1:length(sig)) {
  common_genes <- fread(paste0(sig[i],"merged_PRS_MAGMA_COMMON_GENES.txt"))
  common_genes <- unique(common_genes)
  names(common_genes) <- "GENE"
  Annotation_file <- fread("./output/MAGMA_Gene_regions_for_python_script_chr_22.txt")
  a <- merge(common_genes,Annotation_file, by = 'GENE')
  SNP_to_use <- a$SNP
  write(SNP_to_use, paste0(sig[i], "update_Bed_files_for_MAGMA_analysis.txt"))
}
