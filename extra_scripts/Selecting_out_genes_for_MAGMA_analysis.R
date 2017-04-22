NCBI_gene_regions <- fread("/Volumes/HD-PCU2/Stationary_data/NCBI37.3.gene.loc")
PGC_data <- fread("/Volumes/HD-PCU2/PGC_noCLOZUK_data/PGC_table14_new.txt")
NCBI_gene_regions <- NCBI_gene_regions[V2 == 14]
NCBI_gene_regions$CHR <- as.numeric(NCBI_gene_regions$CHR)
setnames(NCBI_gene_regions,c("GENE_ID", "CHR", "BP_START","BP_END","STRAND","GENE_NAME"))
integers <- which(PGC_data$BP > NCBI_gene_regions$BP_START[100] & PGC_data$BP < NCBI_gene_regions$BP_END[100])

New_gene_only_data <- PGC_data[integers]
New_gene_only_data$CHR <- as.character(New_gene_only_data$CHR)
setkey(New_gene_only_data,CHR)
New_gene_only_data <- New_gene_only_data[,.(SNP,CHR,BP)]

write.table(New_gene_only_data, "/Volumes/HD-PCU2/PGC_noCLOZUK_data/testing_pathway_and_gene_MAGMA_chr_14_gene_85446.txt",quote = F,row.names = F)
