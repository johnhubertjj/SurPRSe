MAGMA_different_gene_regions <-fread("NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
colnames(MAGMA_different_gene_regions) <- c("ID","Chromosome","BP_start","BP_finish","strand","Gene_symbol")
MAGMA_different_gene_regions$BP_start <- MAGMA_different_gene_regions$BP_start - 35000
MAGMA_different_gene_regions$BP_finish <- MAGMA_different_gene_regions$BP_finish + 10000

write.table(MAGMA_different_gene_regions, file = "NCBI37_-35kb_+10kb.gene.loc",quote = F, col.names = F, row.names = F)
