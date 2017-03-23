## removing Genes with one SNP involved

sig <- c(0.0001, 0.001, 0.01 ,0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

for (i in 1:length(sig)){
  current_table <- fread(paste0("~/Dropbox/PhDwork-PRS/",sig[i], "CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file_whole_genome.txt"))
  current_table <- current_table[is.finite(rowSums(current_table))]
  setnames(current_table, c("pval", "Genes"))
  One_SNP_table <- fread(paste0("~/Dropbox/PhDwork-PRS/whole_Genome_test_",sig[i],"_one_SNP_only_dictionary.txt"))
  setnames(One_SNP_table, c("Genes", "pval_thresh"))
  testat <- merge(current_table,One_SNP_table, by = "Genes", all = T)
  testat <- testat[is.na(testat$pval_thresh)]
  testat[,pval_thresh := NULL]
  write.table(testat,file = paste0("~/Dropbox/PhDwork-PRS/",sig[i], "CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file_whole_genome_more_than_one_SNP.txt"), quote = F, row.names = F, col.names = F)
}


