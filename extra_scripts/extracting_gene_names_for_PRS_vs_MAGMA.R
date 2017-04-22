dataset_for_table <- fread("~/Dropbox/PhD_WORK/PhD meetings/PhD_results/0.5_all_outliers_all_SNPs.csv")
Gene_names <- fread("~/Dropbox/PhD_WORK/PhD meetings/PhD_results/gprofiler_results_genes_for_28_03_17_lab_presentation.csv")
setnames(Gene_names, c("number","GENE", "who_cares1", "whocares2", "Gene_name", "whocares3","whocares4"))

Gene_names2 <- Gene_names[,GENE := as.numeric(gsub("WIKIGENE_ACC:","",Gene_names$GENE,perl = T))]
both_matching_genes <- merge(dataset_for_table, Gene_names2, by = "GENE", all = T)

SNP_table_final <- both_matching_genes[,.(GENE,Gene_name,CHR,START,STOP,NSNPS,NPARAM,P.x,P.y)]
setkey(SNP_table_final,P.x)
setnames(SNP_table_final, old = c("GENE","NSNPS","NPARAM", "P.x", "P.y"), c("Entrezgene ID", "NSNPS_MAGMA", "NPARAM_MAGMA", "P_PRS", "P_MAGMA"))

write.csv(SNP_table_final, file = "~/Dropbox/PhD_WORK/PhD meetings/PhD_results/29_03_17lab_presentation_useful_info_table.csv", row.names = F)
