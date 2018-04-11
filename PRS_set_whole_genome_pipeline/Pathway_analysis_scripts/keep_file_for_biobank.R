library(data.table)
biobank_phenotype <- fread("~/Documents/Biobank_Alspac_work/Image_UK_Irish_with8cov.txt")
Keep_file <- biobank_phenotype[,.(app6553,app6553)]
write.table(Keep_file,file = "~/Documents/Biobank_Alspac_work/biobank_wv_1_3000_brains.keep", row.names = F, col.names = F, quote = T)
