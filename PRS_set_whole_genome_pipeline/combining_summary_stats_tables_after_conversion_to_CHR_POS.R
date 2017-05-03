## Create a larger table ##
library(data.table)

for (i in 1:22){
  assign(paste0("PGC_table",i),fread(paste0("PGC_table",i,"_new.txt")),envir = .GlobalEnv)
}
l = list()
for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0("PGC_table",i)))
}
combined_final_table <- rbindlist(l)
write.table (combined_final_table, file = "combined_PGC_table_with_CHR.POS_identifiers.txt", quote = F,row.names = F)

library(data.table)

for (i in 2:22){
  assign(paste0("Neurot_Assoc_Biobank_table",i),fread(paste0("Neurot_Assoc_Biobank_table",i,"_new.txt")),envir = .GlobalEnv)
}
l = list()

for (i in 2:22) {
  l[[i]] <- eval(parse(text = paste0("Neurot_Assoc_Biobank_table",i)))
}
combined_final_table <- rbindlist(l)
write.table (combined_final_table, file = "combined_Neurot_Assoc_table_with_CHR.POS_identifiers.txt", quote = F,row.names = F)