## Create a larger table for CLOZUK ##

for (i in 1:22){
  assign(paste0("CLOZUK_PGC_CLUMPED_FINAL_chr",i),fread(paste0("./output/CLOZUK_PGC_CLUMPED_FINAL_chr",chromosome.number,".txt")),envir = .GlobalEnv)
}
l = list()

for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0("CLOZUK_PGC_CLUMPED_FINAL_chr",i)))
}

combined_final_table <- rbindlist(l)
write.table (combined_final_table, file = "combined_CLOZUK_PGC_CLUMPED_FINAL.txt", quote = F,row.names = F)