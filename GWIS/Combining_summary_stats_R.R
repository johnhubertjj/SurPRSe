# Analysis of COGS training on CLOZUK_PGC2-COGS
library(data.table)
Chromosomes_to_split <- seq(1:22)

## Read in the latest Summary stats tables after converting to one table
#for (i in Chromosomes_to_split){
#  assign(paste0("Training_table", i), fread(paste0("CLOZUK2_noPGC2_chr",i,".assoc.dosage")),envir = .GlobalEnv)
#}
#l = list()

## print out to one table under a common filename
#for (i in Chromosomes_to_split) {
#  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
#}
#combined_final_table <- rbindlist(l)


## Read in the latest Summary stats tables after converting to one table
for (i in Chromosomes_to_split){
  assign(paste0("Training_table", i), fread(paste0("./RSupdate/CLOZUK_chr",i,".KAVIAR.update")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_split) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}
combined_final_table <- rbindlist(l)


summ_stat_1 <- fread("CLOZUK2_noPGC2.assoc.dosage")
colnames
summ_stat_RS <- merge(combined_final_table)
