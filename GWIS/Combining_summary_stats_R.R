# Analysis of COGS training on CLOZUK_PGC2-COGS
library(data.table)
Chromosomes_to_split <- seq(1:22)

## Read in the latest Summary stats tables after converting to one table
for (i in Chromosomes_to_split){
  assign(paste0("Training_table", i), fread(paste0("CLOZUK2_noPGC2_chr",i,".assoc.dosage")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_split) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}
combined_final_table <- rbindlist(l)

write.table (combined_final_table, file = paste0("CLOZUK2_noPGC2.assoc.dosage"), quote = F, row.names = F)

summ_stat_1 <- fread("CLOZUK2_noPGC2.assoc.dosage")
summ_stat_2 <- fread("daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt")

### Replacing PGC values###
PGC.integers.to.change <- summ_stat_1[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T, invert = T)]]
PGC.divisive.names <- PGC.data.frame[PGC.integers.to.change]
PGC.alternative <- PGC.data.frame[,SNP := paste0(CHR,":",BP)]
PGC.data.frame <- PGC.data.frame[,SNP := paste0(CHR,":",BP)]
#PGC.alternative <- PGC.alternative[,c("CHR","SNP","BP","A1","A2","OR"),with = F]
PGC.alternative <- PGC.alternative[, Column_headers, with = F]
