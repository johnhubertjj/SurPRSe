# Preparation of polygenic risk score for COGS # 

COGS <- fread("~/Documents/CLOZUK_cognition_METAL/COGSv2016_imputed/COGSv2016_IMPUTE2.bim")
CLOZUK_PGC2noCOGS <- fread("~/Dropbox/CLOZUK_PGC2_noCOGS1.tbl")
Extra_info <- fread("~/Dropbox/PGC2noCLOZUK_whole_genome_flipped_RS_SNP.txt")
setnames(CLOZUK_PGC2noCOGS, old = c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value"), new = c("RS_SNP", "A1", "A2", "BETA", "SE", "P"))
rm(COGS)
final_summ_stats <- merge(CLOZUK_PGC2noCOGS,Extra_info, by = c("RS_SNP"))
setnames(final_summ_stats, old = c("RS_SNP", "A1.x", "A2.x","SE.x" ,"P.x"), new = c("SNP", "A1","A2","SE","P"))

final_summ_stats <- final_summ_stats[,.(CHR,BP,SNP,A1,A2,BETA,SE,P)]

write.table(final_summ_stats,file = "~/Documents/CLOZUK_cognition_METAL/CLOZUK_PGC2_noCOGSmetaanalysis_flipped_INFOabove0point9.txt", quote = F, row.names = F, col.names = T)
