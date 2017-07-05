# Set the options so the data is output without exp. notation
options(scipen = 999)

# read in relevant files (previously problems with names, maybe plug in somewhere in CLOZUK_PGC_COMBINE?) # check for rsID's etc...
GenSCOT <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/pgc.scz.full.2012-04.txt")
setnames(GenSCOT,c("SNP","CHR","BP","A1","A2","OR","SE","P","INFO","ngt","CEUaf"))
map_file_for_liftover <- GenSCOT[,.(CHR, SNP, 0, BP)]

# Ran the first part of liftover script to get map-ish file with rs numbers intact
# if scripting, provide an API link to the software in question, will run a lot faster

Bed_file_for_liftover <- GenSCOT[,BP0 := (BP-1)]

# attempted substitution for other software
Bed_file_for_liftover$CHR <- sub("^", "chr", Bed_file_for_liftover$CHR)
Bed_file_for_liftover <- Bed_file_for_liftover[,.(CHR,BP0, BP, SNP)]
write.table(Bed_file_for_liftover, file = "PGC1_hg18_bed_file_for_liftover.txt", quote = F, row.names = F, col.names = F)
write.table(map_file_for_liftover, file = "BIPvsSCZ_map_file_for_liftover.txt", quote = F, row.names = F, col.names = F)

# Some API link to UCSC (probably the best to use here...)

new_bp_positions <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/hglft_PGC1_hg19.bed")
setnames(new_bp_positions, c("CHR", "BP_dont_use", "BP_use", "SNP"))
new_bp_positions$CHR <- sub("chr", "", new_bp_positions$CHR)
new_bp_positions$CHR <- as.integer(new_bp_positions$CHR)
Hg19_BIPvsSCZ_raw_data <- merge(new_bp_positions, GenSCOT, by = c("SNP", "CHR") , all = F)
Hg19_BIPvsSCZ_raw_data <- Hg19_BIPvsSCZ_raw_data[,.(CHR,SNP,BP_use,A1,A2,INFO,OR,SE,P,ngt,CEUaf)]
setnames(Hg19_BIPvsSCZ_raw_data, "BP_use", "BP")
head(Hg19_BIPvsSCZ_raw_data)
write.table(Hg19_BIPvsSCZ_raw_data, file = "HG19_pgc.scz.full.2012-04.txt", quote = F, row.names = F, col.names = T)
