GenSCOT <- fread("BPSCZ.bp_v_scz.results.txt")

map_file_for_liftover <- GenSCOT[,.(CHR, SNP, 0, BP)]

# Ran the first part of liftover script to get map-ish file with rs numbers intact
# if scripting, provide an ATP link to the software in question, will run a lot faster

Bed_file_for_liftover <- GenSCOT[,BP0 := (BP-1)]

# attempted substitution for other software
Bed_file_for_liftover$CHR <- sub("^", "chr", Bed_file_for_liftover$CHR)
Bed_file_for_liftover <- Bed_file_for_liftover[,.(CHR,BP0, BP, SNP, A1, A2)]
setcolorder(Bed_file_for_liftover, c(1,11,3,2,4:10))
write.table(Bed_file_for_liftover, file = "BIPvsSCZ_bed_file_for_liftover.txt", quote = F, row.names = F, col.names = F)
write.table(map_file_for_liftover, file = "BIPvsSCZ_map_file_for_liftover.txt", quote = F, row.names = F, col.names = F)


new_bp_positions <- fread("~/Downloads/remapped_BIPvsSCZ_map_file_for_liftover.txt.bed")
setnames(new_bp_positions, c("CHR", "BP_dont_use", "BP_use", "SNP"))
new_bp_positions$CHR <- sub("chr", "", new_bp_positions$CHR)
new_bp_positions$CHR <- as.integer(new_bp_positions$CHR)
Hg19_BIPvsSCZ_raw_data <- merge(new_bp_positions, GenSCOT, by = c("SNP", "CHR") , all = F)
Hg19_BIPvsSCZ_raw_data <- Hg19_BIPvsSCZ_raw_data[,.(CHR,SNP,BP_use,A1,A2,INFO,OR,SE,P,ngt)]
setnames(Hg19_BIPvsSCZ_raw_data, "BP_use", "BP")
head(Hg19_BIPvsSCZ_raw_data)
write.table(Hg19_BIPvsSCZ_raw_data, file = "BPSCZ.bp_v_scz.results_GR37.p13.txt", quote = F, row.names = F, col.names = T)
