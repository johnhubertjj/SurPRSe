testing_GeM <- fread("GeM_FULL_GENOME_consensus_with_PGCnoClozuk_flipped_alleles_no_duplicates.bim")
pgc_full <- fread("PGC_table_full.txt")
setkey(pgc_full,SNP)
