# Combining COGS phentypes #

Famfile <- fread("~/Documents/CLOZUK_cognition_METAL/COGSv2016_imputed/COGSv2016_IMPUTE2.fam")
Phenotype_file <- fread("~/Documents/Copy_of_Cardiff_g_calculation_Original_Data.csv")
setnames(Famfile,old = "V2", new = "STUDY.ID")
new_phenotype <- merge(Famfile,Phenotype_file, by = "STUDY.ID")
