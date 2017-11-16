# Cognition and treatment resistance in CLOZUK ##
library(data.table)

Samples_to_keep <- fread("~/Documents/Cardiff_cogs_treatment_resis_data/John_pheno_data.txt")
fam_file <- fread("~/Dropbox/Stationary_data/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam_file) <- c("FID","FAM_ID","GD","Pheno0.5","PHeno1", "Pheno2")

rows_to_keep_SCZ_unique <- Samples_to_keep[DSM_diagnosis_new == 1]
rows_to_keep_SCZ_total <- Samples_to_keep[DSM_diagnosis_new == 1 | DSM_diagnosis_new == 2 | DSM_diagnosis_new == 3 | DSM_diagnosis_new == 4]
rows_to_keep_schizoaffective <- Samples_to_keep[DSM_diagnosis_new == 2 | DSM_diagnosis_new == 3]

keep_file_SCZ_unique <- merge(fam_file, rows_to_keep_SCZ_unique, by ="FAM_ID", all = F)
setcolorder(keep_file_SCZ_unique,c(2,1,3:12))
keep_file_SCZ_unique_cognition <- keep_file_SCZ_unique[,.(FID,FAM_ID,GD,Pheno0.5,PHeno1,Composite_Zimpute)]
keep_file_SCZ_unique_treatment_resistance <-keep_file_SCZ_unique[,.(FID,FAM_ID,GD,Pheno0.5,PHeno1,FINAL_TRS)]

keep_file_SCZ_unique_cognition_rm_missing <- keep_file_SCZ_unique_cognition[Composite_Zimpute != 999]
keep_file_SCZ_unique_treatment_resistance_rm_missing <- keep_file_SCZ_unique_treatment_resistance[FINAL_TRS != 9]

keep_file_SCZ_unique <- keep_file_SCZ_unique[,.(FID,FAM_ID)]

#write.table(keep_file_SCZ_unique_cognition_rm_missing, file = "~/Dropbox/Stationary_data/CLOZUK.r7.GWAS_IDs_phenotype_cognition_COGS_ONLY_UNIQUE_SCZ_ONLY_DSM_1.fam",row.names = F, col.names = F, quote = F)

load("~/Dropbox/Separating_BP_and_SCZ_risk.RData")

setnames(keep_file_SCZ_unique_cognition_rm_missing, old = "FAM_ID", new = "IID")
setnames(keep_file_SCZ_unique_treatment_resistance_rm_missing, old = "FAM_ID", new = "IID")
setnames(PCA_matrix_2, old = "Individuals", new = "FID")


cognition_test <- merge(PCA_matrix_2, keep_file_SCZ_unique_cognition_rm_missing, by = "FID")
trtment_test <- merge(PCA_matrix_2, keep_file_SCZ_unique_treatment_resistance_rm_missing, by="FID")

residuals_to_save_cognition <- which(PCA_matrix_SCZ_BP$FID %in% cognition_test$FID)
residuals_to_save_trtment_resist <- which(PCA_matrix_SCZ_BP$FID %in% trtment_test$FID)

residuals_to_test_cog <- residuals_PRS_only[residuals_to_save_cognition]
residuals_to_test_trt <- residuals_PRS_only[residuals_to_save_trtment_resist]

alt_cog_test <- lm(cognition_test$Composite_Zimpute ~ residuals_to_test_cog)
null_cog_test <- lm(cognition_test$Composite_Zimpute ~ cognition_test$Schizophrenia)

alt_trt_test <- glm(trtment_test$FINAL_TRS ~ residuals_to_save_trtment_resist, family = binomial)
null_trt_test <- glm(trtment_test$FINAL_TRS ~ trtment_test$Schizophrenia, family = binomial)

print(summary(alt_cog_test), summary(null_cog_test), summary(alt_trt_test), summary(null_trt_test))
