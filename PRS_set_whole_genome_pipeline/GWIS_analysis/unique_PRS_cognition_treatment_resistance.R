# Cognition and treatment resistance in CLOZUK ##
library(data.table)

Fam_IDs <- fread("~/Documents/CLOZUK_cognition_METAL/COGS_IDs_V2.txt")
Phenotypes <- fread("~/Documents/CLOZUK_cognition_METAL/John_Phenotypes.txt")
fam_file <- fread("~/Documents/CLOZUK_cognition_METAL/COGSv2016_imputed/COGSv2016_IMPUTE2.fam")
colnames(fam_file) <- c("FID","IID","PID","MID","Sex","PHENO")

setnames(Phenotypes, old = "ID", new = "Study_ID")
New_phenotypes <- merge(Phenotypes, Fam_IDs, by = "Study_ID")

rows_to_keep_SCZ_unique <- New_phenotypes[DSM_diagnosis_new == 1 | DSM_diagnosis_new == 2 ]
rows_to_keep_SCZ_unique_SCZ_affective <- New_phenotypes[DSM_diagnosis_new == 1 | DSM_diagnosis_new == 2 | DSM_diagnosis_new == 3 | DSM_diagnosis_new == 4]


keep_file_SCZ_unique <- merge(fam_file, rows_to_keep_SCZ_unique, by =c("FID","IID"), all = F)
keep_file_SCZ_affective_unique <- merge(fam_file, rows_to_keep_SCZ_unique_SCZ_affective, by =c("FID","IID"), all = F)

setcolorder(keep_file_SCZ_unique,c(2,1,3:12))
setcolorder(keep_file_SCZ_affective_unique,c(2,1,3:12))

keep_file_SCZ_unique_cognition <- keep_file_SCZ_unique[,.(FID,IID,PID,MID,Sex,Composite_Zimpute)]
keep_file_SCZ_unique_treatment_resistance <-keep_file_SCZ_unique[,.(FID,IID,PID,MID,Sex,FINAL_TRS)]

keep_file_SCZ_unique_affective_cognition <- keep_file_SCZ_affective_unique[,.(FID,IID,PID,MID,Sex,Composite_Zimpute)]
keep_file_SCZ_unique_affective_treatment_resistance <-keep_file_SCZ_affective_unique[,.(FID,IID,PID,MID,Sex,FINAL_TRS)]

keep_file_SCZ_unique_cognition_rm_missing <- keep_file_SCZ_unique_cognition[Composite_Zimpute != 999]
keep_file_SCZ_unique_treatment_resistance_rm_missing <- keep_file_SCZ_unique_treatment_resistance[FINAL_TRS != 9]

keep_file_SCZ_unique_affective_cognition_rm_missing <- keep_file_SCZ_unique_affective_cognition[Composite_Zimpute != 999]
keep_file_SCZ_unique_affective_treatment_resistance_rm_missing <- keep_file_SCZ_unique_affective_treatment_resistance[FINAL_TRS != 9]

#keep_file_SCZ_unique <- keep_file_SCZ_unique[,.(FID,FAM_ID)]



setnames(PCA_matrix_2, old = "Individuals", new = "FID")


cognition_test <- merge(PCA_matrix_2, keep_file_SCZ_unique_cognition_rm_missing, by = "FID")
trtment_test <- merge(PCA_matrix_2, keep_file_SCZ_unique_treatment_resistance_rm_missing, by="FID")

residuals_to_save_cognition <- which(PCA_matrix_2$FID %in% cognition_test$FID)
residuals_to_save_trtment_resist <- which(PCA_matrix_2$FID %in% trtment_test$FID)

residuals_to_test_cog <- residuals_PRS_only[residuals_to_save_cognition]
residuals_to_test_trt <- residuals_PRS_only[residuals_to_save_trtment_resist]

alt_cog_test <- lm(cognition_test$Composite_Zimpute ~ residuals_to_test_cog)
null_cog_test <- lm(cognition_test$Composite_Zimpute ~ cognition_test$Schizophrenia)

alt_trt_test <- glm(trtment_test$FINAL_TRS ~ residuals_to_test_trt, family = binomial)
null_trt_test <- glm(trtment_test$FINAL_TRS ~ trtment_test$Schizophrenia, family = binomial)

##############################
cognition_test_affective <- merge(PCA_matrix_2, keep_file_SCZ_unique_affective_cognition_rm_missing, by = "FID")
trtment_test_affective <- merge(PCA_matrix_2, keep_file_SCZ_unique_affective_treatment_resistance_rm_missing, by="FID")

residuals_to_save_cognition_affective <- which(PCA_matrix_2$FID %in% cognition_test_affective$FID)
residuals_to_save_trtment_resist_affective <- which(PCA_matrix_2$FID %in% trtment_test_affective$FID)

residuals_to_test_cog_affective <- residuals_PRS_only[residuals_to_save_cognition_affective]
residuals_to_test_trt_affective <- residuals_PRS_only[residuals_to_save_trtment_resist_affective]

alt_cog_test_affective <- lm(cognition_test_affective$Composite_Zimpute ~ residuals_to_test_cog_affective)
null_cog_test_affective <- lm(cognition_test_affective$Composite_Zimpute ~ cognition_test_affective$Schizophrenia)

alt_trt_test_affective <- glm(trtment_test_affective$FINAL_TRS ~ residuals_to_test_trt_affective, family = binomial)
null_trt_test_affective <- glm(trtment_test_affective$FINAL_TRS ~ trtment_test_affective$Schizophrenia, family = binomial)

###############################
# PCA as a residual #

####### AFFECTIVE #######
residuals_to_test_cog_affective_PCA <- residuals_PRS[residuals_to_save_cognition_affective]
residuals_to_test_trt_affective_PCA <- residuals_PRS[residuals_to_save_trtment_resist_affective]

alt_cog_test_affective_PCA <- lm(cognition_test_affective$Composite_Zimpute ~ residuals_to_test_cog_affective_PCA)
null_cog_test_affective_PCA <- lm(cognition_test_affective$Composite_Zimpute ~ cognition_test_affective$Schizophrenia)

alt_trt_test_affective_PCA <- glm(trtment_test_affective$FINAL_TRS ~ residuals_to_test_trt_affective_PCA, family = binomial)
null_trt_test_affective_PCA <- glm(trtment_test_affective$FINAL_TRS ~ trtment_test_affective$Schizophrenia, family = binomial)

###### SCZ #######
residuals_to_test_cog_PCA <- residuals_PRS[residuals_to_save_cognition]
residuals_to_test_trt_PCA <- residuals_PRS[residuals_to_save_trtment_resist]

alt_cog_test_PCA <- lm(cognition_test$Composite_Zimpute ~ residuals_to_test_cog_PCA)
null_cog_test_PCA <- lm(cognition_test$Composite_Zimpute ~ cognition_test$Schizophrenia)

alt_trt_test_PCA <- glm(trtment_test$FINAL_TRS ~ residuals_to_test_trt_PCA, family = binomial)
null_trt_test_PCA <- glm(trtment_test$FINAL_TRS ~ trtment_test$Schizophrenia, family = binomial)


