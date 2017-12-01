
library(data.table)

# dont' change obj it makes the loop output a dataframe with the regression results
obj <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
obj2 <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0)
obj3 <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0)

tests <- c("Unique_alt_BP","Unique_null_BP", "Affective_alt_BP", "Affective_null_BP", "Unique_alt_PCA", "Unique_null_PCA", "Affective_alt_PCA", "Affective_null_PCA")

setwd("~/Documents/CLOZUK_cognition_METAL/")
significance_thresholds <- c(5e-08,1e-06,1e-04,0.05,0.01,0.1,0.2,0.5,1)

#Schizophrenia <- fread("/Volumes/PhD_storage/PGC_CLOZUK_output/PRS_scoring/PGC_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
for (i in significance_thresholds){
Schizophrenia <- fread(paste0("PRS_scoring_SCZ_CLOZUK_PGC2nocogs/CLOZUK_PGC2noCOGS_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_",i,".profile"))

#Bipolar <- fread("/Volumes/PhD_storage/BIP_CLOZUK_output/PRS_scoring/BIP_CLOZUK_whole_genome_significance_threshold_at_0.4.profile")
Bipolar <- fread(paste0("PRS_scoring_BIP/BIP_PGC2_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_",i,".profile"))

SCZ_BP_datasets <- list(Schizophrenia, Bipolar)


#Neuropsychiatric_datasets <- list(Schizophrenia,Bipolar,Educational_attainment, PGC_MDD, BIPvsSCZ, Neuroticism)

#covariates <- fread("~/Dropbox/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")
fam2 <- fread("~/Documents/CLOZUK_cognition_METAL/COGSv2016_imputed/COGSv2016_IMPUTE2.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")

#Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","POBI","T1DGC")
#Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","TWINSUK","1958BC")

for (l in 1:2){
  setkey(SCZ_BP_datasets[[l]],FID)
  # Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
  
  
  
  PRS.profiles.1 <- SCZ_BP_datasets[[l]]
  
  #PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles.1, by.x="FID", by.y="FID", all = F)
  #PRS.Profiles.with.covariates <- merge(PRS.Profiles.with.covariates, fam2, by.x = "FID", by.y = "FID", all = F)
  
  #  res$model[i]<-sig[i]
  
  # Calculate model including covariates against the polygenic risk score
  # model0<-glm(SCORE~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
  #m1<-mean(residuals(model0))
  #sd1<-sd(residuals(model0))
  
  # Calculate the Normalised score
  #PRS.Profiles.with.covariates$NORMSCORE<-(residuals(model0)-m1)/sd1
  PRS.profiles.1$NORMSCORE <- scale(PRS.profiles.1$SCORE, center = T, scale = T)
  #             hist(PRS.profiles$NORMSCORE)
  
  # change the Phenotypes so that they will work in a binary model
  PRS.profiles.1$PHENO.y <- PRS.profiles.1$PHENO.y - 1
  
  SCZ_BP_datasets[[l]] <- PRS.profiles.1
  
  #Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
}



PCA_matrix_2 <- data.frame(SCZ_BP_datasets[[1]]$FID, SCZ_BP_datasets[[1]]$PHENO, SCZ_BP_datasets[[1]]$NORMSCORE, SCZ_BP_datasets[[2]]$NORMSCORE)



names(PCA_matrix_2) <- c("Individuals","PHENOTYPE", "Schizophrenia", "Bipolar")
PCA_matrix_2$Colours <- "NA"




#write.csv(PCA_matrix, file = "/Volumes/PhD_storage/PRS_cross_disorder_table_optimised_thresholds.csv", col.names = T, row.names = F)



testing_reduce_controls <- prcomp(PCA_matrix_2[3:4], center = T)


# Cognition and treatment resistance in CLOZUK ##

Fam_IDs <- fread("~/Documents/CLOZUK_cognition_METAL/COGS_IDs_V2.txt")
Phenotypes <- fread("~/Documents/CLOZUK_cognition_METAL/John_Phenotypes.txt")


setnames(Phenotypes, old = "ID", new = "Study_ID")
New_phenotypes <- merge(Phenotypes, Fam_IDs, by = "Study_ID")

keep_file_SCZ <- merge(fam2, New_phenotypes, by =c("FID","IID"), all = F)
setnames(PCA_matrix_2, old = "Individuals", new = "FID")


cognition_test <- merge(PCA_matrix_2, keep_file_SCZ, by = "FID")

PRS_cognition_model_SZ <- lm(Composite_Zimpute ~ Schizophrenia, data = cognition_test)
PRS_cognition_model_BP <- lm(Composite_Zimpute ~ Bipolar, data = cognition_test)
PRS_treatment_resistance_model_SZ <- lm(FINAL_TRS ~ Schizophrenia, data = cognition_test)
PRS_treatment_resistance_model_BP <- lm(FINAL_TRS ~ Bipolar, data = cognition_test)


# this example if for a linear regression (the phenotype of interest is a quantiative trait)
# is using a discrete phenotype, a logistic regression needs to be run, and the code altered from 'lm' to 'glm' including the argument of 'family = binomial'
# alterations for the calculation of R2 will also need to be made using the command highlighted above



    tmp <- coef(summary(PRS_cognition_model_SZ))
    tmp2 <- summary(PRS_cognition_model_SZ)
   #hold <- summary(fit1)
    true_r2 <- tmp2$r.squared
    CI <- confint(PRS_cognition_model_SZ, level=0.95)
    CI <- CI[2,]
    tmp3 <- c("Schizophrenia",i,tmp[2,], true_r2, CI)
    obj <- rbind(obj, tmp3)
  
    tmp <- coef(summary(PRS_cognition_model_BP))
    tmp2 <- summary(PRS_cognition_model_BP)
    #hold <- summary(fit1)
    true_r2 <- tmp2$r.squared
    CI <- confint(PRS_cognition_model_BP, level=0.95)
    CI <- CI[2,]
    tmp3 <- c("Bipolar",i,tmp[2,], true_r2, CI)
    obj <- rbind(obj, tmp3)

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
    

    cognition_test_unique <- merge(PCA_matrix_2, keep_file_SCZ_unique_cognition_rm_missing, by = "FID")
    trtment_test_unique <- merge(PCA_matrix_2, keep_file_SCZ_unique_treatment_resistance_rm_missing, by="FID")
    

    
    #keep_file_SCZ_unique <- keep_file_SCZ_unique[,.(FID,FAM_ID)]
    
    
# Save residuals of the PCA
logistic_regress_3 <- glm(PCA_matrix_2$Schizophrenia ~ testing_reduce_controls$x[,"PC1"], family = gaussian)
residuals_PCA <- logistic_regress_3$residuals   

# Save residuals of Bipolar
logistic_regress_7 <- glm(PCA_matrix_2$Schizophrenia ~ PCA_matrix_2$Bipolar, family = gaussian )
residuals_PRS_only <- logistic_regress_7$residuals




# Include affective samples
##############################
cognition_test_affective <- merge(PCA_matrix_2, keep_file_SCZ_unique_affective_cognition_rm_missing, by = "FID")
trtment_test_affective <- merge(PCA_matrix_2, keep_file_SCZ_unique_affective_treatment_resistance_rm_missing, by="FID")

residuals_to_save_cognition <- which(PCA_matrix_2$FID %in% cognition_test_unique$FID)
residuals_to_save_trtment_resist <- which(PCA_matrix_2$FID %in% trtment_test_unique$FID)

residuals_to_save_cognition_affective <- which(PCA_matrix_2$FID %in% cognition_test_affective$FID)
residuals_to_save_trtment_resist_affective <- which(PCA_matrix_2$FID %in% trtment_test_affective$FID)

residuals_to_test_cog <- residuals_PRS_only[residuals_to_save_cognition]
residuals_to_test_trt <- residuals_PRS_only[residuals_to_save_trtment_resist]

residuals_to_test_cog_affective <- residuals_PRS_only[residuals_to_save_cognition_affective]
residuals_to_test_trt_affective <- residuals_PRS_only[residuals_to_save_trtment_resist_affective]

# PCA as a residual #
####### AFFECTIVE #######
residuals_to_test_cog_affective_PCA <- residuals_PCA[residuals_to_save_cognition_affective]
residuals_to_test_trt_affective_PCA <- residuals_PCA[residuals_to_save_trtment_resist_affective]

###### SCZ #######
residuals_to_test_cog_PCA <- residuals_PCA[residuals_to_save_cognition]
residuals_to_test_trt_PCA <- residuals_PCA[residuals_to_save_trtment_resist]

# Regressions ###

# unique_BP_residual#
alt_cog_test <- lm(cognition_test_unique$Composite_Zimpute ~ residuals_to_test_cog)
tmp <- coef(summary(alt_cog_test))
obj2 <- rbind(obj2, c(tests[1],i,tmp[2,]))
              
null_cog_test <- lm(cognition_test_unique$Composite_Zimpute ~ cognition_test_unique$Schizophrenia)
tmp <- coef(summary(null_cog_test))
obj2 <- rbind(obj2,c(tests[2],i,tmp[2,]))

alt_trt_test <- glm(trtment_test$FINAL_TRS ~ residuals_to_test_trt, family = binomial)
tmp <- coef(summary(alt_cog_test))
obj3 <- rbind(obj3, c(tests[1],i,tmp[2,]))

null_trt_test <- glm(trtment_test$FINAL_TRS ~ trtment_test$Schizophrenia, family = binomial)
tmp <- coef(summary(null_trt_test))
obj3 <- rbind(obj3,c(tests[2],i,tmp[2,]))


#Affective_BP_residual
alt_cog_test_affective <- lm(cognition_test_affective$Composite_Zimpute ~ residuals_to_test_cog_affective)
tmp <- coef(summary(alt_cog_test_affective))
obj2 <- rbind(obj2, c(tests[3],i,tmp[2,]))

null_cog_test_affective <- lm(cognition_test_affective$Composite_Zimpute ~ cognition_test_affective$Schizophrenia)
tmp <- coef(summary(null_cog_test_affective))
obj2 <- rbind(obj2, c(tests[4],i,tmp[2,]))

alt_trt_test_affective <- glm(trtment_test_affective$FINAL_TRS ~ residuals_to_test_trt_affective, family = binomial)
tmp <- coef(summary(alt_trt_test_affective))
obj3 <- rbind(obj3, c(tests[3],i,tmp[2,]))

null_trt_test_affective <- glm(trtment_test_affective$FINAL_TRS ~ trtment_test_affective$Schizophrenia, family = binomial)
tmp <- coef(summary(null_trt_test_affective))
obj3 <- rbind(obj3, c(tests[4],i,tmp[2,]))

#Unique_PCA
alt_cog_test_PCA <- lm(cognition_test_unique$Composite_Zimpute ~ residuals_to_test_cog_PCA)
tmp <- coef(summary(alt_cog_test_PCA))
obj2 <- rbind(obj2, c(tests[5],i,tmp[2,]))

null_cog_test_PCA <- lm(cognition_test_unique$Composite_Zimpute ~ cognition_test_unique$Schizophrenia)
tmp <- coef(summary(null_cog_test_PCA))
obj2 <- rbind(obj2, c(tests[6],i,tmp[2,]))

alt_trt_test_PCA <- glm(trtment_test$FINAL_TRS ~ residuals_to_test_trt_PCA, family = binomial)
tmp <- coef(summary(alt_trt_test_PCA))
obj3 <- rbind(obj3, c(tests[5],i,tmp[2,]))

null_trt_test_PCA <- glm(trtment_test$FINAL_TRS ~ trtment_test$Schizophrenia, family = binomial)
tmp <- coef(summary(null_trt_test_PCA))
obj3 <- rbind(obj3, c(tests[6],i,tmp[2,]))

# Affective_PCA
alt_cog_test_affective_PCA <- lm(cognition_test_affective$Composite_Zimpute ~ residuals_to_test_cog_affective_PCA)
tmp <- coef(summary(alt_cog_test_affective_PCA))
obj2 <- rbind(obj2, c(tests[7],i,tmp[2,]))

null_cog_test_affective_PCA <- lm(cognition_test_affective$Composite_Zimpute ~ cognition_test_affective$Schizophrenia)
tmp <- coef(summary(null_cog_test_affective_PCA))
obj2 <- rbind(obj2, c(tests[8],i,tmp[2,]))

alt_trt_test_affective_PCA <- glm(trtment_test_affective$FINAL_TRS ~ residuals_to_test_trt_affective_PCA, family = binomial)
tmp <- coef(summary(alt_trt_test_affective_PCA))
obj3 <- rbind(obj3, c(tests[7],i,tmp[2,]))

null_trt_test_affective_PCA <- glm(trtment_test_affective$FINAL_TRS ~ trtment_test_affective$Schizophrenia, family = binomial)
tmp <- coef(summary(null_trt_test_affective_PCA))
obj3 <- rbind(obj3, c(tests[8],i,tmp[2,]))
}

results_cognition <- as.data.table(obj2)
results_trt <- as.data.table(obj3)

# Cognition #
Results_PCA_cognition_0.05 <- results_cognition[grep(pattern = "PCA",results_cognition$test)]
Results_PCA_cognition_0.05 <- Results_PCA_cognition_0.05[score == "0.05"]

Results_BPresid_cognition_0.05 <- results_cognition[grep(pattern = "BP",results_cognition$test)]
Results_BPresid_cognition_0.05 <- Results_BPresid_cognition_0.05[score == "0.05"]

Results_PCA_cognition_0.5 <- results_cognition[grep(pattern = "PCA",results_cognition$test)]
Results_PCA_cognition_0.5 <- Results_PCA_cognition_0.5[score == "0.5"]

Results_BPresid_cognition_0.5 <- results_cognition[grep(pattern = "BP",results_cognition$test)]
Results_BPresid_cognition_0.5 <- Results_BPresid_cognition_0.5[score == "0.5"]

# Treatment_resistance #
Results_PCA_trt_0.05 <- results_trt[grep(pattern = "PCA",results_trt$test)]
Results_PCA_trt_0.05 <- Results_PCA_trt_0.05[score == "0.05"]

Results_BPresid_trt_0.05 <- results_trt[grep(pattern = "BP",results_trt$test)]
Results_BPresid_trt_0.05 <- Results_BPresid_trt_0.05[score == "0.05"]

Results_PCA_trt_0.5 <- results_trt[grep(pattern = "PCA",results_trt$test)]
Results_PCA_trt_0.5 <- Results_PCA_trt_0.5[score == "0.5"]

Results_BPresid_trt_0.5 <- results_trt[grep(pattern = "BP",results_trt$test)]
Results_BPresid_trt_0.5 <- Results_BPresid_trt_0.5[score == "0.5"]

# this is a clean-up step - do not change
results <- obj[which(obj$score %in% significance_thresholds),]


#write.table(results,file = "PRS_threshold_selection_on_cognition_for_SCZandBP_COGS_CLOZUKPGC2noCOGS.txt", quote = F,row.names = F)


##### Testing models in likelihood test #####

logistic_regress_1 <- glm (cognition_test$Composite_Zimpute ~ testing_reduce_controls$x[,"PC1"], family = gaussian)
logistic_regress_2 <- glm (cognition_test$Composite_Zimpute ~ cognition_test$Schizophrenia + testing_reduce_controls$x[,"PC1"], family = gaussian)

library(epiDisplay)
lrtest(logistic_regress_1,logistic_regress_2)

logistic_regress_5 <- glm(cognition_test$Composite_Zimpute ~ cognition_test$Bipolar, family = gaussian)
logistic_regress_6 <- glm(cognition_test$Composite_Zimpute ~ cognition_test$Schizophrenia + cognition_test$Bipolar, family = gaussian)
lrtest(logistic_regress_5,logistic_regress_6)


logistic_regress_4 <- glm(PCA_matrix_2$PHENOTYPE ~ residuals_PRS , family = binomial)




logistic_regress_7 <- glm(PCA_matrix_2$Schizophrenia ~ PCA_matrix_2$Bipolar, family = gaussian )
residuals_PRS_only <- logistic_regress_7$residuals

logistic_regress_8 <- glm(PCA_matrix_2$PHENOTYPE ~ residuals_PRS_only, family = binomial )


}

residualPlots(logistic_regress_1)


type2 <- data.frame(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,"PC1"], residuals)
names(type2) <- c("Phenotype", "PC1", "residuals")

logistic_regress_2 <- lm(Phenotype ~ PC1 + residuals, data = type2) 

model1 <- logistic_regress_1
model2 <- logistic_regress_2

model1 <- logistic_regress_5
model2 <- logistic_regress_6

l0 <- deviance(model1)
df0 <- df.residual(model1)

l1 <- deviance(model2)
df1 <- df.residual(model2)

degf <- df0-df1
if (degf > 0) pvalue <- 1-pchisq(l0-l1, df0-df1)

library(epiDisplay)
lrtest(model1,model2)

logistic_regress_3 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Bipolar, binomial(link = 'logit'))
logistic_regress_5 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Schizophrenia , binomial(link = 'logit'))
logistic_regress_4 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Bipolar + PCA_matrix_2$Schizophrenia , binomial(link = 'logit'))




for (i in 1:8){
  t <- summary(eval(parse(text = paste0("logistic_regress_",i))))$coefficients
  print(t)
}

