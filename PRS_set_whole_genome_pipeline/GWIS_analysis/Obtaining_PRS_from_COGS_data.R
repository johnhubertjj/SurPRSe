library(data.table)
significance_thresholds <- c(5e-08,1e-06,0.05,0.01,0.1,0.2,0.3,0.5,0.1)

setwd("~/Documents/CLOZUK_cognition_METAL/")

for (i in 1:length(significance_thresholds)){
#Schizophrenia <- fread("/Volumes/PhD_storage/PGC_CLOZUK_output/PRS_scoring/PGC_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
Schizophrenia <- fread(paste0("PRS_scoring_SCZ_CLOZUK_PGC2nocogs/CLOZUK_PGC2noCOGS_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_",significance_thresholds[i],".profile"))

#Bipolar <- fread("/Volumes/PhD_storage/BIP_CLOZUK_output/PRS_scoring/BIP_CLOZUK_whole_genome_significance_threshold_at_0.4.profile")
Bipolar <- fread(paste0("PRS_scoring_BIP/BIP_PGC2_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_",significance_thresholds[i]".profile"))

SCZ_BP_datasets <- list(Schizophrenia, Bipolar)

#Neuropsychiatric_datasets <- list(Schizophrenia,Bipolar,Educational_attainment, PGC_MDD, BIPvsSCZ, Neuroticism)

# covariates <- fread("~/Dropbox/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")
fam2 <- fread("~/Documents/CLOZUK_cognition_METAL/COGSv2016_imputed/COGSv2016_IMPUTE2.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")

#Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","POBI","T1DGC")
#Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","TWINSUK","1958BC")

for (i in 1:2){
  setkey(SCZ_BP_datasets[[i]],FID)
  # Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
  
  
  
  PRS.profiles.1 <- SCZ_BP_datasets[[i]]
  
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
  
  SCZ_BP_datasets[[i]] <- PRS.profiles.1
  
  #Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
}



PCA_matrix_2 <- data.frame(SCZ_BP_datasets[[1]]$FID, SCZ_BP_datasets[[1]]$PHENO, SCZ_BP_datasets[[1]]$NORMSCORE, SCZ_BP_datasets[[2]]$NORMSCORE)



names(PCA_matrix_2) <- c("Individuals","PHENOTYPE", "Schizophrenia", "Bipolar")
PCA_matrix_2$Colours <- "NA"




#write.csv(PCA_matrix, file = "/Volumes/PhD_storage/PRS_cross_disorder_table_optimised_thresholds.csv", col.names = T, row.names = F)



testing_reduce_controls <- prcomp(PCA_matrix_2[3:4], center = T)




##### Separating out risk of SCZ and BIP #####

#logistic_regress_1 <- glm (PCA_matrix_2$PHENOTYPE ~ testing_reduce_controls$x[,"PC1"], family = binomial)
#logistic_regress_2 <- glm (PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Schizophrenia + testing_reduce_controls$x[,"PC1"], family = binomial)

logistic_regress_3 <- glm(PCA_matrix_2$Schizophrenia ~ testing_reduce_controls$x[,"PC1"], family = gaussian)
residuals_PRS <- logistic_regress_3$residuals 

#logistic_regress_4 <- glm(PCA_matrix_2$PHENOTYPE ~ residuals_PRS , family = binomial)


#logistic_regress_5 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Bipolar, family = binomial)
#logistic_regress_6 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Schizophrenia + PCA_matrix_2$Bipolar, family = binomial)

logistic_regress_7 <- glm(PCA_matrix_2$Schizophrenia ~ PCA_matrix_2$Bipolar, family = gaussian )
residuals_PRS_only <- logistic_regress_7$residuals

#logistic_regress_8 <- glm(PCA_matrix_2$PHENOTYPE ~ residuals_PRS_only, family = binomial )




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

