
library(data.table)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

setwd("~/Documents/CLOZUK_cognition_METAL/")

#Schizophrenia <- fread("/Volumes/PhD_storage/PGC_CLOZUK_output/PRS_scoring/PGC_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
Schizophrenia_point_zero_five <- fread("PRS_scoring_SCZ_CLOZUK_PGC2nocogs/CLOZUK_PGC2noCOGS_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_0.05.profile")
Schizophrenia_point_five <- fread("PRS_scoring_SCZ_CLOZUK_PGC2nocogs/CLOZUK_PGC2noCOGS_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_0.5.profile")

#Bipolar <- fread("/Volumes/PhD_storage/BIP_CLOZUK_output/PRS_scoring/BIP_CLOZUK_whole_genome_significance_threshold_at_0.4.profile")
Bipolar_point_zero_five <- fread("PRS_scoring_BIP/BIP_PGC2_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_0.05.profile")
Bipolar_point_five <- fread("PRS_scoring_BIP/BIP_PGC2_COGSv2016_IMPUTE2_whole_genome_significance_threshold_at_0.5.profile")

SCZ_BP_datasets_point_zero_five <- list(Schizophrenia_point_zero_five, Bipolar_point_zero_five)
SCZ_BP_datasets_point_five <- list (Schizophrenia_point_five, Bipolar_point_five)

#Neuropsychiatric_datasets <- list(Schizophrenia,Bipolar,Educational_attainment, PGC_MDD, BIPvsSCZ, Neuroticism)

covariates <- fread("~/Dropbox/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")
fam2 <- fread("~/Documents/CLOZUK_cognition_METAL/COGSv2016_imputed/COGSv2016_IMPUTE2.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")

Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","POBI","T1DGC")
#Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","TWINSUK","1958BC")

for (i in 1:2){
  setkey(SCZ_BP_datasets_point_zero_five[[i]],FID)
  # Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
  
  
  
  PRS.profiles.1 <- SCZ_BP_datasets_point_zero_five[[i]]
  
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
  
  SCZ_BP_datasets_point_zero_five[[i]] <- PRS.profiles.1
  
  #Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
}
for (i in 1:2){
  setkey(SCZ_BP_datasets_point_five[[i]],FID)
  # Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
  
  
  
  PRS.profiles.1 <- SCZ_BP_datasets_point_five[[i]]
  
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
  
  SCZ_BP_datasets_point_five[[i]] <- PRS.profiles.1
  
  #Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
}


PCA_matrix_2 <- data.frame(SCZ_BP_datasets_point_five[[1]]$FID, SCZ_BP_datasets_point_five[[1]]$PHENO.y, SCZ_BP_datasets_point_five[[1]]$NORMSCORE, SCZ_BP_datasets_point_five[[2]]$NORMSCORE)

PCA_matrix_2 <- data.frame(SCZ_BP_datasets_point_zero_five[[1]]$FID, SCZ_BP_datasets_point_zero_five[[1]]$PHENO.y, SCZ_BP_datasets_point_zero_five[[1]]$NORMSCORE, SCZ_BP_datasets_point_zero_five[[2]]$NORMSCORE)



names(PCA_matrix_2) <- c("Individuals","PHENOTYPE", "Schizophrenia", "Bipolar")
PCA_matrix_2$Colours <- "NA"




#write.csv(PCA_matrix, file = "/Volumes/PhD_storage/PRS_cross_disorder_table_optimised_thresholds.csv", col.names = T, row.names = F)



testing_reduce_controls <- prcomp(PCA_matrix_2[3:4], center = T)




##### Separating out risk of SCZ and BIP #####

logistic_regress_1 <- glm (PCA_matrix_2$PHENOTYPE ~ testing_reduce_controls$x[,"PC1"], family = binomial)
logistic_regress_2 <- glm (PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Schizophrenia + testing_reduce_controls$x[,"PC1"], family = binomial)

logistic_regress_3 <- glm(PCA_matrix_2$Schizophrenia ~ testing_reduce_controls$x[,"PC1"], family = gaussian)

residuals_PRS <- logistic_regress_3$residuals 

logistic_regress_4 <- glm(PCA_matrix_2$PHENOTYPE ~ residuals_PRS , family = binomial)


logistic_regress_5 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Bipolar, family = binomial)
logistic_regress_6 <- glm(PCA_matrix_2$PHENOTYPE ~ PCA_matrix_2$Schizophrenia + PCA_matrix_2$Bipolar, family = binomial)

logistic_regress_7 <- glm(PCA_matrix_2$Schizophrenia ~ PCA_matrix_2$Bipolar, family = gaussian )
residuals_PRS_only <- logistic_regress_7$residuals

logistic_regress_8 <- glm(PCA_matrix_2$PHENOTYPE ~ residuals_PRS_only, family = binomial )




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

