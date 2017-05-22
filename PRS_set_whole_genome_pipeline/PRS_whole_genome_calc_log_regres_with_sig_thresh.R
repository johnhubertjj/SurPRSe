### Start Timer
ptm <- proc.time()

### Library
library(data.table)

## set wd
setwd(".")

########################################
# adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
getwd()
print(args)

# Specify the different input tables #
Training_name <- args[3]
Validation_name <- args [4]
significance_thresholds <- as.numeric(args[c(5:length(args))])
print(significance_thresholds)

## Read in covariates and fam file !may need to put fam file earlier
#covariates <- fread("/Users/JJ/Documents/PhD_clumping/Profiles/CLOZUK2.r7.select2PC.eigenvec")
covariates <- fread("/Volumes/HD-PCU2/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")
fam2 <- fread("/Volumes/HD-PCU2/Stationary_data/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")


## Create list for storage of p-values from logistic regression
residuals <- matrix(data=rep(NA,5*length(significance_thresholds)), ncol = 5, nrow = length(significance_thresholds))

## Calculate PRS using plink file format
for (i in 1:length(significance_thresholds)) {
  PRS.profiles <- fread(paste0("./",Training_name,"_",Validation_name,"_output/PRS_scoring/",Training_name,"_",Validation_name,"_whole_genome_significance_threshold_at_",significance_thresholds[i],".profile"))
    
  PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles,by.x="FID", by.y="FID", all = F)
  PRS.Profiles.with.covariates <- merge(PRS.Profiles.with.covariates, fam2, by.x = "FID", by.y = "FID", all = F)
    
  #  res$model[i]<-sig[i]
    
  # Calculate model including covariates against the polygenic risk score
  model0<-glm(SCORE~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
  m1<-mean(residuals(model0))
  sd1<-sd(residuals(model0))
    
  # Calculate the Normalised score
  PRS.Profiles.with.covariates$NORMSCORE<-(residuals(model0)-m1)/sd1
    
  #             hist(PRS.profiles$NORMSCORE)
    
  # change the Phenotypes so that they will work in a binary model
  PRS.Profiles.with.covariates$PHENO.y <- PRS.Profiles.with.covariates$PHENO.y - 1
  model<-glm(PHENO.y~PRS.Profiles.with.covariates$NORMSCORE, family=binomial, data=PRS.Profiles.with.covariates)
    
    #             res$Effect[i]<-summary(model)$coefficients[2, "Estimate"]
    #             res$SE[i]<-summary(model)$coefficients[2, "Std. Error"]
    residuals[i,] <- c(significance_thresholds[i], summary(model)$coefficients[2, "Pr(>|z|)"], summary(model)$coefficients[2,"z value"], summary(model)$coefficients[2,"Std. Error"], summary(model)$coefficients[2,"Estimate"])
    #             res$NSNPs[i]<-max(dat$CNT)/2
  }#i

colnames(residuals) <- c("Pval_threshold", paste0(Validation_name,"_", Training_name,"_PRS_p_value"), "Z_value", "Std.Error")

write.csv(residuals, file = paste0("./",Training_name,"_",Validation_name,"_output/PRS_scoring/",Training_name,"_",Validation_name,"_PRS_residuals_using_fam_file.csv"), row.names = F)

