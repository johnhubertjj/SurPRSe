
# logistic regression calculation #
sig <-c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

library(data.table)
library(base)

## Read in covariates and fam file !may need to put fam file earlier
covariates <- fread("/Users/johnhubert/Dropbox/CLOZUK2.r7.select2PC.eigenvec")
fam2 <- fread("/Users/johnhubert/Dropbox/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")
Useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning_memory_conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities_coordination_movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")


## Create list for storage of p-values from logistic regression
x <- matrix(rep(NA,180), ncol = 12, nrow=15)
colnames(x) <- sig
rownames(x) <- Useful_pathways


## Calculate PRS using plink file format
for (l in 1:length(Useful_pathways)){
  for (i in 1:length(sig)) {
    
    
    # read in profiles
    PRS.profiles <-fread(paste0(Useful_pathways[l],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[l],"_",sig[i],".profile"))
    PRS.profiles2 <- merge(PRS.profiles,fam2,by="FID")
    PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles2,by.x="FID", by.y="FID", all = F)
    #PRS.Profiles.with.covariates_new_fam <- merge(PRS.Profiles.with.covariates, fam2, by.x = "FID", by.y = "FID", all = F)
    
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
    x[l,i] <- c(summary(model)$coefficients[2, "Pr(>|z|)"])
    #             res$NSNPs[i]<-max(dat$CNT)/2
  }#i
}

write.table(file=fout, res, col.names=T, row.names=F, quote=F, sep="\t")
x2 <- x
x2[,1:12] <- format.pval(x[,2],digits = 3,eps = 1e-80)
x2[,1] <- format.pval(x[,1],digits = 3)

write.csv(x3, file = "/Users/johnhubert/Dropbox/Pathway_analysis_background_selection_range.csv")


library(reshape2)
x3 <- melt(x)
x3$rowid <- c(1:15)
library(ggplot2)
ggplot(x3, aes(Var2,value, group = factor(Var1))) + geom_line(aes(color = factor(Var1)))
