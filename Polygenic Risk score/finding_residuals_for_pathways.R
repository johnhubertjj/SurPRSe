
# logistic regression calculation #
sig <-c(1,0.05)

library(data.table)
library(base)

## Read in covariates and fam file !may need to put fam file earlier
covariates <- fread("/Users/JJ/Documents/PhD_clumping/Profiles/CLOZUK2.r7.select2PC.eigenvec")
fam2 <- fread("/Users/JJ/Documents/PhD_clumping/Profiles/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")
Useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")


## Create list for storage of p-values from logistic regression
x <- matrix(rep(NA,26), ncol = 2, nrow=13)
colnames(x) <- sig
rownames(x) <- Useful_pathways


memory.limit(20000)
## Calculate PRS using plink file format
for (l in 1:length(Useful_pathways)){
  for (i in 1:length(sig)) {
    
    
    # read in profiles
    PRS.profiles <-fread(paste0("./pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[l],"_",sig[i],".profile"))
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

for (i in 1:length(sig)){
  write.table(residuals[[i]], file = paste0(sig[i],"CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file.txt"), quote = F, row.names = F)
}
mtext("FT4-Replication. MAF>0.1", side = 1, line = -1, outer = TRUE)
write.table(file=fout, res, col.names=T, row.names=F, quote=F, sep="\t")
x2 <- x
x2[,2] <- format.pval(x[,2],digits = 3,eps = 1e-80)
x2[,1] <- format.pval(x[,1],digits = 3)
write.csv(x2, file = "Pathway_analysis_background_selection.csv")
