data_plot <- fread("~/Documents/CLOZUK_PGC2_pathway_analysis/SCHIZ_gene-sets_analysis/FINAL_PATHWAY_RESULTS_PRS_PROFILESPGC2_no_CLOZUK_CLOZUK_no_PGC2_Selected_Pocklington_plus_GO_pathways_SCHIZ_extended.txt")
data_plot_2_whole_genome_full <-fread("~/Documents/CLOZUK_PGC2_pathway_analysis/CLOZUKminusPGC2_PGC2minusCLOZUK_whole_genome_PRS/Full_genome/PRS_scoring/PGC2minusCLOZUK_CLOZUKminusPGC2_whole_genome_significance_threshold_at_0.05.profile")
data_plot_2_whole_genome_genic <- fread("~/Documents/CLOZUK_PGC2_pathway_analysis/CLOZUKminusPGC2_PGC2minusCLOZUK_whole_genome_PRS/Genic_PRS/PRS_scoring/PGC2minusCLOZUK_CLOZUKminusPGC2_extended_gene_regions_Clumped_whole_genome_final_significance_threshold_at_1.profile")

fam2 <- fread("~/Dropbox/Stationary_data/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")
covariates <- fread("~/Dropbox/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")

setkey(data_plot,FID)

# Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
#Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]



#PRS.profiles.1 <- data_plot
#PRS.profiles.1 <- data_plot_2_whole_genome_full 
PRS.profiles.1 <- data_plot_2_whole_genome_genic

PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles.1, by.x="FID", by.y="FID", all = F)
PRS.Profiles.with.covariates <- merge(PRS.Profiles.with.covariates, fam2, by.x = "FID", by.y = "FID", all = F)

#  res$model[i]<-sig[i]

# Calculate model including covariates against the polygenic risk score
#model0<-glm(`SCORE_Lek2015_LoFintolerant_90_5e-08`~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
#model0<-glm(SCORE_Lek2015_LoFintolerant_90_0.05~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
#model0<-glm(SCORE_Lek2015_LoFintolerant_90_1~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
model0<-glm(SCORE~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
m1<-mean(residuals(model0))
sd1<-sd(residuals(model0))

# Calculate the Normalised score
#PRS.Profiles.with.covariates$NORMSCORE_lek_fiveeminus8<-(residuals(model0)-m1)/sd1
#PRS.Profiles.with.covariates$NORMSCORE_lek_zeropointzero5<-(residuals(model0)-m1)/sd1
#PRS.Profiles.with.covariates$NORMSCORE_lek_one<-(residuals(model0)-m1)/sd1
PRS.Profiles.with.covariates$NORMSCORE<-(residuals(model0)-m1)/sd1
#             hist(PRS.profiles$NORMSCORE)

# change the Phenotypes so that they will work in a binary model
#PRS.Profiles.with.covariates$PHENO.y <- PRS.Profiles.with.covariates$PHENO - 1
PRS.Profiles.with.covariates$PHENO.y <- PRS.Profiles.with.covariates$PHENO.y - 1
PRS.Profiles.with.covariates$PHENO.y <- as.factor(PRS.Profiles.with.covariates$PHENO.y)


   #ggplot(PRS.Profiles.with.covariates, aes(x=NORMSCORE_lek_fiveeminus8)) + geom_density(aes(group=PHENO.y, colour = PHENO.y, fill = PHENO.y), alpha = 0.3)

   #ggplot(PRS.Profiles.with.covariates, aes(x=NORMSCORE_lek_zeropointzero5)) + geom_density(aes(group=PHENO.y, colour = PHENO.y, fill = PHENO.y), alpha = 0.3)
  # ggplot(PRS.Profiles.with.covariates, aes(x=NORMSCORE_lek_one)) + geom_density(aes(group=PHENO.y, colour = PHENO.y, fill = PHENO.y), alpha = 0.3)
   ggplot(PRS.Profiles.with.covariates, aes(x=NORMSCORE)) + geom_density(aes(group=PHENO.y, colour = PHENO.y, fill = PHENO.y), alpha = 0.3)
   