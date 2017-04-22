# Library
library(mgcv)
library(testit)
library(data.table)
library(parallel)


## set new environment
e <- new.env()

## logistic regression calculation #
assign("sig" , c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5), envir = e)

Calculating_PRS <- function(i) {
  
  # Use the same table and select out the pvalue matching all the genes
  Files_to_parse <- e$Files_to_parse_total[which(e$Files_to_parse_total$pval == e$sig[i]), ]
  
  # Calculate PRS per gene per individual for each pval threshold
  for (l in 1:nrow(Files_to_parse)){
    
    if( l == 1){
      cat("threshold", e$sig[i], "has begun")
    }
    
    # Check if you can read in the file from reference guide 
    #if(has_error(fread(paste0("./output/Profiles/whole_Genome_test_", Files_to_parse$Genes[l], "_", Files_to_parse$pval[l], "_a.profile"))) == T){
      #unrecorded_pvals_genes <- Files_to_parse[l,]
      #erroneous_genes <- rbind(erroneous_genes, unrecorded_pvals_genes)
      #next()
    #}
    
    # read in profiles
    PRS.profiles <-fread(paste0("./output/Profiles/whole_Genome_test_", Files_to_parse$Genes[l], "_", Files_to_parse$pval[l], "_a.profile"))
    
    if (length(which(PRS.profiles$PHENO == -9)) >= 1){
      rows_to_remove <- which(PRS.profiles$PHENO == -9)
      PRS.profiles <- PRS.profiles[!rows_to_remove]
    }
    
    # check if anyone actually has a score
    if(all(PRS.profiles$SCORE==0)){
      no_score_values <- Files_to_parse[l,]
      e$no_score <- rbind(e$no_score, no_score_values)
      next()
    }
    
    PRS.Profiles.with.covariates <- merge(e$covariates,PRS.profiles,by.x="FID", by.y="FID", all = F)
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
    PRS.Profiles.with.covariates$PHENO <- PRS.Profiles.with.covariates$PHENO - 1
    model<-glm(PHENO~PRS.Profiles.with.covariates$NORMSCORE, family=binomial, data=PRS.Profiles.with.covariates)
    
    #             res$Effect[i]<-summary(model)$coefficients[2, "Estimate"]
    #             res$SE[i]<-summary(model)$coefficients[2, "Std. Error"]
    
    # print(c(summary(model)$coefficients[2, "Pr(>|z|)"],Files_to_parse$Genes[l]))
    
    e$residuals[[i]][l,] <- c(summary(model)$coefficients[2, "Pr(>|z|)"],Files_to_parse$Genes[l])
    
    if( l == 10){
      stop()
    }
    #             res$NSNPs[i]<-max(dat$CNT)/2
  }#i
}


## Read in covariates and fam file !may need to put fam file earlier
assign("covariates", fread("./output/CLOZUK2.r7.select2PC.eigenvec.txt"), envir = e)
assign("fam2", fread("./output/CLOZUK.r7.GWAS_IDs.fam"), envir = e)
colnames(e$fam2) <- c("FID","IID","PID","MID","Sex","PHENO")

## Read in the file identifiers and check for duplicates
Files_to_parse_total <- fread("./output/Index_of_genes_and_pval_1.txt",stringsAsFactors = F)
names(Files_to_parse_total) <- c("Genes", "pval")
Files_to_parse_total <- Files_to_parse_total[,unique(Genes), by = pval]
Files_to_parse_total <- Files_to_parse_total[,.(V1,pval)]
names(Files_to_parse_total) <- c("Genes", "pval")

assign("Files_to_parse_total", Files_to_parse_total, envir = e)
duplicate_rows <- uniquecombs(e$Files_to_parse_total, ordered = T)

## Create list for storage of p-values from logistic regression
x <- matrix(data=rep(NA,2*length(unique(e$Files_to_parse_total$Genes))), ncol = 2, nrow = 10407)
assign("residuals", list(x,x,x,x,x,x,x,x,x), envir = e)
names(e$residuals) <- e$sig

## Stop if duplicates exist
if (nrow(duplicate_rows) > nrow(Files_to_parse_total)) {
  stop("duplicate genes found from profile analysis")
}

# Create a matrix recording files which do not exist
assign("erroneous_genes", matrix(ncol = 2), envir = e)
colnames(e$erroneous_genes) <- c("Genes", "pval")

#create a matric recording instances where no score was found
assign("no_score", matrix(ncol=2), envir = e)
colnames(e$no_score) <- c("Genes", "pval")

# memory.limit(20000)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

# Export the environment to the cluster
clusterExport(cl, "e")
parLapply(cl, 1:length(e$sig), Calculating_PRS)
stopCluster(cl)


## Calculate PRS using plink file format

Calculating_PRS <- function(i) {
  
  # Use the same table and select out the pvalue matching all the genes
  Files_to_parse <- Files_to_parse_total[which(Files_to_parse_total$pval == sig[i]), ]
  
  # Calculate PRS per gene per individual for each pval threshold
  for (l in 1:nrow(Files_to_parse)){
    
    if( l == 1){
      cat("threshold",sig[i],"has begun")
    }
    
    # Check if you can read in the file from reference guide 
    if(has_error(fread(paste0("./output/Profiles/whole_Genome_test_", Files_to_parse$Genes[l], "_", Files_to_parse$pval[l], "_a.profile"))) == T){
      unrecorded_pvals_genes <- Files_to_parse[l,]
      erroneous_genes <- rbind(erroneous_genes, unrecorded_pvals_genes)
      next()
    }
    
    # read in profiles
    PRS.profiles <-fread(paste0("./output/Profiles/whole_Genome_test_", Files_to_parse$Genes[l], "_", Files_to_parse$pval[l], "_a.profile"))
     
     if (length(which(PRS.profiles$PHENO == -9)) >= 1){
    rows_to_remove <- which(PRS.profiles$PHENO == -9)
    PRS.profiles <- PRS.profiles[!rows_to_remove]
    }
    # check if anyone actually has a score
    if(all(PRS.profiles$SCORE==0)){
      no_score_values <- Files_to_parse[l,]
      no_score <- rbind(no_score, no_score_values)
      next()
    }
    
    PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles,by.x="FID", by.y="FID", all = F)
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
  PRS.Profiles.with.covariates$PHENO <- PRS.Profiles.with.covariates$PHENO - 1
  model<-glm(PHENO~PRS.Profiles.with.covariates$NORMSCORE, family=binomial, data=PRS.Profiles.with.covariates)
  
  #             res$Effect[i]<-summary(model)$coefficients[2, "Estimate"]
  #             res$SE[i]<-summary(model)$coefficients[2, "Std. Error"]
  residuals[[i]][l,] <- c(summary(model)$coefficients[2, "Pr(>|z|)"],Files_to_parse$Genes[l])
  #             res$NSNPs[i]<-max(dat$CNT)/2
}#i
}
}
for (i in 1:length(sig)){
  write.table(residuals[[i]], file = paste0(sig[i],"CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file.txt"), quote = F, row.names = F)
}
mtext("FT4-Replication. MAF>0.1", side = 1, line = -1, outer = TRUE)
write.table(file=fout, res, col.names=T, row.names=F, quote=F, sep="\t")
