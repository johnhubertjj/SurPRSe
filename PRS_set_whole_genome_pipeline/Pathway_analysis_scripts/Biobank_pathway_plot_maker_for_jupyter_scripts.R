
list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly", "devtools", "gridExtra", "cowplot", "repr", "knitr", "kableExtra", "IRdisplay")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


#library(ggplot2) 

# loads the required packages
lapply(list.of.packages, require, character.only = TRUE)

list_of_datasets <- list()

#21003-2.0 - age at scan
#25010-2.0 Volume brain (White+grey matter)
#25878-2.0 Left Thalamus
#25879-2.0 Right Thalamus
#25880-2.0 Left Caudate
#25881-2.0 Right Caudate
#25882-2.0 Left Putamen
#25883-2.0 Right Putamen
#25884-2.0 Left Pallidum
#25885-2.0 Right Pallidum
#25886-2.0 Left Hippocampus
#25887-2.0 Right Hippocampus
#25888-2.0 Left Amygdala
#25889-2.0 Right Amygdala
#25890-2.0 Left Ventral striatum
#25892-2.0 Right Ventral Striatum
Xaviers_dataset <- fread("/Users/johnhubert/Dropbox/ALSPAC_brain_regions_paper_and_extra_information/Ant_data.csv")
covariates_and_samples <- fread("/Users/johnhubert/Dropbox/ALSPAC_work_Biobank_whole_genome_PRS/Image_UK_Irish_with8cov.txt")

setnames(Xaviers_dataset, old = "eid", new = "app17044")

setnames(covariates_and_samples, old = colnames(covariates_and_samples[,4:190]), new = substring(names(covariates_and_samples[,4:190]), 2))
setnames(Xaviers_dataset, old = colnames(Xaviers_dataset[,2:ncol(Xaviers_dataset)]), new = gsub(x = colnames(Xaviers_dataset[,2:ncol(Xaviers_dataset)]), pattern = "\\-", replacement = "\\."))

vector_to_combine <- colnames(Xaviers_dataset)
new_combined_samples <- merge(x = Xaviers_dataset, y = covariates_and_samples, by = vector_to_combine, all.y = T)



setnames(new_combined_samples, old = "app6553", new = "FID")
setnames(x = new_combined_samples, old = c("21003.2.0","25010.2.0","25878.2.0", "25879.2.0","25880.2.0","25881.2.0","25882.2.0","25883.2.0","25884.2.0","25885.2.0","25886.2.0","25887.2.0","25888.2.0","25889.2.0","25890.2.0", "25892.2.0"),
         new = c("age", "ICV", "Lthal_grey", "Rthal_grey", "Lcaud_grey", "Rcaud_grey", "Lput_grey", "Rput_grey", "Lpal_grey","Rpal_grey","Lhippo_grey","Rhippo_grey",
                 "Lamyg_grey","Ramyg_grey","Lventralstria_grey","Rventralstria_grey"))  

setnames(x = new_combined_samples, old = c("25011.2.0", "25012.2.0", "25013.2.0","25014.2.0","25015.2.0","25016.2.0","25017.2.0","25018.2.0","25019.2.0","25020.2.0","25021.2.0","25022.2.0","25023.2.0","25024.2.0"),
         new = c("Lthal", "Rthal", "Lcaud", "Rcaud", "Lput", "Rput", "Lpal","Rpal","Lhippo","Rhippo",
                 "Lamyg","Ramyg","Laccumbens","Raccumbens"))

columns_to_standardise <- c("ICV","Lthal_grey", "Rthal_grey", "Lcaud_grey", "Rcaud_grey", "Lput_grey", "Rput_grey", "Lpal_grey","Rpal_grey","Lhippo_grey","Rhippo_grey",
                            "Lamyg_grey","Ramyg_grey","Lventralstria_grey","Rventralstria_grey","Lthal", "Rthal", "Lcaud", "Rcaud", "Lput", "Rput", "Lpal","Rpal","Lhippo","Rhippo",
                            "Lamyg","Ramyg","Laccumbens","Raccumbens")

########################################
# adding in arguments from BASH script #
########################################
#args <- commandArgs(trailingOnly = T)
#getwd()
#print(args)

# Specify the different input tables #
#Training_name <- args[3]
#Validation_name <- args [4]
#Scale_phenotypes <- args [5]
#Scale_PRS <- args[6]
#Average_brain_regions <-args[7]
#significance_thresholds <- as.numeric(args[c(8:length(args))])
#print(significance_thresholds)

print("Attempting to automate figure creation...optimism at an all time high...")

#Training_names <- c("HG19_pgc.scz.full.2012-04", "scz.swe.pgc1.results.v3","PGC2","CLOZUK_PGC2noclo")

Training_names <- "CLOZUK_PGC2noclo"
Validation_name <- "Biobank_3000brains"
Scale_phenotypes <- T
Scale_PRS <- F
significance_thresholds <- c(5e-08, 1e-06, 1e-04, 0.01, 0.05, 0.1, 0.2, 0.5, 1)


ALL_PRS_DATA_2 <- read.table("Documents/test_biobank/FINAL_PATHWAY_RESULTS_PRS_PROFILESCLOZUK_PGC2noclo_Biobank_3000brains_Selected_Pocklington_plus_GO_pathways_SCHIZ.txt", header = T)


for (Training_name in Training_names){ 
  #Training_name <- "PGC2"
  #Training_name <- "CLOZUK_PGC2noclo"
  #Training_name <- "scz.swe.pgc1.results.v3"
  #Training_name <- "HG19_pgc.scz.full.2012-04"
  
  for (l in 3:(ncol(ALL_PRS_DATA_2))){
    current_polygenic_score <- ALL_PRS_DATA_2[,c(1,l)]
    test_one <- merge(current_polygenic_score,new_combined_samples,  by= "FID")
    model0 <- glm(test_one[,2]~c01+c02+c03+c04+c05+c06+c08+array, family = gaussian, data = test_one)
    m1<-mean(residuals(model0))
    sd1<-sd(residuals(model0))
    current_polygenic_score$NORMSCORE <- (residuals(model0)-m1)/sd1
    ALL_PRS_DATA_2[,l] <- current_polygenic_score$NORMSCORE
  }
  
  if (Scale_phenotypes == T){
    DATA_pheno_scaled <- cbind(new_combined_samples[,.(FID,gender,age)],(scale(new_combined_samples[,columns_to_standardise,with=FALSE], center = T, scale = T)))
    DATA_pheno <- DATA_pheno_scaled
  }
  
  # Add names of phenotypes corresponding to the table
  subcortical_grey <-colnames(DATA_pheno[,c(5:18)])
  subcortical_grey_left <- subcortical_grey[c(T,F)]
  subcortical_grey_right <- subcortical_grey[c(F,T)]
  subcortical_grey_names <- c("avrg_LR_thal_grey","avrg_LR_caud_grey","avrg_LR_put_grey","avrg_LR_pal_grey","avrg_LR_hippo_grey","avrg_LR_amyg_grey","avrg_LR_ventralstria_grey")
  
  # Add names of phenotypes corresponding to the table
  subcortical <- colnames(DATA_pheno[,c(19:32)])
  subcortical_left <- subcortical[c(T,F)]
  subcortical_right <- subcortical[c(F,T)]
  subcortical_names <- c("avrg_LR_thal","avrg_LR_caud","avrg_LR_put","avrg_LR_pal","avrg_LR_hippo","avrg_LR_amyg","avrg_LR_accumbens")
  
  for (i in 1:length(subcortical_left)){
    names <- subcortical_names[i]
    new_df <- data.frame(DATA_pheno[[subcortical_left[i]]],DATA_pheno[[subcortical_right[i]]])
    current_column <- rowMeans(new_df,na.rm = T)
    DATA_pheno[[names]] <- current_column
  } 
  
  if (Scale_PRS == T){
    ALL_PRS_DATA_2[-1:-2] <- scale(ALL_PRS_DATA_2[,-1:-2], center = T, scale = T)
  }
  
  
  DATA <- join_all(list(ALL_PRS_DATA_2, DATA_pheno), by=c("FID"), type='inner')
  
  # create a score list
  score_list <- colnames(ALL_PRS_DATA_2[,-1:-2])
  phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])
  DATA$gender <- as.factor(DATA$gender)
  
  # Remove individuals with missing values
  #remove_ICV <- which(is.na(DATA$ICV))
  #DATA <- DATA[-remove_ICV,]
  
  # dont' change obj it makes the loop output a dataframe with the regression results
  obj <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
  obj_ICV <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
  
  # this example if for a linear regression (the phenotype of interest is a quantiative trait)
  # is using a discrete phenotype, a logistic regression needs to be run, and the code altered from 'lm' to 'glm' including the argument of 'family = binomial'
  # alterations for the calculation of R2 will also need to be made using the command highlighted above
  
  for (i in score_list) {
    for (j in phenotypes) {
      fit <- lm(DATA[,j] ~ DATA[,i] + gender + age + ICV,  data=DATA)
      fit1 <- lm(DATA[,j] ~ gender + age + ICV, data=DATA)
      tmp <- coef(summary(fit))
      tmp2 <- summary(fit)
      hold <- summary(fit1)
      true_r2 <- tmp2$r.squared - hold$r.squared
      CI <- confint(fit, level=0.95)
      CI <- CI[2,]
      tmp3 <- c(j,i,tmp[2,], true_r2, CI)
      obj <- rbind(obj, tmp3)
    }
  }
  
  for (i in score_list){
    fit <- lm(ICV ~ DATA[,i] + gender + age,  data=DATA)
    fit1 <- lm(ICV ~ gender + age, data=DATA)
    tmp <- coef(summary(fit))
    tmp2 <- summary(fit)
    hold <- summary(fit1)
    true_r2 <- tmp2$r.squared - hold$r.squared
    CI <- confint(fit, level=0.95)
    CI <- CI[2,]
    tmp3 <- c("ICV",i,tmp[2,], true_r2, CI)
    obj_ICV <- rbind(obj_ICV, tmp3)
  }
  
  # this is a clean-up step - do not change
  results <- obj[which(obj$score %in% score_list),]
  results2 <- obj_ICV[which(obj_ICV$score %in% score_list),]
  
  assign(paste0("results_",Training_name), results)
  assign(paste0("results_ICV",Training_name), results2)
}

Training_names <-  "CLOZUK_PGC2noclo"
Validation_name <- "Biobank_3000brains"
Scale_phenotypes <- T
Scale_PRS <- F #The PRS for biobank are already scaled in the for l loop using covariates
significance_thresholds <- c(5e-08, 1e-06, 1e-04, 0.01, 0.05, 0.1, 0.2, 0.5, 1)




for (Training_name in Training_names){ 
  
  my_files <- paste0("Documents/test_biobank/", Training_name, "_", Validation_name, "_normal_gene_regions_Clumped_whole_genome_final_significance_threshold_at_", significance_thresholds, ".profile")
  
  my_data <- lapply(my_files, read.table, header=TRUE) 
  names(my_data) <- str_replace(my_files, pattern = ".profile", replacement = "")
  
  # iterate through significance thresholds 
  for (i in 1:length(significance_thresholds)) {
    my_data[[i]] <- my_data[[i]][,c(1,2,6)]
    colnames(my_data[[i]]) <- c("FID", "IID", paste("SCORE_", significance_thresholds[i], sep=""))
  }
  
  all_prs <- join_all(my_data, by=c("FID", "IID"), type = 'left')
  
  
  for (l in 3:(length(significance_thresholds) + 2)){
    current_polygenic_score <- all_prs[,c(1,l)]
    test_one <- merge(current_polygenic_score,new_combined_samples,  by= "FID")
    model0 <- glm(test_one[,2]~c01+c02+c03+c04+c05+c06+c08+array, family = gaussian, data = test_one)
    m1<-mean(residuals(model0))
    sd1<-sd(residuals(model0))
    current_polygenic_score$NORMSCORE <- (residuals(model0)-m1)/sd1
    all_prs[,l] <- current_polygenic_score$NORMSCORE
  }
  
  if (Scale_phenotypes == T){
    DATA_pheno_scaled <- cbind(new_combined_samples[,.(FID,gender,age)],(scale(new_combined_samples[,columns_to_standardise,with=FALSE], center = T, scale = T)))
    DATA_pheno <- DATA_pheno_scaled
  }
  
  # Add names of phenotypes corresponding to the table
  subcortical_grey <-colnames(DATA_pheno[,c(5:18)])
  subcortical_grey_left <- subcortical[c(T,F)]
  subcortical_grey_right <- subcortical[c(F,T)]
  subcortical_grey_names <- c("avrg_LR_thal_grey","avrg_LR_caud_grey","avrg_LR_put_grey","avrg_LR_pal_grey","avrg_LR_hippo_grey","avrg_LR_amyg_grey","avrg_LR_ventralstria_grey")
  
  # Add names of phenotypes corresponding to the table
  subcortical <- colnames(DATA_pheno[,c(19:32)])
  subcortical_left <- subcortical[c(T,F)]
  subcortical_right <- subcortical[c(F,T)]
  subcortical_names <- c("avrg_LR_thal","avrg_LR_caud","avrg_LR_put","avrg_LR_pal","avrg_LR_hippo","avrg_LR_amyg","avrg_LR_accumbens")
  
  for (i in 1:length(subcortical_left)){
    names <- subcortical_names[i]
    new_df <- data.frame(DATA_pheno[[subcortical_left[i]]],DATA_pheno[[subcortical_right[i]]])
    current_column <- rowMeans(new_df,na.rm = T)
    DATA_pheno[[names]] <- current_column
  } 
  
  DATA <- join_all(list(all_prs, DATA_pheno), by=c("FID"), type='inner')
  
  # create a score list
  score_list <- colnames(all_prs[,-1:-2])
  phenotypes <- colnames(DATA_pheno[,19:39])
  
  #Training_names <- c("HG19_pgc.scz.full.2012-04", "scz.swe.pgc1.results.v3","PGC2","CLOZUK_PGC2noclo")
  # dont' change obj it makes the loop output a dataframe with the regression results
  obj <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
  obj_ICV <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
  
  # this example if for a linear regression (the phenotype of interest is a quantiative trait)
  # is using a discrete phenotype, a logistic regression needs to be run, and the code altered from 'lm' to 'glm' including the argument of 'family = binomial'
  # alterations for the calculation of R2 will also need to be made using the command highlighted above
  
  for (i in score_list) {
    for (j in phenotypes) {
      fit <- lm(DATA[,j] ~ DATA[,i] + gender + age + ICV,  data=DATA)
      fit1 <- lm(DATA[,j] ~ gender + age + ICV, data=DATA)
      tmp <- coef(summary(fit))
      tmp2 <- summary(fit)
      hold <- summary(fit1)
      true_r2 <- tmp2$r.squared - hold$r.squared
      CI <- confint(fit, level=0.95)
      CI <- CI[2,]
      tmp3 <- c(j,i,tmp[2,], true_r2, CI)
      obj <- rbind(obj, tmp3)
    }
  }
  
  for (i in score_list){
    fit <- lm(ICV ~ DATA[,i] + gender + age,  data=DATA)
    fit1 <- lm(ICV ~ gender + age, data=DATA)
    tmp <- coef(summary(fit))
    tmp2 <- summary(fit)
    hold <- summary(fit1)
    true_r2 <- tmp2$r.squared - hold$r.squared
    CI <- confint(fit, level=0.95)
    CI <- CI[2,]
    tmp3 <- c("ICV",i,tmp[2,], true_r2, CI)
    obj_ICV <- rbind(obj_ICV, tmp3)
  }
  
  # this is a clean-up step - do not change
  results <- obj[which(obj$score %in% score_list),]
  results2 <- obj_ICV[which(obj_ICV$score %in% score_list),]
  

  assign(paste0("results_",Validation_name), results)
  assign(paste0("results_ICV",Validation_name), results2)
  
  phenotypes <- colnames(DATA_pheno[,5:18])
  
  #Training_names <- c("HG19_pgc.scz.full.2012-04", "scz.swe.pgc1.results.v3","PGC2","CLOZUK_PGC2noclo")
  # dont' change obj it makes the loop output a dataframe with the regression results
  obj <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
  obj_ICV <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
  
  # this example if for a linear regression (the phenotype of interest is a quantiative trait)
  # is using a discrete phenotype, a logistic regression needs to be run, and the code altered from 'lm' to 'glm' including the argument of 'family = binomial'
  # alterations for the calculation of R2 will also need to be made using the command highlighted above
  
  for (i in score_list) {
    for (j in phenotypes) {
      fit <- lm(DATA[,j] ~ DATA[,i] + gender + age + ICV,  data=DATA)
      fit1 <- lm(DATA[,j] ~ gender + age + ICV, data=DATA)
      tmp <- coef(summary(fit))
      tmp2 <- summary(fit)
      hold <- summary(fit1)
      true_r2 <- tmp2$r.squared - hold$r.squared
      CI <- confint(fit, level=0.95)
      CI <- CI[2,]
      tmp3 <- c(j,i,tmp[2,], true_r2, CI)
      obj <- rbind(obj, tmp3)
    }
  }
  
  
  # this is a clean-up step - do not change
  results <- obj[which(obj$score %in% score_list),]
  
  assign(paste0("results_",Training_name,"_grey_volume"), results)
}

ALSPAC_all_pathways <- as.data.table(results_CLOZUK_PGC2noclo)


ALSPAC_full_brain_regions <- results_Biobank_3000brains

testa <- ALSPAC_full_brain_regions
ALSPAC_full_brain_regions <- as.data.table(ALSPAC_full_brain_regions)

# Pathways
str(ALSPAC_all_pathways)

ALSPAC_only_0.5 <- ALSPAC_all_pathways[grepl("_0.5",ALSPAC_all_pathways$score)] 
ALSPAC_only_0.05 <- ALSPAC_all_pathways[grepl("_0.05",ALSPAC_all_pathways$score)] 

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("abnormal_grooming_behavior", ALSPAC_only_0.05$score)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("abnormal_grooming_behavior", ALSPAC_only_0.5$score)]

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("Lek2015_LoFintolerant_90", ALSPAC_only_0.05$score)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("Lek2015_LoFintolerant_90", ALSPAC_only_0.5$score)]

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("grey", ALSPAC_only_0.05$test)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("grey", ALSPAC_only_0.5$test)]

split_data_frame_0.05 <- split(ALSPAC_only_0.05, by = "test")
split_data_frame_0.5 <- split(ALSPAC_only_0.5, by = "test")
#setkey(split_data_frame_0.05,p)
#setkey(split_data_frame_0.5,p)

#Whole brain
ALSPAC_all_brain_0.05 <- ALSPAC_full_brain_regions[grepl("_0.05",ALSPAC_full_brain_regions$score)]
ALSPAC_all_brain_0.5 <- ALSPAC_full_brain_regions[grepl("_0.5",ALSPAC_full_brain_regions$score)] 
gglist_pathways <- list()
gglist_beta_pathways <- list()
gglist_r2_dir_pathways <- list()

for (i in 1:length(split_data_frame_0.05)){
  
  # Obtaining the log p-values
  top_five_pathways <- split_data_frame_0.05[[i]][order(split_data_frame_0.05[[i]]$p)]
  top_five_pathways <- top_five_pathways[1:5]
  whole_brain_p_value_0.05 <- ALSPAC_all_brain_0.05[which(top_five_pathways$test[1] == ALSPAC_all_brain_0.05)]
  top_five_pathways <- rbind(top_five_pathways,whole_brain_p_value_0.05)
  top_five_pathways$P_threshold <- 0.05
  top_five_pathways$Type <- c(rep("Pathway",5),"Whole_genome")
  top_five_pathways_0.05 <- which(top_five_pathways$score == "SCORE_0.05")
  
  # By Factoring the p-values, you can get them to trend decreasingly without affecting the whole genome p-value
  top_five_pathways$p <- as.numeric(top_five_pathways$p)
  top_five_pathways$score <- factor(top_five_pathways$score, levels = top_five_pathways$score[c(top_five_pathways_0.05, order(top_five_pathways$p[-top_five_pathways_0.05]))])
  
  
  # Do the same for the 0.5 p-thresholds
  top_five_pathways2 <- split_data_frame_0.5[[i]][order(split_data_frame_0.5[[i]]$p)]
  top_five_pathways2 <- top_five_pathways2[1:5]
  whole_brain_p_value_0.5 <- ALSPAC_all_brain_0.5[which(top_five_pathways2$test[1] == ALSPAC_all_brain_0.5)]
  top_five_pathways2 <- rbind(top_five_pathways2,whole_brain_p_value_0.5)
  top_five_pathways2$P_threshold <- 0.5
  top_five_pathways2$Type <- c(rep("Pathway",5),"Whole_genome")
  
  top_five_pathways_0.5 <- which(top_five_pathways2$score == "SCORE_0.5")
  top_five_pathways2$p <- as.numeric(top_five_pathways2$p)
  top_five_pathways2$score <- factor(top_five_pathways2$score, levels = top_five_pathways2$score[c(top_five_pathways_0.5, order(top_five_pathways2$p[-top_five_pathways_0.5]))])
  
  
  top_five_pathways_all <- rbind(top_five_pathways,top_five_pathways2)
  top_five_pathways_all$p <- as.numeric(top_five_pathways_all$p)
  top_five_pathways_all$logp <- -log10(top_five_pathways_all$p)
  
  top_five_pathways_all$estimate <- as.numeric(top_five_pathways_all$estimate)
  top_five_pathways_all$upper <- as.numeric(top_five_pathways_all$upper)
  top_five_pathways_all$lower <- as.numeric(top_five_pathways_all$lower)
  top_five_pathways_all$SE <- as.numeric(top_five_pathways_all$SE)
  
  top_five_pathways_all$SE_higher <- top_five_pathways_all$estimate + top_five_pathways_all$SE
  top_five_pathways_all$SE_lower <- top_five_pathways_all$estimate - top_five_pathways_all$SE
  
  top_five_pathways_all$r2_dir <- 100 * (as.numeric(top_five_pathways_all$r.squared) *
                                           (sign(as.numeric(top_five_pathways_all$estimate))))
  
  top_five_pathways_all$p_value_text <- paste("p =", scientific(top_five_pathways_all$p, digits = 2), sep = " ")
  assign(paste0(top_five_pathways_all$test[1],"_top_five_pathways_all_thresholds"), top_five_pathways_all, envir = .GlobalEnv)
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=logp, fill = Type, group=P_threshold))
  
  p <- p +
    geom_point(aes(colour = Type))
  
  p <- p + scale_x_discrete(labels=c("SCORE_0.05" = "Whole Genome PRS", "SCORE_0.5" = "Whole Genome PRS"))
  #p <- p + scale_y
  p <- p + facet_grid(. ~ P_threshold,scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = expression(-'log'[10]*'(p)'))
  p <- p + xlab(label = "Polygenic risk score")
  p
  gglist_pathways[[i]] <- p
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=estimate, fill = Type, group=P_threshold))
  
  p <- p +
    geom_errorbar(aes(ymin = SE_higher, ymax = SE_lower), position = "dodge", width = 0.25) +
    geom_point(aes(colour = Type))
  
  p <- p + scale_x_discrete(labels=c("SCORE_0.05" = "Whole Genome PRS", "SCORE_0.5" = "Whole Genome PRS"))
  p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "BETA")
  p <- p + xlab(label = "Polygenic risk score")
  
  gglist_beta_pathways[[i]] <- p
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=r2_dir, fill = Type, group=P_threshold))
  p <- p +
    geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
    geom_text(data=subset(top_five_pathways_all, p < 0.05),
              aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
  
  p <- p + scale_x_discrete(labels=c("SCORE_0.05" = "Whole Genome PRS", "SCORE_0.5" = "Whole Genome PRS"))
  p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
  p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "R2_dir (%)")
  p <- p + xlab(label = "Polygenic risk score")
  
  gglist_r2_dir_pathways[[i]] <- p
  
}

pdf("~/Documents/Biobank_CLOZUK_pathway_whole_genome_genic_comparison_R2_dir_9_pathways.pdf", onefile = TRUE)
invisible(lapply(gglist_r2_dir_pathways,print))
dev.off()
