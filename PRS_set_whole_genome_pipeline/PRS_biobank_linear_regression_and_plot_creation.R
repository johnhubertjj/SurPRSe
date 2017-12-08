# Biobank linear regression #
library(data.table)

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

setnames(x = covariates_and_samples, old = c("X21003.2.0","X25010.2.0","X25878.2.0", "X25879.2.0","X25880.2.0","X25881.2.0","X25882.2.0","X25883.2.0","X25884.2.0","X25885.2.0","X25886.2.0","X25887.2.0","X25888.2.0","X25889.2.0","X25890.2.0", "X25892.2.0"),
                                     new = c("age", "ICV", "Lthal", "Rthal", "Lcaud", "Rcaud", "Lput", "Rput", "Lpal","Rpal","Lhippo","Rhippo",
                                             "Lamyg","Ramyg","Lventralstria","Rventralstria"))  

columns_to_standardise <- c("Lthal", "Rthal", "Lcaud", "Rcaud", "Lput", "Rput", "Lpal","Rpal","Lhippo","Rhippo",
"Lamyg","Ramyg","Lventralstria","Rventralstria")

covariates_and_samples <- fread("/home/johnhubert/Dropbox/ALSPAC_work_Biobank_whole_genome_PRS/Image_UK_Irish_with8cov.txt")
setnames(covariates_and_samples, old = "app6553", new = "FID")

Training_names <- c("PGC1", "PGC1_sweden", "PGC2noCLOZUK", "CLOZUK_PGC2noclo")
Validation_name <- "Biobank_3000brains"
Scale_phenotypes <- T
Scale_PRS <- F
significance_thresholds <- c(5e-08, 1e-06, 1e-04, 0.01, 0.05, 0.1, 0.2, 0.5, 1)

library(ggplot2)
library(gridExtra)
library(data.table)
library(cowplot)

list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# loads the required packages
lapply(list.of.packages, require, character.only = TRUE)

for (Training_name in Training_names){ 

  
  my_files <- paste0("/home/johnhubert/Dropbox/ALSPAC_work_Biobank_whole_genome_PRS/PRS_scoring/", Training_name, "_", Validation_name, "_whole_genome_significance_threshold_at_", significance_thresholds, ".profile")
  
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
   test_one <- merge(current_polygenic_score,covariates_and_samples,  by= "FID")
   model0 <- glm(test_one[,2]~c01+c02+c03+c04+c05+c06+c08+array, family = gaussian, data = test_one)
   m1<-mean(residuals(model0))
   sd1<-sd(residuals(model0))
   current_polygenic_score$NORMSCORE <- (residuals(model0)-m1)/sd1
   all_prs[,l] <- current_polygenic_score$NORMSCORE
  }
  
  if (Scale_phenotypes == T){
    DATA_pheno_scaled <- cbind(covariates_and_samples[,.(FID,gender,age,ICV)],(scale(covariates_and_samples[,columns_to_standardise,with=FALSE], center = T, scale = T)))
    DATA_pheno <- DATA_pheno_scaled
  }
  
  # Add names of phenotypes corresponding to the table
  subcortical <-colnames(DATA_pheno[,c(5:18)])
  subcortical_left <- subcortical[c(T,F)]
  subcortical_right <- subcortical[c(F,T)]
  subcortical_names <- c("avrg_LR_thal","avrg_LR_caud","avrg_LR_put","avrg_LR_pal","avrg_LR_hippo","avrg_LR_amyg","avrg_LR_ventralstria")
  
  for (i in 1:length(subcortical_left)){
    names <- subcortical_names[i]
    new_df <- data.frame(DATA_pheno[[subcortical_left[i]]],DATA_pheno[[subcortical_right[i]]])
    current_column <- rowMeans(new_df,na.rm = T)
    DATA_pheno[[names]] <- current_column
  } 
  
  
  
  DATA <- join_all(list(all_prs, DATA_pheno), by=c("FID"), type='inner')
  
  # create a score list
  score_list <- colnames(all_prs[,-1:-2])
  phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])
  
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
  
  assign(paste0("results_",Training_name), results)
  assign(paste0("results_ICV",Training_name), results2)
}



for (Training_name in Training_names){ 
  #Training_name <- "PGC2noCLOZUK"
  #Training_name <- "CLOZUK_PGC2noclo"
  #Training_name <- "scz.swe.pgc1.results.v3"
  #Training_name <- "HG19_pgc.scz.full.2012-04"
  
  # Write in the phenotype files
  Pheno_brain_thickness <- fread("~/Dropbox/ALSPAC_thickness_27january2017.txt")
  Pheno_brain_volume <-fread("~/Dropbox/ALSPAC_volumes_27january2017.txt") 
  
  # Collate all the scores together into one table_normal:
  
  #my_files <- paste0(Training_name, "_", Validation_name,"_output/PRS_scoring/", Training_name, "_", Validation_name, "_normal_gene_regions_Clumped_whole_genome_final_significance_threshold_at_", significance_thresholds, ".profile")
  my_files <- paste0(Training_name, "_", Validation_name,"_output/PRS_scoring/", Training_name, "_", Validation_name, "_extended_gene_regions_Clumped_whole_genome_final_significance_threshold_at_", significance_thresholds, ".profile")
  
  my_data <- lapply(my_files, read.table, header=TRUE) 
  names(my_data) <- str_replace(my_files, pattern = ".profile", replacement = "")
  
  # iterate through significance thresholds 
  for (i in 1:length(significance_thresholds)) {
    my_data[[i]] <- my_data[[i]][,c(1,2,6)]
    colnames(my_data[[i]]) <- c("FID", "IID", paste("SCORE_", significance_thresholds[i], sep=""))
  }
  
  all_prs <- join_all(my_data, by=c("FID", "IID"), type='left')
  
  if (Scale_PRS == T){
    all_prs[,-1:-2] <- scale(all_prs[,-1:-2], center = T, scale = T)
  }
  
  
  # Join all the tables together
  setnames(Pheno_brain_volume,old = "SubjID", new = "FID")
  setnames(Pheno_brain_thickness,old = "SubjID", new = "FID")
  
  DATA_pheno <- join_all(list(Pheno_brain_thickness, Pheno_brain_volume), by=c("FID", "kz021", "age", "ICV"), type='inner')
  
  if (Scale_phenotypes == T){
    DATA_pheno_scaled <- cbind(DATA_pheno[,1:3],(scale(DATA_pheno[,4:ncol(DATA_pheno)], center = T, scale = T)))
    DATA_pheno <- DATA_pheno_scaled
  }
  
  # Add names of phenotypes corresponding to the table
  phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])
  subcortical <-colnames(DATA_pheno[,c(25:40)])
  subcortical_left <- subcortical[c(T,F)]
  subcortical_right <- subcortical[c(F,T)]
  subcortical_names <- c("avrg_LR_LatVent","avrg_LR_thal","avrg_LR_caud","avrg_LR_put","avrg_LR_pal","avrg_LR_hippo","avrg_LR_amyg","avrg_LR_accumb")
  
  for (i in 1:length(subcortical_left)){
    names <- subcortical_names[i]
    new_df <- data.frame(DATA_pheno[[subcortical_left[i]]],DATA_pheno[[subcortical_right[i]]])
    current_column <- rowMeans(new_df,na.rm = T)
    DATA_pheno[[names]] <- current_column
  } 
  
  
  
  DATA <- join_all(list(all_prs, DATA_pheno), by=c("FID"), type='inner')
  
  # create a score list
  score_list <- colnames(all_prs[,-1:-2])
  phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])
  DATA$kz021 <- as.factor(DATA$kz021)
  
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
      fit <- lm(DATA[,j] ~ DATA[,i] + kz021 + age + ICV,  data=DATA)
      fit1 <- lm(DATA[,j] ~ kz021 + age + ICV, data=DATA)
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
    fit <- lm(ICV ~ DATA[,i] + kz021 + age,  data=DATA)
    fit1 <- lm(ICV ~ kz021 + age, data=DATA)
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

brain_volume_progression <- list(results_PGC1, results_PGC1_sweden, results_PGC2noCLOZUK, results_CLOZUK_PGC2noclo)
ICV_progression <- list(`results_ICVHG19_pgc.scz.full.2012-04`, results_ICVscz.swe.pgc1.results.v3, results_ICVPGC2, results_ICVCLOZUK_PGC2noclo)

#results_scaled_1 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_PGC1_no_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results_scaled_2 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_PGC1_with_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results")
#results_scaled_3 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_PGC2_daner_file_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results_scaled_4 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_CLOZUK_PGC2noclo_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")

#brain_volume_progression_scaled <- list(results_scaled_1, results_scaled_2, results_scaled_3, results_scaled_4)


## PLOTS OF BETA AND SE (AND CI) ## 
newdataframe_beta <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$estimate, brain_volume_progression[[2]]$estimate, brain_volume_progression[[3]]$estimate, brain_volume_progression[[4]]$estimate)
names(newdataframe_beta) <- c("test", "PGC1", "PGC1swe", "PGC2", "CLOZUK_PGC2noclo")
test_df <- t(newdataframe_beta)
colnames(test_df) <- as.character(test_df[1,])
test_df_beta <- test_df[-1,]
test_df_beta <- cbind(test_df_beta,c("PGC1", "PGC1swe", "PGC2", "CLOZUK_PGC2noclo"))
test_df_beta <- as.data.frame(test_df_beta)



number_of_phenotypes <- length(phenotypes)
number_of_thresholds <- length(significance_thresholds)

newdataframe_SE <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$SE, brain_volume_progression[[2]]$SE, brain_volume_progression[[3]]$SE, brain_volume_progression[[4]]$SE)
names(newdataframe_SE) <- c("test", "PGC1", "PGC1swe", "PGC2", "CLOZUK_PGC2noclo")
test_df_SE <- t(newdataframe_SE)
colnames(test_df_SE) <- as.character(test_df_SE[1,])
test_df_SE <- test_df_SE[-1,]
test_df_SE <- cbind(test_df_SE,c("PGC1", "PGC1swe", "PGC2", "CLOZUK_PGC2noclo"))
test_df_SE <- as.data.frame(test_df_SE)

# create the first dataset for one significance threshold
if(number_of_thresholds > 1){
  BETA_dataframe_current <- newdataframe_beta[1:number_of_phenotypes,]
  test_beta <- melt(BETA_dataframe_current, id.vars="test")
  names(test_beta) <- c("Brain_region", "Dataset","BETA")
  
  SE_dataframe_current <- newdataframe_SE[1:number_of_phenotypes,]
  test_SE <- melt(SE_dataframe_current, id.vars="test")
  names(test_SE) <- c("Brain_region", "Dataset","SE")
  test_SE$SE <- as.numeric(test_SE$SE)
  test_beta$BETA <- as.numeric(test_beta$BETA)
  Beta_se_plots <- cbind(test_beta,test_SE$SE)
  Beta_se_plots$upper <- Beta_se_plots$BETA + Beta_se_plots$`test_SE$SE`
  Beta_se_plots$lower <- Beta_se_plots$BETA - Beta_se_plots$`test_SE$SE`
  number_of_datasets <- nrow(test_df_beta)
  Beta_se_plots$P_Value_threshold <- c(rep(significance_thresholds[1],number_of_phenotypes*number_of_datasets))  
  names(Beta_se_plots) <- c("Brain_region", "Dataset","BETA", "SE", "upper", "lower", "Pvaluethreshold")
  all_plots <- Beta_se_plots
  
  # Loop through the rest based on the number of significance thresholds you have
  
  
  for (i in 2:number_of_thresholds){
    # This looks complicated but all i am essentially saying is this:
    # for my current dataframe, find out the number of phenotypes (in this case 44)
    # multiply that by the significance threshold we are currently on (because we are going to have 44 rows of that significance threshold) eg 44 * 1  + 1(because we have 2 signifance thresholds)
    # That row (row 45) is the row which starts with a new significance threshold of 0.5
    # then take all the rows up to the last phenotype (44*2) = 88 (in this case our last row)
    # now do the same as above on this subsection of data...
    indexone <- (((number_of_phenotypes)*(i-1)) + 1)
    indextwo <- (number_of_phenotypes*i)
    BETA_dataframe_current <- newdataframe_beta[indexone:indextwo,]
    test_beta_to_add <- melt(BETA_dataframe_current, id.vars="test")
    names(test_beta_to_add) <- c("Brain_region", "Dataset","BETA")
    
    SE_dataframe_current <- newdataframe_SE[indexone:indextwo,]
    test_SE_to_add <- melt(SE_dataframe_current, id.vars="test")
    names(test_SE_to_add) <- c("Brain_region", "Dataset","SE")
    test_SE_to_add$SE <- as.numeric(test_SE_to_add$SE)
    test_beta_to_add$BETA <- as.numeric(test_beta_to_add$BETA)
    Beta_se_plots_to_add <- cbind(test_beta_to_add,test_SE_to_add$SE)
    
    Beta_se_plots_to_add$upper <- Beta_se_plots_to_add$BETA + Beta_se_plots_to_add$`test_SE_to_add$SE`
    Beta_se_plots_to_add$lower <- Beta_se_plots_to_add$BETA - Beta_se_plots_to_add$`test_SE_to_add$SE`
    # repeat the same significance threshold over and over again
    Beta_se_plots_to_add$P_Value_threshold <- c(rep(significance_thresholds[i],number_of_phenotypes*number_of_datasets))  
    names(Beta_se_plots_to_add) <- c("Brain_region", "Dataset","BETA", "SE", "upper", "lower", "Pvaluethreshold")
    
    all_plots <- rbind(all_plots,Beta_se_plots_to_add)
  }
}

all_plots$Pvaluethreshold <- as.factor(all_plots$Pvaluethreshold)
gglist_0.05 <- list() 

all_plots2<- all_plots[all_plots$Pvaluethreshold %in% c(5e-08, 1e-06, 0.05, 0.5),]
all_plots2 <- all_plots2[all_plots2$Brain_region %in% c(subcortical_left,subcortical_right),] 

subcortical_lR <- c(subcortical_left,subcortical_right)
subcortical_phenotypes <- c("Thalamus", "Caudate Nucleus", " Putamen", "Pallidum", "Hippocampus", "Amygdala", "Nucleus Accumbens")

gglist_final <- list()

## Scatter plots with SE error bars in paper format ##
for (i in 1:length(subcortical_left)){
  # Plots for Betas and SE based on files #
  all_plots_current <- all_plots2
  all_plots_current$Brain_region <- as.character(all_plots_current$Brain_region)
  
  all_plots_current <- as.data.table(all_plots_current)
  setkey(all_plots_current,"Brain_region")
  
  all_plots_current1 <- all_plots_current[subcortical_left[i]]
  all_plots_current2 <- all_plots_current[subcortical_right[i]]
  all_plots_currentLR <- rbind(all_plots_current1, all_plots_current2)
  
  p <- ggplot(all_plots_currentLR, aes(x=Dataset, y=BETA, fill = Pvaluethreshold, group=Dataset))
  
  p <- p +
    geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
    geom_point()
  
  p <- p + scale_x_discrete(labels=c("PGC1" = "PGC1", "PGC1swe" = "PGC1swe", "PGC2" = "PGC2", "CLOZUK_PGC2noclo" = "CLOZUK"))
  p <- p + facet_grid(Brain_region ~ Pvaluethreshold) +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_discrete(guide=FALSE)
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(subcortical_phenotypes[i])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  
  gglist_final[[i]] <- p
}

Final_diagram <- plot_grid(gglist_final[[1]], gglist_final[[2]],gglist_final[[3]], gglist_final[[4]], gglist_final[[5]], gglist_final[[6]], gglist_final[[7]], ncol = 4, nrow = 2, labels = "AUTO")
