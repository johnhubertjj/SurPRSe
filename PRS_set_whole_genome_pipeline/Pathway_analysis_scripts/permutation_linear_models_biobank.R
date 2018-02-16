# Linear regression for pathways #
### Start Timer
ptm <- proc.time()

#Print R version
getRversion()    

# devtools::install_github('hadley/ggplot2')

library(ggplot2)
library(gridExtra)
library(data.table)

# required packages -- will only install if they are not already installed
list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# loads the required packages
lapply(list.of.packages, require, character.only = TRUE)

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
Scale_phenotypes <- args [5]
Scale_PRS <- args[6]
Average_brain_regions <-args[7]
significance_thresholds <- as.numeric(args[c(8:length(args))])
print(significance_thresholds)

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


  
# Training_names <- c("HG19_pgc.scz.full.2012-04", "scz.swe.pgc1.results.v3","PGC2","CLOZUK_PGC2noclo")

Training_names <- "CLOZUK_PGC2noclo"
Validation_name <- "Biobank_3000brains"
Scale_phenotypes <- T
Scale_PRS <- F
significance_thresholds <- c(0.05, 0.5)


set_names <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268")

next_name <- "Calcium_ion_import_GO0070509"

for (next_name in set_names){
  
  
  ALL_PRS_DATA_2 <- fread(paste0("./FINAL_PATHWAY_RESULTS_PRS_PROFILESCLOZUK_PGC2noclo_Biobank_3000brains_",next_name,"randomised_gene_sets.txt"),header = T)
  
  for (Training_name in Training_names){ 

    
    for (l in 3:(ncol(ALL_PRS_DATA_2))){
      current_polygenic_score <- ALL_PRS_DATA_2[,c(1,l), with=F]
      test_one <- merge(current_polygenic_score,new_combined_samples,  by= "FID")
      model0 <- glm(test_one[[2]]~c01+c02+c03+c04+c05+c06+c08+array, family = gaussian, data = test_one)
      m1<-mean(residuals(model0))
      sd1<-sd(residuals(model0))
      current_polygenic_score$NORMSCORE <- (residuals(model0)-m1)/sd1
      ALL_PRS_DATA_2[[l]] <- current_polygenic_score$NORMSCORE
    }
    
    if (Scale_phenotypes == T){
      DATA_pheno_scaled <- cbind(new_combined_samples[,.(FID,gender,age)],(scale(new_combined_samples[,columns_to_standardise,with=FALSE], center = T, scale = T)))
      DATA_pheno <- DATA_pheno_scaled
    }
    
    # Add names of phenotypes corresponding to the table
    subcortical_grey <- colnames(DATA_pheno[,c(5:18)])
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
    

    DATA <- join_all(list(ALL_PRS_DATA_2, DATA_pheno), by=c("FID"), type='inner')
    
    # create a score list
    score_list <- colnames(ALL_PRS_DATA_2[,c(3:ncol(ALL_PRS_DATA_2)),with=F])
    phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])
    DATA$gender <- as.factor(DATA$gender)
    
    
    # Remove individuals with missing values
    #remove_ICV <- which(is.na(DATA$ICV))
    #DATA <- DATA[-remove_ICV,]
    # dont' change obj it makes the loop output a dataframe with the regression results
    
    phenotype_pvalue_list <- matrix(NA,nrow = 2000,ncol = length(phenotypes))
    #obj_ICV <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0, lower = 0, upper = 0)
    
    # this example if for a linear regression (the phenotype of interest is a quantiative trait)
    # is using a discrete phenotype, a logistic regression needs to be run, and the code altered from 'lm' to 'glm' including the argument of 'family = binomial'
    # alterations for the calculation of R2 will also need to be made using the command highlighted above
    
    for (j in phenotypes) {
      obj <- NULL
      for (i in score_list) {
        
        fit <- lm(DATA[[j]] ~ DATA[[i]] + gender + age + ICV, data=DATA)
        tmp <- coef(summary(fit))
        tmp <- tmp[2,4]
        obj <- c(obj, tmp)
      }
      phenotype_pvalue_list[,which(phenotypes == j)]<- obj
    }
    colnames(phenotype_pvalue_list) <- phenotypes
    phenotype_pvalue_list <- cbind(phenotype_pvalue_list, score_list)
    write.table(phenotype_pvalue_list,file = paste0(next_name,"_table_of_random_p_values.txt"), col.names = T,row.names = F ,quote=F)
  }
}

real_p_values <- fread("CLOZUK_PGC2noclo_Biobank_Selected_SCZ_pathways.txt")
real_p_values <- real_p_values[,.(test,score,p)]
setkey(real_p_values,score)


for (sig in significance_thresholds){
  
  final_results <- matrix(NA,nrow = 35,ncol=9)
  rownames(final_results) <- phenotypes
  colnames(final_results) <- set_names
  
  for (next_name in set_names){
    matching <- paste0("SCORE_",next_name,"_",sig)
    relevant_table <- real_p_values[score == matching]
    setkey(relevant_table,test)
    
    for (j in phenotypes){
      relevant_p_value <- as.numeric(relevant_table[test == j,.(p)])
      relevant_random_sets <- fread(paste0(next_name,"_table_of_random_p_values.txt"), header = T)
      pattern <- paste0("\\_",sig,"$")
      permutation_p_value <- relevant_random_sets[grep(score_list, pattern = pattern, perl = T), sum(.SD < relevant_p_value)/1000 ,.SDcol = j]
      
      final_results[j,next_name] <- permutation_p_value
    }
  }
  assign(paste("final_results_",sig), final_results,envir = .GlobalEnv)
}


write.table(`final_results_ 0.05`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.05_OH_GOD_ITS_DONE_biobank.txt",col.names = T,row.names = T,quote = F)
write.table(`final_results_ 0.5`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.5_OH_GOD_ITS_DONE_biobank.txt",col.names = T,row.names = T,quote = F)

write.csv(`final_results_ 0.05`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.05_biobank.csv",col.names = T,row.names = T)
write.csv(`final_results_ 0.5`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.5_biobank.csv",col.names = T,row.names = T)

write.csv(Permutation_ALSPAC_0.05, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.05_ALSPAC.csv",col.names = T,row.names = T)
write.csv(Permutation_ALSPAC_0.5, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.5_ALSPAC.csv",col.names = T,row.names = T)

testing_4 <- grep(x = row.names(`final_results_ 0.05`),pattern = "thickavg")
new_table <- `final_results_ 0.05`[-testing_4,]

permutation_p_value

set.seed(1)


