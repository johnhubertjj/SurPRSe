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


#Training_names <- c("HG19_pgc.scz.full.2012-04", "scz.swe.pgc1.results.v3","PGC2","CLOZUK_PGC2noclo")

Training_names <- "CLOZUK_PGC2noclo"
Validation_name <- "ALSPAC"
Scale_phenotypes <- F
Scale_PRS <- T
significance_thresholds <- c(0.05, 0.5)


set_names <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268")

next_name <- "Calcium_ion_import_GO0070509"

for (next_name in set_names){
  

ALL_PRS_DATA <- fread(paste0("./FINAL_PATHWAY_RESULTS_PRS_PROFILESCLOZUK_PGC2noclo_ALSPAC_",next_name,"randomised_gene_sets.txt"),header = T)

for (Training_name in Training_names){ 
  #Training_name <- "PGC2"
  #Training_name <- "CLOZUK_PGC2noclo"
  #Training_name <- "scz.swe.pgc1.results.v3"
  #Training_name <- "HG19_pgc.scz.full.2012-04"
  
  # Write in the phenotype files
  
  Pheno_brain_thickness <- fread("/Users/johnhubert/Dropbox/ALSPAC_thickness_27january2017.txt")
  Pheno_brain_volume <-fread("/Users/johnhubert/Dropbox/ALSPAC_volumes_27january2017.txt") 
  
  
  if (Scale_PRS == T){
    current_dataframe <- ALL_PRS_DATA[,lapply(.SD, scale, center = T, scale = T) ,.SDcols = !c("FID","IID")]
    Identifiers <- ALL_PRS_DATA[,c("FID", "IID"), with = F]
    ALL_PRS_DATA <- cbind(Identifiers,current_dataframe)
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
  
  
  
  DATA <- join_all(list(ALL_PRS_DATA, DATA_pheno), by=c("FID"), type='inner')
  
  # create a score list
  score_list <- colnames(ALL_PRS_DATA[,-1:-2])
  
  phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])
  DATA$kz021 <- as.factor(DATA$kz021)
  
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
    
      fit <- lm(DATA[[j]] ~ DATA[[i]] + kz021 + age + ICV, data=DATA)
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

real_p_values <- fread("CLOZUK_PGC2noclo_Selected_SCZ_pathways.txt")
real_p_values <- real_p_values[,.(test,score,p)]
setkey(real_p_values,score)


for (sig in significance_thresholds){

final_results <- matrix(NA,nrow = 44,ncol=9)
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


write.table(`final_results_ 0.05`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.05_OH_GOD_ITS_DONE.txt",col.names = T,row.names = T,quote = F)
write.table(`final_results_ 0.5`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.5_OH_GOD_ITS_DONE.txt",col.names = T,row.names = T,quote = F)

write.csv(`final_results_ 0.05`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.05_OH_GOD_ITS_DONE.csv",col.names = T,row.names = T)
write.csv(`final_results_ 0.5`, "FINAL_RESULTS_OF_PERMUTATION_TESTS_OF_RANDOM_SNPS_0.5_OH_GOD_ITS_DONE.csv",col.names = T,row.names = T)

testing_4 <- grep(x = row.names(`final_results_ 0.05`),pattern = "thickavg")
new_table <- `final_results_ 0.05`[-testing_4,]

permutation_p_value

set.seed(1)


