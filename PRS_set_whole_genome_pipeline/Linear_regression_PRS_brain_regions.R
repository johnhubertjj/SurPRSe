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
Training_names <- c("PGC1_only_more_sig_thresh", "PGC1_plussweden_more_sig_thresh", "PGC2_more_sig_thresh", "CLOZUK_PGC2noclo_more_sig_thresh")
Validation_name <- "ALSPAC"
Scale_phenotypes <- T
Scale_PRS <- T
significance_thresholds <- c(5e-08, 1e-06, 1e-04, 0.01, 0.05, 0.1, 0.2, 0.5, 1)

for (Training_name in Training_names){ 
#Training_name <- "PGC2"
#Training_name <- "CLOZUK_PGC2noclo"
#Training_name <- "scz.swe.pgc1.results.v3"
#Training_name <- "HG19_pgc.scz.full.2012-04"
 
# Write in the phenotype files
Pheno_brain_thickness <- fread("ALSPAC_thickness_27january2017.txt")
Pheno_brain_volume <-fread("ALSPAC_volumes_27january2017.txt") 

# Collate all the scores together into one table
my_files <- paste0(Training_name, "_", Validation_name,"_output/PRS_scoring/", Training_name, "_", Validation_name, "_whole_genome_significance_threshold_at_", significance_thresholds, ".profile")
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



# Write the file name relevant to the results you have produced
#write.table(results, file = "all_scaled_PGC1_no_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt", row.names = F, quote = F)


### I would stop here Eleonora ####
#results1 <- fread("PGC1_no_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results2 <- fread("PGC1_with_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results3 <- fread ("PGC2_daner_file_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results4 <- fread("CLOZUK_PGC2noclo_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")


# Place the results into a list
brain_volume_progression <- list(`results_HG19_pgc.scz.full.2012-04`, results_scz.swe.pgc1.results.v3, results_PGC2, results_CLOZUK_PGC2noclo)
brain_volume_progression <- list(results_PGC1_only_more_sig_thresh, results_PGC1_plussweden_more_sig_thresh, results_PGC2_more_sig_thresh, results_CLOZUK_PGC2noclo_more_sig_thresh)
ICV_progression <- list(`results_ICVHG19_pgc.scz.full.2012-04`, results_ICVscz.swe.pgc1.results.v3, results_ICVPGC2, results_ICVCLOZUK_PGC2noclo)

#results_scaled_1 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_PGC1_no_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results_scaled_2 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_PGC1_with_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results")
#results_scaled_3 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_PGC2_daner_file_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
#results_scaled_4 <- fread("~/Documents/CLOZUK_PGC2noclo.METAL/all_scaled_CLOZUK_PGC2noclo_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")

#brain_volume_progression_scaled <- list(results_scaled_1, results_scaled_2, results_scaled_3, results_scaled_4)


## PLOTS OF BETA AND SE (AND CI) ## 
newdataframe_beta <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$estimate, brain_volume_progression[[2]]$estimate, brain_volume_progression[[3]]$estimate, brain_volume_progression[[4]]$estimate)
names(newdataframe_beta) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df <- t(newdataframe_beta)
colnames(test_df) <- as.character(test_df[1,])
test_df_beta <- test_df[-1,]
test_df_beta <- cbind(test_df_beta,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df_beta <- as.data.frame(test_df_beta)



number_of_phenotypes <- length(phenotypes)
number_of_thresholds <- length(significance_thresholds)

newdataframe_SE <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$SE, brain_volume_progression[[2]]$SE, brain_volume_progression[[3]]$SE, brain_volume_progression[[4]]$SE)
names(newdataframe_SE) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df_SE <- t(newdataframe_SE)
colnames(test_df_SE) <- as.character(test_df_SE[1,])
test_df_SE <- test_df_SE[-1,]
test_df_SE <- cbind(test_df_SE,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
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

## Scatter Plots with SE error bars ##
for (i in 1:number_of_phenotypes){
  # Plots for Betas and SE based on files #
  all_plots_current <- all_plots
  all_plots_current$Brain_region <- as.character(all_plots_current$Brain_region)
  
  all_plots_current <- as.data.table(all_plots_current)
  setkey(all_plots_current,"Brain_region")
  all_plots_current <- all_plots_current[phenotypes[i],]
  
  p <- ggplot(all_plots_current, aes(x=Dataset, y=BETA, fill = Pvaluethreshold, group=Dataset))
  
  p <- p +
    geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
    geom_point()
  
  p <- p + scale_x_discrete(labels=c("PGC1" = "PGC1", "PGC1swe" = "PGC1swe", "PGC2" = "PGC2", "QCLOZUK_PGC2noclo" = "CLOZUK"))
  p <- p + facet_grid(. ~ Pvaluethreshold) +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_discrete(guide=FALSE)
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(all_plots$Brain_region[i])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  
  gglist_0.05[[i]] <- p
}


setwd("~/Desktop/")


pdf("ALSPAC_PGC_all_brain_regions_all_sig_thresh_all_scaled.pdf", onefile = TRUE,paper = "A4")
invisible(lapply(gglist_0.05,print))
dev.off()



for (i in 1:44){
  # Plots for Betas and SE based on files #
  p <- ggplot(all_plots[c(i,i+44,i+88,i+132, i+176, i+220, i+264, i+308),], aes(x=Dataset, y=BETA, fill = Pvaluethreshold, group=Dataset))
  
  p <- p +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25)
  p <- p + scale_x_discrete(labels=c("PGC1" = "PGC1", "PGC1swe" = "PGC1swe", "PGC2" = "PGC2", "QCLOZUK_PGC2noclo" = "CLOZUK"))
  p <- p + facet_grid(. ~ Pvaluethreshold) +
    theme(strip.text.x = element_text(size = 10))
  p <- p + scale_fill_discrete(guide=FALSE)
  p <- p + theme(axis.text.x = element_text(size = 7))
  p <- p + ggtitle(all_plots$Brain_region[i])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  
  gglist_0.05[[i]] <- p
}







# BETA AND SE PLOTS SCALED ###
newdataframe_beta <- data.frame(brain_volume_progression_scaled[[1]]$test,brain_volume_progression_scaled[[1]]$estimate, brain_volume_progression_scaled[[2]]$estimate, brain_volume_progression_scaled[[3]]$estimate, brain_volume_progression_scaled[[4]]$estimate)

names(newdataframe_beta) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df <- t(newdataframe_beta)
colnames(test_df) <- as.character(test_df[1,])
test_df_beta <- test_df[-1,]
test_df_beta <- cbind(test_df_beta,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df_beta <- as.data.frame(test_df_beta)

BETA_dataframe1 <- newdataframe_beta[1:36,]
BETA_dataframe2 <- newdataframe_beta[37:72,]
test_beta <- stack(BETA_dataframe1)
names(test_beta) <- c("BETA", "Dataset")
test_beta$v3 <- rep(BETA_dataframe1$test, 4)

test_beta_0.5 <- stack(BETA_dataframe2)
names(test_beta_0.5) <- c("BETA", "Dataset")
test_beta_0.5$v3 <- rep(BETA_dataframe2$test, 4)


newdataframe_SE <- data.frame(brain_volume_progression_scaled[[1]]$test,brain_volume_progression_scaled[[1]]$SE, brain_volume_progression_scaled[[2]]$SE, brain_volume_progression_scaled[[3]]$SE, brain_volume_progression_scaled[[4]]$SE)
names(newdataframe_SE) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df_SE <- t(newdataframe_SE)
colnames(test_df_SE) <- as.character(test_df_SE[1,])
test_df_SE <- test_df_SE[-1,]
test_df_SE <- cbind(test_df_SE,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df_SE <- as.data.frame(test_df_SE)

SE_dataframe1 <- newdataframe_SE[1:36,]
SE_dataframe2 <- newdataframe_SE[37:72,]
test_SE <- stack(SE_dataframe1)
names(test_SE) <- c("SE", "Dataset")
test_SE$v3 <- rep(SE_dataframe1$test, 4)

test_SE_0.5 <- stack(SE_dataframe2)
names(test_SE_0.5) <- c("SE", "Dataset")
test_SE_0.5$v3 <- rep(SE_dataframe2$test, 4)

Beta_se_plots_0.5 <- cbind(test_beta_0.5,test_SE_0.5$SE)
Beta_se_plots_0.5$upper <- Beta_se_plots_0.5$BETA + Beta_se_plots_0.5$`test_SE_0.5$SE`
Beta_se_plots_0.5$lower <- Beta_se_plots_0.5$BETA - Beta_se_plots_0.5$`test_SE_0.5$SE`

Beta_se_plots <- cbind(test_beta,test_SE$SE)
Beta_se_plots$upper <- Beta_se_plots$BETA + Beta_se_plots$`test_SE$SE`
Beta_se_plots$lower <- Beta_se_plots$BETA - Beta_se_plots$`test_SE$SE`


#pdf("brain_area_plots_pthresh_0.5.pdf", onefile = TRUE)
#invisible(lapply(gglist_0.5,print))
#dev.off()

## PLOTS OF ICV ## 
newdataframe_beta <- data.frame(ICV_progression[[1]]$test,ICV_progression[[1]]$estimate, ICV_progression[[2]]$estimate, ICV_progression[[3]]$estimate, ICV_progression[[4]]$estimate)
names(newdataframe_beta) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df <- t(newdataframe_beta)
colnames(test_df) <- as.character(test_df[1,])
test_df_beta <- test_df[-1,]
test_df_beta <- cbind(test_df_beta,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df_beta <- as.data.frame(test_df_beta)

BETA_dataframe1 <- newdataframe_beta[1,]
BETA_dataframe2 <- newdataframe_beta[2,]
test_beta <- melt(BETA_dataframe1, id.vars="test")
names(test_beta) <- c("Brain_region", "Dataset","BETA")

test_beta_0.5 <- melt(BETA_dataframe2,id.vars="test")
names(test_beta_0.5) <- c("Brain_region", "Dataset","BETA")



newdataframe_SE <- data.frame(ICV_progression[[1]]$test,ICV_progression[[1]]$SE, ICV_progression[[2]]$SE, ICV_progression[[3]]$SE, ICV_progression[[4]]$SE)
names(newdataframe_SE) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df_SE <- t(newdataframe_SE)
colnames(test_df_SE) <- as.character(test_df_SE[1,])
test_df_SE <- test_df_SE[-1,]
test_df_SE <- cbind(test_df_SE,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df_SE <- as.data.frame(test_df_SE)

SE_dataframe1 <- newdataframe_SE[1,]
SE_dataframe2 <- newdataframe_SE[2,]
test_SE <- melt(SE_dataframe1, id.vars="test")
names(test_SE) <- c("Brain_region", "Dataset","SE")


test_SE_0.5 <- melt(SE_dataframe2, id.vars="test")
names(test_SE_0.5) <- c("Brain_region", "Dataset","SE")

test_SE_0.5$SE <- as.numeric(test_SE_0.5$SE)
test_SE$SE <- as.numeric(test_SE$SE)
test_beta$BETA <- as.numeric(test_beta$BETA)
test_beta_0.5$BETA <- as.numeric(test_beta_0.5$BETA)

Beta_se_plots_0.5 <- cbind(test_beta_0.5,test_SE_0.5$SE)
Beta_se_plots_0.5$upper <- Beta_se_plots_0.5$BETA + Beta_se_plots_0.5$`test_SE_0.5$SE`
Beta_se_plots_0.5$lower <- Beta_se_plots_0.5$BETA - Beta_se_plots_0.5$`test_SE_0.5$SE`

Beta_se_plots <- cbind(test_beta,test_SE$SE)
Beta_se_plots$upper <- Beta_se_plots$BETA + Beta_se_plots$`test_SE$SE`
Beta_se_plots$lower <- Beta_se_plots$BETA - Beta_se_plots$`test_SE$SE`

Beta_se_plots$P_Value_threshold <- c(rep(significance_thresholds[1],4))
Beta_se_plots_0.5$P_Value_threshold <- c(rep(significance_thresholds[2],4))

names(Beta_se_plots) <- c("Brain_region", "Dataset","BETA", "SE", "upper", "lower", "Pvaluethreshold")
names(Beta_se_plots_0.5) <- c("Brain_region", "Dataset","BETA", "SE", "upper", "lower", "Pvaluethreshold")

all_plots <- rbind(Beta_se_plots,Beta_se_plots_0.5)
all_plots$Pvaluethreshold <- as.factor(all_plots$Pvaluethreshold)
gglist_0.05 <- list() 

  # Plots for Betas and SE based on files #
  p <- ggplot(all_plots, aes(Dataset, BETA, fill = Pvaluethreshold, group= Dataset))
  
  p <- p +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25)
  p <- p + scale_x_discrete(labels=c("PGC1" = "PGC1", "PGC1swe" = "PGC1swe", "PGC2" = "PGC2", "QCLOZUK_PGC2noclo" = "CLOZUK"))
  p <- p + facet_grid(. ~ Pvaluethreshold) +
    theme(strip.text.x = element_text(size = 10))
  p <- p + scale_fill_discrete(guide=FALSE)
  p <- p + theme(axis.text.x = element_text(size = 7))
  p <- p + ggtitle(all_plots$Brain_region)
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  
  p
}


