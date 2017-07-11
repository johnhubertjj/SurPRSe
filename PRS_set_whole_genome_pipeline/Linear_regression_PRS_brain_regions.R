### Start Timer
ptm <- proc.time()

#Print R version
getRversion()    

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
significance_thresholds <- as.numeric(args[c(5:length(args))])
print(significance_thresholds)

# Write in the phenotype files
Pheno_brain_thickness <- fread("ALSPAC_thickness_27january2017.txt")
Pheno_brain_volume <-fread("ALSPAC_volumes_27january2017.txt") 

# Collate all the scores together into one table
my_files <- paste0(Training_name,"_",Validation_name,"_output/PRS_scoring/",Training_name,"_",Validation_name,"_whole_genome_significance_threshold_at_",significance_thresholds,".profile")
my_data <- lapply(my_files, read.table, header=TRUE) 
names(my_data) <- str_replace(my_files, pattern = ".profile", replacement = "")

# iterate through significance thresholds 
for (i in 1:length(significance_thresholds)) {
  my_data[[i]] <- my_data[[i]][,c(1,2,6)]
  colnames(my_data[[i]]) <- c("FID", "IID", paste("SCORE_", significance_thresholds[i], sep=""))
}
all_prs <- join_all(my_data, by=c("FID", "IID"), type='left')

# Join all the tables together
setnames(Pheno_brain_volume,old = "SubjID", new = "FID")
setnames(Pheno_brain_thickness,old = "SubjID", new = "FID")

DATA_pheno <- join_all(list(Pheno_brain_thickness, Pheno_brain_volume), by=c("FID", "kz021", "age", "ICV"), type='inner')
DATA <- join_all(list(all_prs, DATA_pheno), by=c("FID"), type='inner')

# create a score list
score_list <- colnames(all_prs[,-1:-2])

# Add names of phenotypes corresponding to the table
phenotypes <- colnames(DATA_pheno[,5:ncol(DATA_pheno)])



DATA$kz021 <- as.factor(DATA$kz021)

# Remove individuals with missing values
#remove_ICV <- which(is.na(DATA$ICV))
#DATA <- DATA[-remove_ICV,]

# dont' change obj it makes the loop output a dataframe with the regression results
obj <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0)

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
    tmp3 <- c(j,i,tmp[2,], true_r2)
    obj <- rbind(obj, tmp3)
  }
}



# this is a clean-up step - do not change
results <- obj[which(obj$score %in% score_list),]

# Write the file name relevant to the results you have produced
write.table(results, file = "CLOZUK_PGC2noclo_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt", row.names = F, quote = F)


### I would stop here Eleon
results1 <- fread("PGC1_no_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
results2 <- fread("PGC1_with_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
results3 <- fread ("PGC2_daner_file_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
results4 <- fread("CLOZUK_PGC2noclo_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
brain_volume_progression <- list(results1, results2, results3, results4)



newdataframe <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$p, brain_volume_progression[[2]]$p, brain_volume_progression[[3]]$p, brain_volume_progression[[4]]$p)
names(newdataframe) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df <- t(newdataframe)
colnames(test_df) <- as.character(test_df[1,])
test_df <- test_df[-1,]
test_df <- cbind(test_df,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df <- as.data.frame(test_df)

library(plotly)
newdataframe1 <- newdataframe[1:36,]
newdataframe2 <- newdataframe[37:72,]
test <- stack(newdataframe1)
names(test) <- c("P_value", "Dataset")
test$v3 <- rep(newdataframe1$test, 4)

p <- plot_ly(test,x=~Dataset, y = ~P_value, type = 'scatter')

  for (i in 1:36){
    p <- add_trace(p, x = test_df$V73, y = test_df[,i], name = phenotypes[i], mode = "lines",type = 'scatter')
  }
p

test <- stack(newdataframe2)

names(test) <- c("P_value", "Dataset")
test$v3 <- rep(newdataframe2$test, 4)
p2 <- plot_ly(test,x=~Dataset, y = ~P_value, type = 'scatter')

for (i in 37:72){
  p2 <- add_trace(p2, x = test_df$V73, y = test_df[,i], name = phenotypes[i-36], mode = "lines",type = 'scatter')
}
p2



r2_dataframe <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$r.squared, brain_volume_progression[[2]]$r.squared, brain_volume_progression[[3]]$r.squared, brain_volume_progression[[4]]$r.squared)
names(r2_dataframe) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df <- t(r2_dataframe)
colnames(test_df) <- as.character(test_df[1,])
test_df <- test_df[-1,]
test_df <- cbind(test_df,c("PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo"))
test_df <- as.data.frame(test_df)


r2dataframe1 <- r2_dataframe[1:36,]
r2dataframe2 <- r2_dataframe[37:72,]
test <- stack(r2dataframe1)
names(test) <- c("rsquared", "Dataset")
test$v3 <- rep(r2dataframe1$test, 4)

r2_1 <- plot_ly(test,x=~Dataset, y = ~rsquared, type = 'scatter')

for (i in 1:36){
  r2_1 <- add_trace(r2_1, x = test_df$V73, y = test_df[,i], name = phenotypes[i], mode = "lines",type = 'scatter')
}
r2_1

test <- stack(r2dataframe2)

names(test) <- c("rsquared", "Dataset")
test$v3 <- rep(r2dataframe2$test, 4)
r2_2 <- plot_ly(test,x=~Dataset, y = ~rsquared, type = 'scatter')

for (i in 37:72){
  r2_2 <- add_trace(r2_2, x = test_df$V73, y = test_df[,i], name = phenotypes[i-36], mode = "lines",type = 'scatter')
}
r2_2

plotly_POST(p, filename = "brain_regions_PRS_pthres0.05_pvalue", sharing = "secret")
plotly_POST(p2, filename = "brain_regions_PRS_pthres0.5_pvalue", sharing = "secret")

plotly_POST(r2_1, filename = "brain_regions_PRS_pthres0.05_r2", sharing = "secret")
plotly_POST(r2_2, filename = "brain_regions_PRS_pthres0.5_r2", sharing = "secret")

newdataframe_beta <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$estimate, brain_volume_progression[[2]]$estimate, brain_volume_progression[[3]]$estimate, brain_volume_progression[[4]]$estimate)
names(newdataframe_beta) <- c("test", "PGC1", "PGC1swe", "PGC2", "QCLOZUK_PGC2noclo")
test_df <- t(newdataframe_beta)
colnames(test_df) <- as.character(test_df[1,])
test_df_beta <- test_df_beta[-1,]
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


newdataframe_SE <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$SE, brain_volume_progression[[2]]$SE, brain_volume_progression[[3]]$SE, brain_volume_progression[[4]]$SE)
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

gglist_0.05 <- list() 

for (i in 1:36){
# Plots for Betas and SE based on files #
p <- ggplot(Beta_se_plots[c(i,i+36,i+72,i+108),], aes(Dataset, BETA, fill = v3))
p <- p +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25)

gglist_0.05[[i]] <- p
}

gglist_0.5 <- list()
for (i in 1:36){

p <- ggplot(Beta_se_plots_0.5[c(i,i+36,i+72,i+108),], aes(Dataset, BETA, fill = v3))
p <- p +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25)
gglist_0.5[[i]] <- p
}

setwd("~/Desktop/")

library(gridExtra)

pdf("brain_area_plots_pthresh_0.05.pdf", onefile = TRUE)
invisible(lapply(gglist_0.05,print))
dev.off()

pdf("brain_area_plots_pthresh_0.5.pdf", onefile = TRUE)
invisible(lapply(gglist_0.5,print))
dev.off()




l <- mget(gglist_0.05)
lapply(ggplots,function(x){ggsave(file=paste(x,"pdf",sep="."),get(x))})
