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


# dont' change obj it makes the loop output a dataframe with the regression results
obj <- data.frame(test=0, score=0, estimate=0, SE=0, tvalue=0, p=0, r.squared=0)

# this example if for a linear regression (the phenotype of interest is a quantiative trait)
# is using a discrete phenotype, a logistic regression needs to be run, and the code altered from 'lm' to 'glm' including the argument of 'family = binomial'
# alterations for the calculation of R2 will also need to be made using the command highlighted above
for (i in score_list) {
  for (j in phenotypes) {
    fit <- lm(DATA[,j] ~ DATA[,i] + DATA$kz021 + DATA$age + DATA$ICV, data=DATA)
    fit1 <- lm(DATA[,j] ~ DATA$kz021 + DATA$age + DATA$ICV, data=DATA)
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

write.table(results, file = "PGC2_daner_file_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt", row.names = F, quote = F)

results1 <- fread("PGC1_no_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
results2 <- fread("PGC1_with_sweden_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")
results3 <- fread ("PGC2_daner_file_HG19_UCSC_ALSPAC_brain_regions_association_PRS_results.txt")

brain_volume_progression <- list(results1, results2, results3)

trace_0 <- rnorm(100, mean = 5)
trace_1 <- rnorm(100, mean = 0)
trace_2 <- rnorm(100, mean = -5)
x <- c(1:100)

newdataframe <- data.frame(brain_volume_progression[[1]]$test,brain_volume_progression[[1]]$p, brain_volume_progression[[2]]$p, brain_volume_progression[[3]]$p)
names(newdataframe) <- c("test", "PGC1", "PGC1swe", "PGC2")
test_df <- t(newdataframe)
colnames(test_df) <- as.character(test_df[1,])
test_df <- test_df[-1,]
test_df <- cbind(test_df,c("PGC1", "PGC1swe", "PGC2"))
test_df <- as.data.frame(test_df)

library(plotly)
newdataframe1 <- newdataframe[1:36,]
newdataframe2 <- newdataframe[37:72,]
test <- stack(newdataframe1)
test$v3 <- rep(newdataframe1$test, 3)

p <- plot_ly(test,x=~ind, y = ~values, type = 'scatter')

  for (i in 1:36){
    p <- add_trace(p, x = test_df$V73, y = test_df[,i], name = phenotypes[i], mode = "lines",type = 'scatter')
  }
p

test <- stack(newdataframe2)
test$v3 <- rep(newdataframe1$test, 3)

p2 <- plot_ly(test,x=~ind, y = ~values, type = 'scatter')

for (i in 1:36){
  p2 <- add_trace(p2, x = test_df$V73, y = test_df[,i], name = phenotypes[i], mode = "lines",type = 'scatter')
}
p


