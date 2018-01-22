### Start Timer
ptm <- proc.time()

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

##################################################
# Checking location for serial or batch analysis #
##################################################

library (data.table)
library (parallel)
library (base)

library(ggplot2)
library(gridExtra)
library(data.table)

# required packages -- will only install if they are not already installed
list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly","tools")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# loads the required packages
lapply(list.of.packages, require, character.only = TRUE)

# Set independent environment when needed #
e <- new.env()

# specify the different input tables #
Training_name <- args[3]
Validation_name <- args[4]
Gene_output_directory <- args[5]
Pathway_output_directory <- args[6]
Randomised_output_directory <- args[7]
gene_loc_file_name <- args[8]
Gene_regions <- args[9] # Whether to include/exclude the regulatory regions of a gene
rand_n = args[10]; # Number of random sets to generate for each gene-set
Pathway_file_name <- args[11] # The input file annotating genes to pathways
Pathways <- as.character(args[c(12:length(args))])

## Variables_useful_for_testing ## 

# Training_name <- "CLOZUK_PGC2noclo"
# Validation_name <- "ALSPAC"
# Gene_output_directory <- paste0(Training_name,"_", Validation_name, "_output/Genes/")
# Pathway_output_directory <- paste0(Training_name,"_",Validation_name,"_output/Pathways/")
# Randomised_output_directory <- paste0(Training_name,"_", Validation_name, "_output/Pathways/Randomised_gene_sets_analysis/Scores/")
# gene_loc_file_name <- "/home/c1020109/Stationary_data/NCBI37.3.gene.loc"
# Gene_regions <- "both"
# rand_n = 1000; # Number of random sets to generate for each gene-set
# Pathway_file_name <- "/home/c1020109/Stationary_data/Selected_Pocklington_plus_GO_pathways_SCHIZ.txt"
# Pathways <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268") 

# Read in array variables from file and rearrage to a vector#
significance_thresholds <- fread(paste0(Training_name,"_", Validation_name,"_plink_significance_thresholds_arguments_file_tmp.txt"))      
significance_thresholds <- as.character(unlist(significance_thresholds$V3))
print(significance_thresholds)

significance_thresholds <- c(0.05,0.5)

Pathways_used <- file_path_sans_ext(basename(Pathway_file_name))

##### read in pathway sets and standardise column names
pathway_sets <- fread(Pathway_file_name)
setnames(pathway_sets, c("Pathway", "Gene"))

# extract names of each pathway
factorise_column_1<- as.factor(pathway_sets$Pathway)
pathway_names <- levels(factorise_column_1)
number_of_pathways_to_analyse <- length(pathway_names)

Randomised_output_directory <- paste0(Pathway_output_directory,"Randomised_gene_sets_analysis/Scores/")
final_scoring_output_file <- NULL
order_of_output_for_sig_thresh <- NULL

Collate_all_profiles <- function(i, pathway_names, Randomised_output_directory, significance_thresholds, rand_n){

# Select the pathway #  
pathway <- pathway_names[i]

# add in an error catch and write to directory for p-value thresholds that do not have score files! ###
# alter the significance threshold file accordingly (aka write over it or use a new one from now onwards ###

# Create the column order for the large data.table #
for (w in 1:length(significance_thresholds)){
  for( col in 1:rand_n){
  scoring_output_file <- paste0(Randomised_output_directory,pathway,"_random_",col,"_with_",significance_thresholds[w],".profile")
  final_scoring_output_file <- c(final_scoring_output_file,scoring_output_file)
  order_of_output_for_sig_thresh <- c(order_of_output_for_sig_thresh, paste0(pathway,"_",rep(col,length(scoring_output_file)),"_threshold_",rep(significance_thresholds[w],length(scoring_output_file))))
}
}

# create a list of all the random SNP sets for each profile #
my_data <- lapply(final_scoring_output_file, fread, header=TRUE) 
names(my_data) <- str_replace(final_scoring_output_file, pattern = ".profile", replacement = "")

# iterate through significance thresholds and random sets, re-name column headings to a simpler format
for (i in 1:length(final_scoring_output_file)) {
  my_data[[i]] <- my_data[[i]][,c(1,2,6), with = F]
  colnames(my_data[[i]]) <- c("FID", "IID", paste("SCORE_", order_of_output_for_sig_thresh[i], sep=""))
}

# Join all list items into one large data.table
all_prs <- join_all(my_data, by=c("FID", "IID"), type='left')

# Write to the output directory
write.table(all_prs,file = paste0("FINAL_PATHWAY_RESULTS_PRS_PROFILES", Training_name, "_", Validation_name, "_", pathway, "randomised_gene_sets.txt"),quote = F,row.names = F)
}

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

# Export the environment to the cluster
clusterExport(cl, "e")

# Create one large data.table for all polygenic risk score profiles for each random set, for each significance threshold
parLapply(cl, 1:number_of_pathways_to_analyse, Collate_all_profiles, pathway_names, Randomised_output_directory, significance_thresholds, rand_n)
stopCluster(cl)

#End Timer
proc.time() - ptm

