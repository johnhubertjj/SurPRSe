### Start Timer
ptm <- proc.time()

#Print R version
getRversion()    

# devtools::install_github('hadley/ggplot2')

library(ggplot2)
library(gridExtra)
library(data.table)

# required packages -- will only install if they are not already installed
list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly","tools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

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
Pathway_directory <- args[5]
Pathway_file_name <- args[6] # the name of the file to be accessed (must be in stationary directory)
Gene_regions <- args[7]

# Define environment
e <- new.env()

### Define functions ###

# This function creates 3 stringed vectors recording the PRS score files that exist, the ones that don't and all together if there was at least one SNP in each

list_of_filenames <- function(Training_name, Validation_name, Pathway_directory, Gene_regions){
  
  Filenames_of_score_files <- NULL
  All_files_for_PRS <- NULL
  FALSE_Filenames_of_score_files <- NULL
  
  Scores_with_SNPs <- fread(paste0(Pathway_directory,"thresholds_to_use_",Gene_regions,".txt"))
  setnames(Scores_with_SNPs, c("Pathway", "Significance_Threshold", "file_exists"))
  Scores_with_SNPs_TRUE <- Scores_with_SNPs[file_exists == TRUE]
  Scores_with_SNPs_FALSE <- Scores_with_SNPs[file_exists == FALSE]
  
  for (i in 1:nrow(Scores_with_SNPs_TRUE)){
    scoring_output_file <- paste0(Pathway_directory,Scores_with_SNPs_TRUE$Pathway[i],"/scoring_",Training_name,"_",Validation_name,"_pathway_",Scores_with_SNPs_TRUE$Pathway[i],"_", Gene_regions,"_with_",Scores_with_SNPs_TRUE$Significance_Threshold[i],".profile")
    Filenames_of_score_files <- c(Filenames_of_score_files,scoring_output_file)
    #order_of_output_for_sig_thresh <- c(order_of_output_for_sig_thresh, paste0(pathway_names,"_",rep(i,length(scoring_output_file))))
  }
  
  for (i in 1:nrow(Scores_with_SNPs_FALSE)){
    scoring_output_file <- paste0(Pathway_directory,Scores_with_SNPs_FALSE$Pathway[i],"/scoring_",Training_name,"_",Validation_name,"_pathway_",Scores_with_SNPs_FALSE$Pathway[i],"_", Gene_regions,"_with_",Scores_with_SNPs_FALSE$Significance_Threshold[i],".profile")
    FALSE_Filenames_of_score_files <- c(FALSE_Filenames_of_score_files, scoring_output_file)
  }
  
  for (i in 1:nrow(Scores_with_SNPs)){
    scoring_output_file <- paste0(Pathway_directory,Scores_with_SNPs$Pathway[i],"/scoring_",Training_name,"_",Validation_name,"_pathway_",Scores_with_SNPs$Pathway[i],"_", Gene_regions,"_with_",Scores_with_SNPs$Significance_Threshold[i],".profile")
    All_files_for_PRS <- c(All_files_for_PRS, scoring_output_file)
  }
  
  
assign("Filenames_of_score_files", Filenames_of_score_files, envir = e)
assign("FALSE_Filenames_of_score_files", FALSE_Filenames_of_score_files, envir = e)
assign("All_files_for_PRS", All_files_for_PRS, envir = e)

}

# This function creates the polygenic risk score table without dummy tables
create_collated_polygenic_risk_scores <- function(Gene_regions){
 
  my_data <- lapply(e$Filenames_of_score_files, read.table, header=TRUE) 
  names(my_data) <- str_replace(e$Filenames_of_score_files, pattern = ".profile", replacement = "")
  
  # iterate through significance thresholds 
  for (i in 1:length(e$Filenames_of_score_files)) {
    my_data[[i]] <- my_data[[i]][,c(1,2,6)]
    name <- names(my_data[i])
    score_significance_thresh <- gsub(".*\\_with\\_(.*$)", "\\1",x = name)
    score_name_pattern <- paste0(".*pathway_(.*?)_",Gene_regions,".*")
    score_name <- gsub(score_name_pattern, "\\1", x = name)
    colnames(my_data[[i]]) <- c("FID", "IID", paste("SCORE_", score_name,"_",score_significance_thresh, sep=""))
  }
  
  all_prs <- join_all(my_data, by=c("FID", "IID"), type='left')
  assign("PRS_profiles_with_SNPS", all_prs, envir = e)
}

# This function creates dummy polygenic risk score profiles in gene-sets where no polygenic risk scores exist
create_dummy_table_for_null_profiles <- function(Gene_regions){
  dummy_table <- matrix(nrow = nrow(e$PRS_profiles_with_SNPS), ncol = (2 + length(e$FALSE_Filenames_of_score_files)))
  dummy_table <- as.data.frame(dummy_table)
  dummy_table[,c("V1", "V2")] <- e$PRS_profiles_with_SNPS[,c("FID", "IID")]
  
  new_false_column_names <- rep(0,length(e$FALSE_Filenames_of_score_files))
  
  for (i in 1:length(e$FALSE_Filenames_of_score_files)){
    name <- e$FALSE_Filenames_of_score_files[i]
    score_significance_thresh <- gsub(".*\\_with\\_(.*?).profile", "\\1",x = name)
    score_name_pattern <- paste0(".*pathway_(.*?)_",Gene_regions,".*")
    score_name <- gsub(score_name_pattern, "\\1", x = name)
    new_false_column_names[i] <- paste("SCORE_", score_name,"_",score_significance_thresh, sep="")
  }
  
  colnames(dummy_table) <- c("FID","IID", new_false_column_names)
  assign("PRS_dummy_profiles_without_SNPS", dummy_table, envir = e)
}

### Create collated PRS profiles for the normal gene regions ###

if (Gene_regions == "both" | Gene_regions == "normal"){
  
  extra_gene_regions <- "normal"
  
  list_of_filenames(Training_name, Validation_name, Pathway_directory, extra_gene_regions)
  create_collated_polygenic_risk_scores(extra_gene_regions)
  create_dummy_table_for_null_profiles(extra_gene_regions)
  
  Final_normal_gene_regions_table <- merge(e$PRS_profiles_with_SNPS,e$PRS_dummy_profiles_without_SNPS, by = c("FID", "IID"))

  write.table(Final_normal_gene_regions_table,file = paste0("FINAL_PATHWAY_RESULTS_PRS_PROFILES",Training_name,"_",Validation_name, "_", extra_gene_regions, ".txt"),quote = F,row.names = F)
}

### Create collated PRS profiles for the extended gene regions ###

if (Gene_regions == "both" | Gene_regions == "extended"){
  
  extra_gene_regions <- "extended"
  
  list_of_filenames(Training_name, Validation_name, Pathway_directory, extra_gene_regions)
  
  create_collated_polygenic_risk_scores(extra_gene_regions)
  
  create_dummy_table_for_null_profiles(extra_gene_regions)
  
  Final_normal_gene_regions_table <- merge(e$PRS_profiles_with_SNPS,e$PRS_dummy_profiles_without_SNPS, by = c("FID", "IID"))
  
  write.table(Final_normal_gene_regions_table,file = paste0("FINAL_PATHWAY_RESULTS_PRS_PROFILES",Training_name,"_",Validation_name, "_", extra_gene_regions, ".txt"),quote = F,row.names = F)
}




