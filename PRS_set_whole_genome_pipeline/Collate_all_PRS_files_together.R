############################################################
###########Collate all PRS profiles from workflow###########
############################################################
# source(paste0(path_to_PRS_scripts,'Pathway_analysis_scripts/Pathway_PRS_scoring.R'))

### Start Timer
ptm <- proc.time()

### Library
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

### environment for functions
e <- new.env()

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
Extra_analyses <- args[5]
Full_genome_PRS_extra_analysis <- args[6]
Gene_regions <- args[7]
whole_genome_genic <- args[8]
Gene_specific_PRS <- args[9]

significance_thresholds <- as.numeric(args[c(10:length(args))])
print(significance_thresholds)



### Write function combining Gene-set PRS into one table

# This checks and reads in the
combine_Gene_set_PRS <- function(Extra_analyses, Gene_regions, Training_name, Validation_name, ALL_PRS){
  if (Extra_analyses == TRUE){

    if (Gene_regions == "both" | Gene_regions == "extended"){

     Gene_regions_extended <- fread(paste0("FINAL_PATHWAY_RESULTS_PRS_PROFILES",Training_name,"_",Validation_name, "_extended.txt"))
     cols <- names(Gene_regions_extended)[4:ncol(Gene_regions_extended)]
     setnames(Gene_regions_extended, old = c(4:ncol(Gene_regions_extended)), new = paste("extended_geneset", cols, sep = "_"))

     if("FID" %in% colnames(ALL_PRS) & "IID" %in% colnames(ALL_PRS) & "PHENO" %in% colnames(ALL_PRS))
     {
        ALL_PRS <- merge(ALL_PRS, Gene_regions_extended, by = c("FID", "IID", "PHENO"))
     }else{
        ALL_PRS <- Gene_regions_extended
     }
    }

    if (Gene_regions == "both" | Gene_regions == "normal"){
      Gene_regions_normal <- fread(paste0("FINAL_PATHWAY_RESULTS_PRS_PROFILES",Training_name,"_",Validation_name, "_normal.txt"))
      cols <- names(Gene_regions_normal)[4:ncol(Gene_regions_normal)]
      setnames(Gene_regions_normal, old = c(4:ncol(Gene_regions_normal)), new = paste("normal_geneset", cols, sep = "_"))

      if("FID" %in% colnames(ALL_PRS) & "IID" %in% colnames(ALL_PRS) & "PHENO" %in% colnames(ALL_PRS))
      {
        ALL_PRS <- merge(ALL_PRS, Gene_regions_normal, by = c("FID", "IID", "PHENO"))
      }else{
        ALL_PRS <- Gene_regions_normal
      }
    }

     assign("ALL_PRS", ALL_PRS, envir = e)

  }else{
    stop("No Gene-Set PRS involved in this analysis")

}


}

combined_whole_genome_genic <- function(whole_genome_genic, Gene_regions, Training_name, Validation_name, significance_thresholds, ALL_PRS){
  if(whole_genome_genic == TRUE){
    if (Gene_regions == "both" | Gene_regions == "extended"){

      Location_of_whole_genome_extended_gene_files <- paste0(Training_name,"_",Validation_name,"_output/PRS_scoring/", Training_name,"_",Validation_name,"_extended_gene_regions_Clumped_whole_genome_final_significance_threshold_at_", significance_thresholds,".profile")

      my_data <- lapply(Location_of_whole_genome_extended_gene_files, fread, header=TRUE)
      names(my_data) <- str_replace(Location_of_whole_genome_extended_gene_files, pattern = ".profile", replacement = "")

      # iterate through significance thresholds
      for (i in 1:length(significance_thresholds)) {
        my_data[[i]] <- my_data[[i]][,c(1:3,6), with = FALSE]
        colnames(my_data[[i]]) <- c("FID", "IID", "PHENO", paste("SCORE_whole_genome_", significance_thresholds[i], sep=""))

      }

      all_prs_extended <- join_all(my_data, by=c("FID", "IID", "PHENO"), type='left')

      cols <- names(all_prs_extended)[4:ncol(all_prs_extended)]
      setnames(all_prs_extended, old = c(4:ncol(all_prs_extended)), new = paste("extended.genic.genome", cols, sep = "_"))


      if("FID" %in% colnames(ALL_PRS) & "IID" %in% colnames(ALL_PRS) & "PHENO" %in% colnames(ALL_PRS))
      {
        ALL_PRS <- merge(ALL_PRS, all_prs_extended, by = c("FID", "IID", "PHENO"))
      }else{
        ALL_PRS <- all_prs_extended
      }
    }

      if (Gene_regions == "both" | Gene_regions == "normal"){

        Location_of_whole_genome_normal_gene_files <- paste0(Training_name,"_",Validation_name,"_output/PRS_scoring/", Training_name,"_",Validation_name,"_normal_gene_regions_Clumped_whole_genome_final_significance_threshold_at_", significance_thresholds,".profile")

        my_data <- lapply(Location_of_whole_genome_normal_gene_files, fread, header=TRUE)
        names(my_data) <- str_replace(Location_of_whole_genome_normal_gene_files, pattern = ".profile", replacement = "")

        # iterate through significance thresholds
        for (i in 1:length(significance_thresholds)) {
          my_data[[i]] <- my_data[[i]][,c(1:3,6)]
          colnames(my_data[[i]]) <- c("FID", "IID","PHENO", paste("SCORE_whole_genome_", significance_thresholds[i], sep=""))
        }

        all_prs_normal <- join_all(my_data, by=c("FID", "IID", "PHENO"), type='left')

        cols <- names(all_prs_normal)[4:ncol(all_prs_normal)]
        setnames(all_prs_normal, old = c(4:ncol(all_prs_normal)), new = paste("normal.genic.genome", cols, sep = "_"))

        if("FID" %in% colnames(ALL_PRS) & "IID" %in% colnames(ALL_PRS) & "PHENO" %in% colnames(ALL_PRS))
        {
          ALL_PRS <- merge(ALL_PRS, all_prs_normal, by = c("FID", "IID", "PHENO"))
        }else{
          ALL_PRS <- all_prs_normal
        }
      }

      assign("ALL_PRS", ALL_PRS, envir = e)

    }else{
      stop("No Gene-Set PRS involved in this analysis")
    }
}

combined_full_genome <- function (Training_name, Validation_name, significance_thresholds, ALL_PRS){

    Location_of_whole_genome_entire_genome_files <-  paste0(Training_name,"_",Validation_name,"_output/PRS_scoring/",Training_name,"_",Validation_name,"_whole_genome_significance_threshold_at_",significance_thresholds,".profile")
    my_data <- lapply(Location_of_whole_genome_entire_genome_files, fread, header=TRUE)
    names(my_data) <- str_replace(Location_of_whole_genome_entire_genome_files, pattern = ".profile", replacement = "")

    # iterate through significance thresholds
    for (i in 1:length(significance_thresholds)) {
      my_data[[i]] <- my_data[[i]][,c(1:3,6)]
      colnames(my_data[[i]]) <- c("FID", "IID", "PHENO", paste("SCORE_whole_genome_", significance_thresholds[i], sep=""))
    }

    all_prs_whole_genome <- join_all(my_data, by=c("FID", "IID", "PHENO"), type='left')
    cols <- names(all_prs_whole_genome)[4:ncol(all_prs_whole_genome)]
    setnames(all_prs_whole_genome, old = c(4:ncol(all_prs_whole_genome)), new = paste("All.genome", cols, sep = "_"))

    if("FID" %in% colnames(ALL_PRS) & "IID" %in% colnames(ALL_PRS) & "PHENO" %in% colnames(ALL_PRS))
    {
      ALL_PRS <- merge(ALL_PRS, all_prs_whole_genome, by = c("FID", "IID", "PHENO"))
    }else{
      ALL_PRS <- all_prs_whole_genome
    }

    assign("ALL_PRS", ALL_PRS, envir = e)
}

ALL_PRS <- data.frame(x = NA, y = NA)
assign("ALL_PRS", ALL_PRS, envir = e)

combine_Gene_set_PRS(Extra_analyses, Gene_regions, Training_name, Validation_name, e$ALL_PRS)
combined_whole_genome_genic(whole_genome_genic, Gene_regions, Training_name, Validation_name, significance_thresholds, e$ALL_PRS)
combined_full_genome(Training_name, Validation_name, significance_thresholds, e$ALL_PRS)

write.table(e$ALL_PRS, file = paste0("Collated_PRS_analysis_for",Training_name,"_", Validation_name,".txt"), quote = F, row.names = F)
