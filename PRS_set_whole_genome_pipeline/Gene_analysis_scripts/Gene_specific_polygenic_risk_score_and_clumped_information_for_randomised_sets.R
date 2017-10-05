#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

# PRS_scoring_and_clumped_SNPs_info_whole_genome_genic
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

# specify the different input tables #

Training_name <- args[3]
Validation_name <- args[4]
Validation_full_name_serial <- args[5]
Gene_file_name <- args[6] # the name of the file to be accessed (must be in stationary directory)
Gene_directory <- args[7]
gene_loc_file_name <- args[8] # the name of the file containing gene locations (magma's version)
Gene_regions <- args[9] # whether to include/exclude the regulatory regions of a gene
chromosomes_to_analyse <- args[10]


setwd(".")

Training_name <- "CLOZUK_PGC2noclo"
Validation_name <- "ALSPAC"
Validation_full_name_serial <- "ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr"
Gene_file_name <- "~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_CLOZUK_PGC2noclo_extended_gene_regions_Clumped_whole_genome_final.bim"
Gene_directory <- "Genes"
gene_loc_file_name <-"/Users/johnhubert/Dropbox/Stationary_data/NCBI37.3.gene.loc"
Gene_regions <- "both"
chromosomes_to_analyse <- 22

#setwd("~/Documents/CLOZUK_ALSPAC_PATHWAY_TESTING")

# Create new variables based on input files to make things easier to read # 
output_directory <- paste0("./", Training_name, "_", Validation_name, "_output/", Gene_directory, "/")
Date <- Sys.Date()

# Reading_in_Pathway_files
library(data.table)
library(reshape2)

### environment for functions
e <- new.env()

### Function which assigns genes to the SNP data ####
### Currently designed for the input from bim files ###

### Function which adds duplicate SNPs which happen to be inside other genes
### Currently designed for extended gene regions ONLY
### requires the run of gene.out R script to obtain the specific gene.annot file for each pathway and chromosome
adding_unread_genes <- function(MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number){
  
  Melted_MAGMA_list <- melt(y)
  names(Melted_MAGMA_list) <- c("SNP", "Gene")
  Melted_MAGMA_list$SNP <- as.character(Melted_MAGMA_list$SNP) 
  clumped_SNPs$SNP <- as.character(clumped_SNPs$SNP)
  merged_table_one <- merge(clumped_SNPs, Melted_MAGMA_list, by = "SNP", all = F)
  merged_table_one$Gene <- as.numeric(merged_table_one$Gene)
  Gene_clumped_SNPs <- merge(MAGMA.gene.regions,merged_table_one, by = "Gene",all.y = T)
  Gene_clumped_SNPs[,"CHR.y" :=NULL]  # remove extra_column
  setcolorder(Gene_clumped_SNPs, c("CHR.x","SNP", "GD", "BP", "A1","A2","GENE_NAME","Gene", "BP_START", "BP_END", "BP_start_extended","BP_end_extended","STRAND"))
  setnames(Gene_clumped_SNPs, c("CHR.x","Gene"), c("CHR","GENE_NUMBER"))
  assign("test_data_frame", Gene_clumped_SNPs, envir = e)
  
}



# input the name of the pathways from the orignal file
#useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning|memory|conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities|coordination|movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")

#read in general MAGMA annotation file
MAGMA.gene.regions <- fread(gene_loc_file_name, colClasses = c("numeric","character", rep("numeric",2), rep("character",2)))
setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))

current_table_name <- MAGMA.gene.regions
setkey(current_table_name,STRAND) 

# equivalent to MAGMAs window option without ignoring strand
current_table_name <- current_table_name[STRAND == "+",BP_start_extended := BP_START - 35000]
current_table_name <- current_table_name[STRAND == "+",BP_end_extended := BP_END + 10000]
current_table_name <- current_table_name[STRAND == "-",BP_start_extended := BP_START - 10000]
current_table_name <- current_table_name[STRAND == "-",BP_end_extended := BP_END + 35000]

# setkey to the chromosome and remove all SNPs in the X chromosome to save time
setkey(current_table_name, CHR)
current_table_name <- current_table_name[!"X"]




if (Gene_regions == "both") {
  
#### Regular gene regions ####

Gene_regions <- "normal"
Gene_file_name <-  paste0(output_directory, Validation_name, "_", Training_name, "_", Gene_regions,"_gene_regions_Clumped_whole_genome_final.bim")
selecting_chromosomes <- fread(Gene_file_name)
names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")

# if there are no genes within that pathway, ignore it
# Ideally want to use the input from last time to check but for now we already have a record, so will check at the end

if(nrow(selecting_chromosomes) == 0) {
  message <- paste0("chromosome ", chromosomes_to_analyse, " contains no genes")
  write(message, file = paste0(output_directory,"remove_gene_analysis_empty_chromosomes_info_file_run2.txt"), append = T)
  next()
}

selecting_chromosomes_BP <- selecting_chromosomes$BP
  
  #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "regular")

  # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes

  GENES_to_snps <- scan(file = paste0(output_directory, Training_name, "_", Validation_name, "_SNPs_", Gene_regions, "_clumped_gene_temp.genes.annot"), what = "", sep = "\n")

  y <- strsplit(GENES_to_snps, "[[:space:]]+")
  names(y) <- sapply(y, '[[', 1)
  y <- lapply(y, '[', -1)
  
  y[[1]] <- NULL
  y[[1]] <- NULL
  
  adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
  
  setkey(e$test_data_frame, GENE_NUMBER)
  
  Markers_per_MB_non_independent_normal <- matrix('NA', nrow = nrow(current_table_name), ncol = 6)
  Markers_per_MB_non_independent_normal <- as.data.table(Markers_per_MB_non_independent_normal)
  #Markers_per_MB_non_independent_normal[, names(Markers_per_MB_non_independent_normal) := lapply(.SD, as.character)]
  Genes_to_sort_through_normal <- unlist(current_table_name[,.(Gene)])
  
  for (k in 1:length(Genes_to_sort_through_normal)){
    rows <- nrow(e$test_data_frame[.(Genes_to_sort_through_normal[k])])
    Length_of_gene <- e$test_data_frame[.(Genes_to_sort_through_normal[k]),.(BP_START,BP_END, BP_start_extended, BP_end_extended)][1]
    #Markers_per_MB_non_independent_normal[k,] <- c(k,rows,Length_of_gene)
    testing_2 <- as.matrix(Length_of_gene[1,])
    test <- unlist(c(Genes_to_sort_through_normal[k],rows,Length_of_gene[1,]))
    test <- as.character(test)
    for (j in seq_len(ncol(Markers_per_MB_non_independent_normal))){
      set(Markers_per_MB_non_independent_normal,k,j,test[j])
    }
  }
  
  names(Markers_per_MB_non_independent_normal) <- c("Gene","Nmarkers_in_Gene","BP_START","BP_END","BP_start_extended","BP_end_extended") 
  Markers_per_MB_non_independent_genes_wo_SNPS_normal <- Markers_per_MB_non_independent_normal[is.na(BP_START),.(Gene)]
  Markers_per_MB_non_independent_normal <- Markers_per_MB_non_independent_normal[!is.na(Markers_per_MB_non_independent_normal$BP_START)]
  
  e$test_data_frame[,c("BP_start_extended","BP_end_extended") := NULL]
  which(duplicated(e$test_data_frame$SNP,fromLast = T))
  
  assign("SNPs_for_clumping_normal", unique(e$test_data_frame$SNP), envir = e)
  assign("Gene_regions_annotation_table_normal", e$test_data_frame, envir = e)
  
  
  
  
  #### Extended gene regions ####
  #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "extended")
  
  Gene_regions <- "extended"
  Gene_file_name <-  paste0(output_directory, Validation_name, "_", Training_name, "_", Gene_regions,"_gene_regions_Clumped_whole_genome_final.bim")
  selecting_chromosomes <- fread(Gene_file_name)
  names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
  
  # if there are no genes within that pathway, ignore it
  # Ideally want to use the input from last time to check but for now we already have a record, so will check at the end
  
  if(nrow(selecting_chromosomes) == 0) {
    message <- paste0("chromosome ", chromosomes_to_analyse, " contains no genes")
    write(message, file = paste0(output_directory,"remove_gene_analysis_empty_chromosomes_info_file_run2.txt"), append = T)
    next()
  }
  
  selecting_chromosomes_BP <- selecting_chromosomes$BP
  
  #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "regular")
  
  # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
    GENES_to_snps <- scan(file = paste0(output_directory,Training_name, "_", Validation_name, "_SNPs_",Gene_regions,"_clumped_gene_temp.genes.annot"), what = "", sep = "\n")

  
  
  y <- strsplit(GENES_to_snps, "[[:space:]]+")
  names(y) <- sapply(y, '[[', 1)
  y <- lapply(y, '[', -1)
  
  
  y[[1]] <- NULL
  y[[1]] <- NULL
  
  adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
  
  setkey(e$test_data_frame, GENE_NUMBER)
  
  Markers_per_MB_non_independent_extended <- matrix('NA', nrow = nrow(current_table_name), ncol = 6)
  Markers_per_MB_non_independent_extended <- as.data.table(Markers_per_MB_non_independent_extended)
  #Markers_per_MB_non_independent_extended[, names(Markers_per_MB_non_independent_extended) := lapply(.SD, as.character)]
  Genes_to_sort_through_extended <- unlist(current_table_name[,.(Gene)])
  
  for (k in 1:length(Genes_to_sort_through_extended)){
    rows <- nrow(e$test_data_frame[.(Genes_to_sort_through_extended[k])])
    Length_of_gene <- e$test_data_frame[.(Genes_to_sort_through_extended[k]),.(BP_START,BP_END, BP_start_extended, BP_end_extended)][1]
    #Markers_per_MB_non_independent_extended[k,] <- c(k,rows,Length_of_gene)
    testing_2 <- as.matrix(Length_of_gene[1,])
    test <- unlist(c(Genes_to_sort_through_extended[k],rows,Length_of_gene[1,]))
    test <- as.character(test)
    for (j in seq_len(ncol(Markers_per_MB_non_independent_extended))){
      set(Markers_per_MB_non_independent_extended,k,j,test[j])
    }
  }
  
  names(Markers_per_MB_non_independent_extended) <- c("Gene","Nmarkers_in_Gene","BP_START","BP_END","BP_start_extended","BP_end_extended") 
  Markers_per_MB_non_independent_genes_wo_SNPS_extended <- Markers_per_MB_non_independent_extended[is.na(BP_start_extended),.(Gene)]
  Markers_per_MB_non_independent_extended <- Markers_per_MB_non_independent_extended[!is.na(Markers_per_MB_non_independent_extended$BP_start_extended)]
  
  e$test_data_frame[,c("BP_START","BP_END") := NULL]
  setnames(e$test_data_frame,old = c("BP_start_extended","BP_end_extended"), new = c("BP_START","BP_END"))
  
  which(duplicated(e$test_data_frame$SNP,fromLast = T))
  
  assign("SNPs_for_clumping_extended", unique(e$test_data_frame$SNP), envir = e)
  assign("Gene_regions_annotation_table_extended", e$test_data_frame, envir = e)
  
  print("all_okay!")
  print(Markers_per_MB_non_independent_normal)
  print(Markers_per_MB_non_independent_genes_wo_SNPS_normal)
  
  write.table(e$Gene_regions_annotation_table_normal, file = paste0(output_directory, Validation_name, "_", Training_name, "_normal_gene_regions_data_table_after_clumping.txt"), quote = F, row.names = F)
  write(e$SNPs_for_clumping_normal, file = paste0(output_directory, Validation_name, "_", Training_name, "_full_genome_SNPs_after_clumping_normal_gene_regions.txt"))
  write.table(Markers_per_MB_non_independent_normal, file = paste0(output_directory, Validation_name, "_", Training_name, "_normal_gene_region_information_for_randomisation_tests_ater_clumping.txt"), quote = F, row.names = F)
  write.table(Markers_per_MB_non_independent_genes_wo_SNPS_normal, file = paste0(output_directory, Validation_name, "_", Training_name,"_Normal_Genes_without_SNPs_annotated_to_them_full_genome_after_clumping.txt"), quote = F,row.names = F,col.names = F)
  
  write.table(e$Gene_regions_annotation_table_extended, file = paste0(output_directory, Validation_name, "_", Training_name, "_extended_gene_regions_data_table_after_clumping.txt"), quote = F, row.names = F)
  write(e$SNPs_for_clumping_extended, file = paste0(output_directory, Validation_name, "_", Training_name,"_full_genome_SNPs_after_clumping_extended_gene_regions.txt"))
  write.table(Markers_per_MB_non_independent_extended, file = paste0(output_directory, Validation_name, "_", Training_name, "_extended_gene_region_information_for_randomisation_tests_after_clumping.txt"), quote = F, row.names = F)
  write.table(Markers_per_MB_non_independent_genes_wo_SNPS_extended, file = paste0(output_directory, Validation_name, "_", Training_name, "_Extended_Genes_without_SNPs_annotated_to_them_full_genome_after_clumping.txt"), quote = F,row.names = F,col.names = F)

  Gene_regions <- "both"
}



if (Gene_regions == "normal") {
  
  Gene_file_name <-  paste0(output_directory, Validation_name, "_", Training_name, "_", Gene_regions,"_gene_regions_Clumped_whole_genome_final.bim")
  selecting_chromosomes <- fread(Gene_file_name)
  names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
  
  # if there are no genes within that pathway, ignore it
  # Ideally want to use the input from last time to check but for now we already have a record, so will check at the end
  
  if(nrow(selecting_chromosomes) == 0) {
    message <- paste0("chromosome ", chromosomes_to_analyse, " contains no genes")
    write(message, file = paste0(output_directory,"remove_gene_analysis_empty_chromosomes_info_file_run2.txt"), append = T)
    next()
  }
  
  selecting_chromosomes_BP <- selecting_chromosomes$BP
  
  #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "regular")
  
  # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
  
  GENES_to_snps <- scan(file = paste0(output_directory, Training_name, "_", Validation_name, "_SNPs_", Gene_regions, "_clumped_gene_temp.genes.annot"), what = "", sep = "\n")
  
  y <- strsplit(GENES_to_snps, "[[:space:]]+")
  names(y) <- sapply(y, '[[', 1)
  y <- lapply(y, '[', -1)
  
  y[[1]] <- NULL
  y[[1]] <- NULL
  
  adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
  
  setkey(e$test_data_frame, GENE_NUMBER)
  
  Markers_per_MB_non_independent_normal <- matrix('NA', nrow = nrow(current_table_name), ncol = 6)
  Markers_per_MB_non_independent_normal <- as.data.table(Markers_per_MB_non_independent_normal)
  #Markers_per_MB_non_independent_normal[, names(Markers_per_MB_non_independent_normal) := lapply(.SD, as.character)]
  Genes_to_sort_through_normal <- unlist(current_table_name[,.(Gene)])
  
  for (k in 1:length(Genes_to_sort_through_normal)){
    rows <- nrow(e$test_data_frame[.(Genes_to_sort_through_normal[k])])
    Length_of_gene <- e$test_data_frame[.(Genes_to_sort_through_normal[k]),.(BP_START,BP_END, BP_start_extended, BP_end_extended)][1]
    #Markers_per_MB_non_independent_normal[k,] <- c(k,rows,Length_of_gene)
    testing_2 <- as.matrix(Length_of_gene[1,])
    test <- unlist(c(Genes_to_sort_through_normal[k],rows,Length_of_gene[1,]))
    test <- as.character(test)
    for (j in seq_len(ncol(Markers_per_MB_non_independent_normal))){
      set(Markers_per_MB_non_independent_normal,k,j,test[j])
    }
  }
  
  names(Markers_per_MB_non_independent_normal) <- c("Gene","Nmarkers_in_Gene","BP_START","BP_END","BP_start_extended","BP_end_extended") 
  Markers_per_MB_non_independent_genes_wo_SNPS_normal <- Markers_per_MB_non_independent_normal[is.na(BP_START),.(Gene)]
  Markers_per_MB_non_independent_normal <- Markers_per_MB_non_independent_normal[!is.na(Markers_per_MB_non_independent_normal$BP_START)]
  
  e$test_data_frame[,c("BP_start_extended","BP_end_extended") := NULL]
  which(duplicated(e$test_data_frame$SNP,fromLast = T))
  
  assign("SNPs_for_clumping_normal", unique(e$test_data_frame$SNP), envir = e)
  assign("Gene_regions_annotation_table_normal", e$test_data_frame, envir = e)
  
  
  
  write.table(e$Gene_regions_annotation_table_normal, file = paste0(output_directory, Validation_name, "_", Training_name, "_normal_gene_regions_data_table_after_clumping.txt"), quote = F, row.names = F)
  write(e$SNPs_for_clumping_normal, file = paste0(output_directory, Validation_name, "_", Training_name, "_full_genome_SNPs_after_clumping_normal_gene_regions.txt"))
  write.table(Markers_per_MB_non_independent_normal, file = paste0(output_directory, Validation_name, "_", Training_name, "_normal_gene_region_information_for_randomisation_tests_ater_clumping.txt"), quote = F, row.names = F)
  write.table(Markers_per_MB_non_independent_genes_wo_SNPS_normal, file = paste0(output_directory, Validation_name, "_", Training_name,"_Normal_Genes_without_SNPs_annotated_to_them_full_genome_after_clumping.txt"), quote = F, row.names = F, col.names = F)
}  

if (Gene_regions == "extended") {
  
  Gene_file_name <-  paste0(output_directory, Validation_name, "_", Training_name, "_", Gene_regions,"_gene_regions_Clumped_whole_genome_final.bim")
  selecting_chromosomes <- fread(Gene_file_name)
  names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
  
  # if there are no genes within that pathway, ignore it
  # Ideally want to use the input from last time to check but for now we already have a record, so will check at the end
  
  if(nrow(selecting_chromosomes) == 0) {
    message <- paste0("chromosome ", chromosomes_to_analyse, " contains no genes")
    write(message, file = paste0(output_directory,"remove_gene_analysis_empty_chromosomes_info_file_run2.txt"), append = T)
    next()
  }
  
  selecting_chromosomes_BP <- selecting_chromosomes$BP
  
  #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "regular")
  
  # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
  GENES_to_snps <- scan(file = paste0(output_directory,Training_name, "_", Validation_name, "_SNPs_",Gene_regions,"_clumped_gene_temp.genes.annot"), what = "", sep = "\n")
  
  
  
  y <- strsplit(GENES_to_snps, "[[:space:]]+")
  names(y) <- sapply(y, '[[', 1)
  y <- lapply(y, '[', -1)
  
  
  y[[1]] <- NULL
  y[[1]] <- NULL
  
  adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
  
  setkey(e$test_data_frame, GENE_NUMBER)
  
  Markers_per_MB_non_independent_extended <- matrix('NA', nrow = nrow(current_table_name), ncol = 6)
  Markers_per_MB_non_independent_extended <- as.data.table(Markers_per_MB_non_independent_extended)
  #Markers_per_MB_non_independent_extended[, names(Markers_per_MB_non_independent_extended) := lapply(.SD, as.character)]
  Genes_to_sort_through_extended <- unlist(current_table_name[,.(Gene)])
  
  for (k in 1:length(Genes_to_sort_through_extended)){
    rows <- nrow(e$test_data_frame[.(Genes_to_sort_through_extended[k])])
    Length_of_gene <- e$test_data_frame[.(Genes_to_sort_through_extended[k]),.(BP_START,BP_END, BP_start_extended, BP_end_extended)][1]
    #Markers_per_MB_non_independent_extended[k,] <- c(k,rows,Length_of_gene)
    testing_2 <- as.matrix(Length_of_gene[1,])
    test <- unlist(c(Genes_to_sort_through_extended[k],rows,Length_of_gene[1,]))
    test <- as.character(test)
    for (j in seq_len(ncol(Markers_per_MB_non_independent_extended))){
      set(Markers_per_MB_non_independent_extended,k,j,test[j])
    }
  }
  
  names(Markers_per_MB_non_independent_extended) <- c("Gene","Nmarkers_in_Gene","BP_START","BP_END","BP_start_extended","BP_end_extended") 
  Markers_per_MB_non_independent_genes_wo_SNPS_extended <- Markers_per_MB_non_independent_extended[is.na(BP_start_extended),.(Gene)]
  Markers_per_MB_non_independent_extended <- Markers_per_MB_non_independent_extended[!is.na(Markers_per_MB_non_independent_extended$BP_start_extended)]
  
  e$test_data_frame[,c("BP_START","BP_END") := NULL]
  setnames(e$test_data_frame,old = c("BP_start_extended","BP_end_extended"), new = c("BP_START","BP_END"))
  
  which(duplicated(e$test_data_frame$SNP,fromLast = T))
  
  assign("SNPs_for_clumping_extended", unique(e$test_data_frame$SNP), envir = e)
  assign("Gene_regions_annotation_table_extended", e$test_data_frame, envir = e)
  
  write.table(e$Gene_regions_annotation_table_extended, file = paste0(output_directory, Validation_name, "_", Training_name, "_extended_gene_regions_data_table_after_clumping.txt"), quote = F, row.names = F)
  write(e$SNPs_for_clumping_extended, file = paste0(output_directory, Validation_name, "_", Training_name,"_full_genome_SNPs_after_clumping_extended_gene_regions.txt"))
  write.table(Markers_per_MB_non_independent_extended, file = paste0(output_directory, Validation_name, "_", Training_name, "_extended_gene_region_information_for_randomisation_tests_after_clumping.txt"), quote = F, row.names = F)
  write.table(Markers_per_MB_non_independent_genes_wo_SNPS_extended, file = paste0(output_directory, Validation_name, "_", Training_name, "_Extended_Genes_without_SNPs_annotated_to_them_full_genome_after_clumping.txt"), quote = F,row.names = F,col.names = F)
}

## End of loops, now onto producing gene_specific_polygenic risk scores...    

#End Timer
proc.time() - ptm

quit()

# Currently Defunct code
unread_pathways_one <- fread(paste0(output_directory,"Pathways_analysis_empty_pathways_info_file.txt"))
unread_pathways_two <- fread(paste0(output_directory,"Pathways_analysis_empty_pathways_info_file_run2.txt"))

setnames(unread_pathways_one,c("pathways","chromosome"))
setnames(unread_pathways_two,c("pathways","chromosome"))

combined_pathways <- merge(x=unread_pathways_one,y=unread_pathways_one,by= c("pathways, chromosome"), all = TRUE)

if (ncol(combined_pathways) != 2){
  stop("different pathways used after magma analysis, check empty pathways file")
}
warnings()



