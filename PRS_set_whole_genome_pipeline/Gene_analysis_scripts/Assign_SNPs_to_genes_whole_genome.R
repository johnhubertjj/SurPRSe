#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

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

#Training_name <- "CLOZUK_PGC2noclo"
#Validation_name <- "ALSPAC"
#Validation_full_name_serial <- "ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr"
#Gene_file_name <- "~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr22_consensus_with_CLOZUK_PGC2noclo_flipped_alleles_no_duplicates.bim"
#Gene_directory <- "Genes"
#gene_loc_file_name <-"/Users/johnhubert/Dropbox/Stationary_data/NCBI37.3.gene.loc"
#Gene_regions <- "both"
#chromosomes_to_analyse <- 22

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
Assigning_genes <- function(pathway_input, BP.clumped.SNPs, clumped_SNPs, outputfilename, gene.regions = c("normal", "extended"), chromosome.number = chromosomes_to_analyse){
  GR <- deparse(substitute(pathway_input))
  
  if(gene.regions == "normal"){
    
    test_vector <- rep(0,nrow(clumped_SNPs))
    gene_ID_vector <- rep(0,nrow(clumped_SNPs))
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((pathway_input$BP_START <= BP.clumped.SNPs[i] & pathway_input$BP_START >= BP.clumped.SNPs[i]) == T)
      if (length(tmp_gene) == 1) {
        test_vector[i] <- pathway_input$GENE_NAME[tmp_gene]
        gene_ID_vector[i] <- pathway_input$Gene[tmp_gene]
        BP1 <- pathway_input$BP1[tmp_gene]
        BP2 <- pathway_input$BP2[tmp_gene]
        
      }else{
        test_vector[i] <- NA
        gene_ID_vector[i] <- NA
        BP1[i] <- NA
        BP2[i] <- NA
      }
    }
    ## This is specialised notation for a data.table R object
    ## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
    ## data.table will automatically create a new column because the column Gene does not exist in the original data.table
    
    clumped_SNPs[, Gene_name := test_vector]
    clumped_SNPs[, Gene_ID := gene_ID_vector]
    clumped_SNPs[, BP_1 := BP1]
    clumped_SNPs[, BP_2 := BP2]
    
    ## Which SNPs are not in genes?
    ## .I prints a vector of the rows that match the expression within the square brackets
    ## .() would do a similar thing but keeps the format of the object as a data.table
    index.non.genic.SNPs <- clumped_SNPs[,.I[which(is.na(Gene_name))]]
    
    ## Remove those SNPs from the table
    Gene_clumped_SNPs <- clumped_SNPs[-index.non.genic.SNPs]
    Gene_clumped_SNPs <- Gene_clumped_SNPs[order(BP, decreasing = F)]
    
    assign("Gene_clumped_SNPs",Gene_clumped_SNPs, envir = e)
  }
  
  if(gene.regions == "extended"){
    
    test_vector <- rep(0,nrow(clumped_SNPs))
    gene_ID <- rep(0,nrow(clumped_SNPs))
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((pathway_input$BP_start_extended <= BP.clumped.SNPs[i] & pathway_input$BP_end_extended >= BP.clumped.SNPs[i]) == T)
      
      if (length(tmp_gene) == 1) {
        test_vector[i] <- pathway_input$GENE_NAME[tmp_gene]
        gene_ID[i] <- pathway_input$Gene[tmp_gene]
        BP1[i] <- pathway_input$BP_start_extended[tmp_gene]
        BP2[i] <- pathway_input$BP_end_extended[tmp_gene]
        
      }else{
        test_vector[i] <- NA
        gene_ID[i] <- NA
        BP1[i] <- NA
        BP2[i] <- NA
      }
    }
    
    ## This is specialised notation for a data.table R object
    ## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
    ## data.table will automatically create a new column because the column Gene does not exist in the original data.table
    clumped_SNPs[, GENE_NAME := test_vector]
    clumped_SNPs[, Gene := gene_ID]
    clumped_SNPs[, BP_start_extended := BP1]
    clumped_SNPs[, BP_end_extended := BP2]
    
    ## Which SNPs are not in genes?
    ## .I prints a vector of the rows that match the expression within the square brackets
    ## .() would do a similar thing but keeps the format of the object as a data.table
    index.non.genic.SNPs <- clumped_SNPs[,.I[which(is.na(GENE_NAME))]]
    
    ## Remove those SNPs from the table
    Gene_clumped_SNPs <- clumped_SNPs[-index.non.genic.SNPs]
    Gene_clumped_SNPs <- Gene_clumped_SNPs[order(BP, decreasing = F)]
    

    assign("extended_Gene_clumped_SNPs", Gene_clumped_SNPs, envir = e)
  }
}

### Function which adds duplicate SNPs which happen to be inside other genes
### Currently designed for extended gene regions ONLY
### requires the run of gene.out R script to obtain the specific gene.annot file for each pathway and chromosome
adding_unread_genes <- function(MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number){
  
  Melted_MAGMA_list <- melt(y)
  names(Melted_MAGMA_list) <- c("SNP", "Gene")
  Melted_MAGMA_list$SNP <- as.character(Melted_MAGMA_list$SNP) 
  clumped_SNPs$SNP <- as.character(clumped_SNPs$SNP)
  merged_table_one <- merge(clumped_SNPs,Melted_MAGMA_list,by = "SNP", all = F)
  merged_table_one$Gene <- as.numeric(merged_table_one$Gene)
  Gene_clumped_SNPs <- merge(MAGMA.gene.regions,merged_table_one, by = "Gene",all.y = T)
  Gene_clumped_SNPs[,"CHR.y" :=NULL]  # remove extra_column
  setcolorder(Gene_clumped_SNPs, c("CHR.x","SNP", "GD", "BP", "A1","A2","GENE_NAME","Gene", "BP_START", "BP_END", "BP_start_extended","BP_end_extended","STRAND"))
  setnames(Gene_clumped_SNPs, c("CHR.x","Gene"), c("CHR","GENE_NUMBER"))
  assign("test_data_frame", Gene_clumped_SNPs, envir = e)
  
}

getwd()

write( paste0("This_file_contains_MAGMA_errors\t",0), file = paste0(output_directory,"MAGMA_empty_files_after_analysis.txt"))
write("This file contains genes with no SNPs", file = paste0(output_directory,"remove_gene_analysis_empty_chromosomes_info_file_run2.txt"))


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
    
    if (Gene_regions == "both") {
    
    #### Regular gene regions ####
    #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "regular")
    
      # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
    error_catching <- tryCatch({
      GENES_to_snps <- scan(file = paste0(output_directory,chromosomes_to_analyse,"_",Training_name, "_", Validation_name, "_SNPs_gene_temp.genes.annot"), what = "", sep = "\n")
      }, error=function(e){cat("all genes on chromosome",chromosomes_to_analyse,"must_be_empty,check with MAGMA errors\n")
      message <- paste0(chromosomes_to_analyse)
      write(message, file = paste0(output_directory,"MAGMA_empty_files_after_analysis.txt"), append = T)
      list(next_loop = T)
      })
      
    y <- strsplit(GENES_to_snps, "[[:space:]]+")
    names(y) <- sapply(y, '[[', 1)
    y <- lapply(y, '[', -1)
    
    y[[1]] <- NULL
    y[[1]] <- NULL
    
    adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
    
    e$test_data_frame[,c("BP_start_extended","BP_end_extended") := NULL]
    which(duplicated(e$test_data_frame$SNP,fromLast = T))
    
    assign("SNPs_for_clumping_normal", unique(e$test_data_frame$SNP), envir = e)
    assign("Gene_regions_annotation_table_normal", e$test_data_frame, envir = e)
    
    
    

    #### Extended gene regions ####
    #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "extended")
   
    # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
    error_catching <- tryCatch({
    GENES_to_snps <- scan(file = paste0(output_directory,chromosomes_to_analyse,"_",Training_name, "_", Validation_name, "_SNPs_extended_gene_temp.genes.annot"), what = "", sep = "\n")
    }, error=function(e){cat("all genes on chromosome",chromosomes_to_analyse,"must_be_empty,check with MAGMA errors\n")
      message <- paste0(chromosomes_to_analyse)
      write(message, file = paste0(output_directory,"MAGMA_empty_files_after_analysis.txt"), append = T)
      list(next_loop = T)
    })
    
    
    y <- strsplit(GENES_to_snps, "[[:space:]]+")
    names(y) <- sapply(y, '[[', 1)
    y <- lapply(y, '[', -1)
    
    
    y[[1]] <- NULL
    y[[1]] <- NULL
    
    adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
    
    e$test_data_frame[,c("BP_START","BP_END") := NULL]
    setnames(e$test_data_frame,old = c("BP_start_extended","BP_end_extended"), new = c("BP_START","BP_END"))
    
    which(duplicated(e$test_data_frame$SNP,fromLast = T))
    
    assign("SNPs_for_clumping_extended", unique(e$test_data_frame$SNP), envir = e)
    assign("Gene_regions_annotation_table_extended", e$test_data_frame, envir = e)
    
    
    write.table(e$Gene_regions_annotation_table_normal, file = paste0(output_directory,Validation_name,"_",Training_name,"_chromosome_", chromosomes_to_analyse, "_normal_gene_regions_data_table.txt"), quote = F, row.names = F)
    write.table(e$Gene_regions_annotation_table_extended, file = paste0(output_directory,Validation_name,"_",Training_name,"_chromosome_", chromosomes_to_analyse, "_extended_gene_regions_data_table.txt"), quote = F, row.names = F)
    write(e$SNPs_for_clumping_normal, file = paste0(output_directory,"chromosome_",chromosomes_to_analyse,"_SNPs_for_clumping_normal_gene_regions.txt"))
    write(e$SNPs_for_clumping_extended, file = paste0(output_directory,"chromosome_",chromosomes_to_analyse,"_SNPs_for_clumping_extended_gene_regions.txt"))
    }
  
 
    
    if (Gene_regions == "normal") {
      
      #### Regular gene regions ####
      #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "regular")
      
      # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
      error_catching <- tryCatch({
        GENES_to_snps <- scan(file = paste0(output_directory,chromosomes_to_analyse,"_",Training_name, "_", Validation_name, "_SNPs_gene_temp.genes.annot"), what = "", sep = "\n")
      }, error=function(e){cat("all genes on chromosome",chromosomes_to_analyse,"must_be_empty,check with MAGMA errors\n")
        message <- paste0(chromosomes_to_analyse)
        write(message, file = paste0(output_directory,"MAGMA_empty_files_after_analysis.txt"), append = T)
        list(next_loop = T)
      })
      
      y <- strsplit(GENES_to_snps, "[[:space:]]+")
      names(y) <- sapply(y, '[[', 1)
      y <- lapply(y, '[', -1)
      
      y[[1]] <- NULL
      y[[1]] <- NULL
      
      adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
      
      e$test_data_frame[,c("BP_start_extended","BP_end_extended") := NULL]
      which(duplicated(e$test_data_frame$SNP,fromLast = T))
      
      assign("SNPs_for_clumping_normal", unique(e$test_data_frame$SNP), envir = e)
      assign("Gene_regions_annotation_table_normal", e$test_data_frame, envir = e)
      
      write.table(e$Gene_regions_annotation_table_normal, file = paste0(output_directory,Validation_name,"_",Training_name,"_chromosome_", chromosomes_to_analyse, "_normal_gene_regions_data_table.txt"), quote = F, row.names = F)
      write(e$SNPs_for_clumping_normal, file = paste0(output_directory,"chromosome_",chromosomes_to_analyse,"_SNPs_for_clumping_normal_gene_regions.txt"))
    } 
    
    if (Gene_regions == "extended") {
     
      #### Extended gene regions ####
      #Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, chromosome.number = l, gene.regions = "extended")
      
      # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
      error_catching <- tryCatch({
        GENES_to_snps <- scan(file = paste0(output_directory,chromosomes_to_analyse,"_",Training_name, "_", Validation_name, "_SNPs_extended_gene_temp.genes.annot"), what = "", sep = "\n")
      }, error=function(e){cat("all genes on chromosome",chromosomes_to_analyse,"must_be_empty,check with MAGMA errors\n")
        message <- paste0(chromosomes_to_analyse)
        write(message, file = paste0(output_directory,"MAGMA_empty_files_after_analysis.txt"), append = T)
        list(next_loop = T)
      })
      
      
      y <- strsplit(GENES_to_snps, "[[:space:]]+")
      names(y) <- sapply(y, '[[', 1)
      y <- lapply(y, '[', -1)
      
      
      y[[1]] <- NULL
      y[[1]] <- NULL
      
      adding_unread_genes(MAGMA.gene.regions.for.chromosome = current_table_name, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = chromosomes_to_analyse)  
      
      e$test_data_frame[,c("BP_START","BP_END") := NULL]
      setnames(e$test_data_frame,old = c("BP_start_extended","BP_end_extended"), new = c("BP_START","BP_END"))
      
      which(duplicated(e$test_data_frame$SNP,fromLast = T))
      
      assign("SNPs_for_clumping_extended", unique(e$test_data_frame$SNP), envir = e)
      assign("Gene_regions_annotation_table_extended", e$test_data_frame, envir = e)

      write.table(e$Gene_regions_annotation_table_extended, file = paste0(output_directory,Validation_name,"_",Training_name,"_chromosome_", chromosomes_to_analyse, "_extended_gene_regions_data_table.txt"), quote = F, row.names = F)
      write(e$SNPs_for_clumping_extended, file = paste0(output_directory,"chromosome_",chromosomes_to_analyse,"_SNPs_for_clumping_extended_gene_regions.txt"))
    }
    ## End of loops    

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



