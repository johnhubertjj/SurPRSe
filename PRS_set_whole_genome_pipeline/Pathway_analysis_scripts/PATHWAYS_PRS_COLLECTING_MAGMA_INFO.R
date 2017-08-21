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
Pathway_directory <- args[6]
Pathway_file_name <- args[7] # the name of the file to be accessed (must be in stationary directory)
gene_loc_file_name <- args[8] # the name of the file containing gene locations (magma's version)
chromosomes_to_analyse <- as.numeric(args[c(9:length(args))])

# Create new variables based on input files to make things easier to read # 
First_half_of_input_genotype_file <- paste0("./", Training_name, "_", Validation_name, "_output/", Pathway_directory, "/", Validation_full_name_serial)
Second_half_of_input_genotype_file <- paste0("_consensus_with_", Training_name, "_flipped_alleles_no_duplicates.bim")
output_directory <- paste0("./", Training_name, "_", Validation_name, "_output/", Pathway_directory, "/")
Date <- Sys.Date()


# Reading_in_Pathway_files
library(data.table)

### environment for functions
e <- new.env()

##### read in pathway sets and standardise column names
pathway_sets <- fread(Pathway_file_name)
setnames(pathway_sets, c("Pathway", "Gene"))

write("This file contains pathways with no SNPs", file = paste0(output_directory,"MAGMA_empty_files_after_analysis.txt"))

# input the name of the pathways from the orignal file
#useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning|memory|conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities|coordination|movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")

#extract names of each pathway
factorise_column_1<- as.factor(pathway_sets$Pathway)
pathway_names <- levels(factorise_column_1)
number_of_pathways_to_analyse <- length(pathway_names)

getwd()

write(pathway_names, file = paste0(output_directory,"Pathway_names.txt"))
# set the key to the pathway name
setkey(pathway_sets, Pathway)

# assign an individual table for each pathway
for (i in 1:number_of_pathways_to_analyse) {
  assign(pathway_names[i], subset(pathway_sets, Pathway == pathway_names[i]), envir = .GlobalEnv)
} 

#read in general MAGMA annotation file
MAGMA.gene.regions <- fread(gene_loc_file_name, colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))


# loop through both chromosomes and pathway tables
for (i in 1:number_of_pathways_to_analyse) {
  
  # merge the original MAGMA annotation file and the pathway file in order to get extra info including BP extensions for regulatory regions 
  assign(paste0("merged",pathway_names[i]), merge(eval(parse(text = paste0("`",pathway_names[i],"`"))), MAGMA.gene.regions, by = "Gene", all = F, sort = F), envir = .GlobalEnv)
  setkey(eval(parse(text = paste0("`","merged",pathway_names[i],"`"))),STRAND) 
  current_table_name <- eval(parse(text = paste0("`","merged",pathway_names[i],"`")))
  
  # equivalent to MAGMAs window option without ignoring strand
  current_table_name <- current_table_name[STRAND == "+",BP_start_extended := BP_START - 35000]
  current_table_name <- current_table_name[STRAND == "+",BP_end_extended := BP_END + 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_start_extended := BP_START - 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_end_extended := BP_END + 35000]
  
  # setkey to the chromosome and remove all SNPs in the X chromosome to save time
  setkey(current_table_name, CHR)
  current_table_name <- current_table_name[!"X"]
  assign(paste0("Gene_regions_all_",pathway_names[i]), current_table_name, envir = .GlobalEnv)
  
  # Now start on the chromosomes
  for (l in chromosomes_to_analyse){
    
    selecting_chromosomes <- fread(paste0(First_half_of_input_genotype_file, l, Second_half_of_input_genotype_file))
    names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
    
    temp_pathway_table <- current_table_name[CHR == l]
    
    if(nrow(temp_pathway_table) == 0) {
      message <- paste0(pathway_names[i],"\t",l)
      write(message, file = paste0(output_directory,"Pathways_analysis_empty_pathways_info_file.txt"), append = T)
      next()
    }
    
    write.table (x = temp_pathway_table[,.(Gene,CHR,BP_START,BP_END,STRAND,GENE_NAME)], file = paste0(output_directory, Training_name,"_", Validation_name,"_", pathway_names[i],"_chromosome_", l, "_temp.gene.loc"), quote = F, col.names = F, row.names = F)
  }
}

# WILL HAVE TO STOP HERE AND REDO ALL OF THIS BELOW BECAUSE YOU NEED TO RUN MAGMA OUTSIDE OF R, PROBLEMS WITH ENVIRONMENT OTHERWISE    

