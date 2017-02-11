# Reading_in_Pathway_files
library(data.table)

##### Setting_up_multi-platform ####
system_information<-Sys.info()

if (system_information[1] == "Windows") fpath <-  "/Users/JJ/" else fpath <-"/Users/johnhubert/"
#####
pathway_sets <- fread(paste0(fpath,"Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Pocklington2015_134sets_LoFi.txt"))
pathway_sets2 <- fread(paste0(fpath, "Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/GeneWide_BGS_strength.txt"))
useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning|memory|conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities|coordination|movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")
setnames(pathway_sets, c("V1", "V2"), c("Pathway", "Gene"))
setnames(pathway_sets2, c("V1", "V2"), c("Pathway", "Gene"))
pathway_sets <- merge(pathway_sets, pathway_sets2, by = c("Pathway","Gene"), all = T)
setkey(pathway_sets, Pathway)

for (i in 1:length(useful_pathways)) {
  assign(useful_pathways[i], subset(pathway_sets, Pathway == useful_pathways[i]), envir = .GlobalEnv)
} 

MAGMA.gene.regions <- fread("NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))

for (i in 1:length(useful_pathways)) {
  assign(paste0("merged",useful_pathways[i]), merge(eval(parse(text = paste0("`",useful_pathways[i],"`"))), MAGMA.gene.regions, by = "Gene", all = F, sort = F), envir = .GlobalEnv)
  setkey(eval(parse(text = paste0("`","merged",useful_pathways[i],"`"))),STRAND) 
  current_table_name <- eval(parse(text = paste0("`","merged",useful_pathways[i],"`")))
  current_table_name <- current_table_name[STRAND == "+",BP_start_extended := BP_START - 35000]
  current_table_name <- current_table_name[STRAND == "+",BP_end_extended := BP_END + 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_start_extended := BP_START - 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_end_extended := BP_END + 35000]
  setkey(current_table_name, CHR)
  current_table_name <- current_table_name[!"X"]
  assign(paste0("Gene_regions_all_",useful_pathways[i]), current_table_name, envir = .GlobalEnv)
  }

selecting_chromosomes <- fread("/Users/JJ/Dropbox/testing_PRS_chromosome_22/output/CLOZUK_GWAS_BGE_chr22_magma_input.bim")
