# Reading_in_Pathway_files
library(data.table)

pathway_sets <- fread("~/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Pocklington2015_134sets_LoFi.txt")
useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning|memory|conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities|coordination|movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90")
setnames(pathway_sets, c("V1", "V2"), c("Pathway", "Gene"))
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

