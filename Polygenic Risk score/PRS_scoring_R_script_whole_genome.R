##################################################
# Checking location for serial or batch analysis #
##################################################
#! /bin/env/R
# arguments = sex chromosomes involved, serial or batch run, pathway or whole genome?

### Start Timer
ptm <- proc.time()

##set new environment
e <- new.env()

##################################################
# Checking location for serial or batch analysis #
##################################################

System_information <- Sys.info()
whereami <- System_information['user']
args <- commandArgs(trailingOnly = T)

# Arguments 
Training_name <- args[3]
Validation_name <- args[4]
# Check if this works properly # 
sig <- as.numeric(args[5])
assign("Gene_regions", args[6], envir = e)

# running the scoring in parallel function
calculating_scores <- function(i) {
  # Each parallel process represents a p-value threshold
  # So for pval threshold you will have a list of SNPs within the specified gene region (both within the gene and regulatory gene regions)
  
  # This function finds the SNPs within any i gene and writes to a score file containing the SNP Identifier (CHR.POS),
  # the Minor allele and the BETA score as previously computed earlier in the pipeline from an OR
  
  # The parallel process then closes and plink is run in serial to calculate each individual profile
  for (l in 1:length(e$Genes_index)){
    a <- copy(e$combined.CLOZUK.PGC.clumped.Genomic.SNPs)    
    SNPs <- a[, .I[which(Gene_ID == e$Genes_index[l] & P <= e$p.value.thresholds[i])]]
    
    if (length(SNPs) != 0){
      a <- a[SNPs, .(SNP,A1,BETA)]
      
      filename <- paste0('./output/PRS_Scoring/score/whole_Genome_test_', e$Genes_index[l],'_', e$p.value.thresholds[i],".score")
      filename_extended <- paste0('./output/PRS_Scoring/score/whole_Genome_test_', e$Genes_index[l],'_', e$p.value.thresholds[i],"_extended.score")
      
      if(length(a$SNP) == 1){
        one_SNP_file <- paste0("./output/PRS_Scoring/score/whole_Genome_test_", e$p.value.thresholds[i],"_one_SNP_only_dictionary.txt")
        one_SNP_file_extended <- paste0("./output/PRS_Scoring/score/whole_Genome_test_", e$p.value.thresholds[i],"_one_SNP_only_dictionary_extended.txt")
        genes_and_pvalue_one_SNP <- data.frame(e$Genes_index[l], e$p.value.thresholds[i])
        
        if(e$Gene_regions == "normal" | e$Gene_regions == "both"){
          write.table(genes_and_pvalue_one_SNP, file = one_SNP_file, quote = F, append = T, col.names = F, row.names = F)
        
        }else if (e$Gene_regions == "expanded" | e$Gene_regions == "both"){
          write.table(genes_and_pvalue_one_SNP, file = one_SNP_file_extended, quote = F, append = T, col.names = F, row.names = F)
    
        }
      }
      
      if(e$Gene_regions == "normal" | e$Gene_regions == "both"){
        write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
        genes_and_pvalue <- data.frame(e$Genes_index[l], e$p.value.thresholds[i])
        write.table(genes_and_pvalue, file = "./output/PRS_scoring/score/Index_of_genes_and_pval_1.txt", quote = F, col.names = F, row.names = F, append = T)
        
      }else if (e$Gene_regions == "expanded" | e$Gene_regions == "both"){
        write.table(file = filename_extended, a, row.names = F, col.names = F, quote = F, sep="\t")        
        genes_and_pvalue <- data.frame(e$Genes_index[l], e$p.value.thresholds[i])
        write.table(genes_and_pvalue, file = "./output/PRS_scoring/score/Index_of_genes_and_pval_1_extended.txt", quote = F, col.names = F, row.names = F, append = T)
      }
      
      rm(a)
      
    }else{
      next()
    }
  }
}

### Library
library(data.table)
library(parallel)
library(base)

### Only if you have specified that you want the normal gene regions in the analysis ###
if (Gene_regions == "normal" | Gene_regions == "both"){

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

### set working directory
setwd(".")

## Read in PGC.data and input data
PGC <- fread(paste0("./output/", Training_name, "_table_for_python.txt"))
test_data_frame <- fread(paste0("./output/MAGMA_Gene_regions_for_",Validation_name, "_", Training_name,".txt"))
#test_data_frame <- fread("./output/MAGMA_Gene_regions_for_python_script.txt")
setnames(test_data_frame, c("CHR", "SNP", "BP", "Gene_ID", "BP_1", "BP_2"))



## Excess info if required ##
#############################
    ## Plink --gene-report tmp1.txt Chromosom_MAGMA_GENE_REGIONS_NCBI37.3.txt
##       --out CLOZUK_BGE_GENE_REGIONS

# Plink --bfile CLOZUK_GWAS_BGE_chr22
#       --make-set Chromosome_22_MAGMA_GENE_REGIONS_NCBI37.3.txt
#       --out testlol
#       --write-set
#############################

# Ideas for getting scores and keeping them within the R environment

# 1 use MAGMA and create the gene.annot file: read as a list into R 
# then create a Data.table with repeated SNP's for each gene in the format SNP, A1, BETA, GENE
# Then parse through each row using the data table notation to create a score and profile wihtin an index
# Problem here will be the output to plink
# Can maybe put into a list of data tables for the scoring? Find out if it is possible...if not then just input, print and delete
# 

# Scoring per gene: rough and undescriptive at the moment

# SCORING 
# sig <- c(0.0001, 0.001,0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)  
assign("p.value.thresholds", sig, envir = e)

## check that it is merging properly here (after analysis is run)
combined.CLOZUK.PGC.clumped.Genomic.SNPs <- merge(test_data_frame, PGC, by.x="SNP", by.y="SNP", all=F, sort=F)
combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1)
combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2)
Genes_index <- unique(combined.CLOZUK.PGC.clumped.Genomic.SNPs$Gene_ID)

assign("Genes_index", Genes_index, envir = e)

Genes_used <- matrix(ncol = 2)
assign("combined.CLOZUK.PGC.clumped.Genomic.SNPs", combined.CLOZUK.PGC.clumped.Genomic.SNPs, envir = e)

clusterExport(cl, "e")
parLapply(cl, 1:length(e$p.value.thresholds), calculating_scores)
stopCluster(cl)
  
Genes_index_for_plink <- fread("./output/Index_of_genes_and_pval_1.txt")
setnames(Genes_index_for_plink,c("Genes", "pval"))
setkey(Genes_index_for_plink,pval)

for (i in 1:length(p.value.thresholds)) {
  Genes_index_for_plink_p_val_thresh <- Genes_index_for_plink[pval == e$p.value.thresholds[i]]
  Genes_index_for_plink_p_val_thresh <- Genes_index_for_plink_p_val_thresh$Genes
  
  if(length(Genes_index_for_plink_p_val_thresh == 0)){
    next()
    
  }else{
    for (l in 1:length(Genes_index_for_plink_p_val_thresh)) {
      filename <- paste0('./output/PRS_Scoring/score/whole_Genome_test_', Genes_index_for_plink_p_val_thresh[l],'_', e$p.value.thresholds[i],".score")
      command <- paste0('plink --bfile ./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS --score ', filename, " --out /output/PRS_scoring/Profiles/whole_Genome_test_", Genes_index_for_plink_p_val_thresh[l],'_', e$p.value.thresholds[i], "_a")
      system(command)
    } 
  } 
}
}

### Only if you wish to include the regulatory regions within the analysis as well ###
if(Gene_regions == "expanded" | Gene_regions == "both"){
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type = "FORK")
  
  ### set working directory
  setwd(".")
  
  ## Read in PGC.data and input data
  PGC <- fread(paste0("./output/", Training_name,"_table_for_python.txt"))
  test_data_frame <- fread(paste0("./output/MAGMA_Gene_regions_for_", Validation_name, "_", Training_name,"_extended.txt"))
  #test_data_frame <- fread("./output/MAGMA_Gene_regions_for_python_script.txt")
  setnames(test_data_frame, c("CHR", "SNP", "BP", "Gene_ID", "BP_1", "BP_2"))
  
  ##set new environment
  e <- new.env()
  
  # SCORING 
  # sig <- c(0.0001, 0.001,0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)  
  assign("p.value.thresholds", sig, envir = e)
  
  ## check that it is merging properly here (after analysis is run)
  combined.CLOZUK.PGC.clumped.Genomic.SNPs <- merge(test_data_frame, PGC, by.x="SNP", by.y="SNP", all=F, sort=F)
  combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1)
  combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2)
  Genes_index <- unique(combined.CLOZUK.PGC.clumped.Genomic.SNPs$Gene_ID)
  
  assign("Genes_index", Genes_index, envir = e)
  
  Genes_used <- matrix(ncol = 2)
  assign("combined.CLOZUK.PGC.clumped.Genomic.SNPs", combined.CLOZUK.PGC.clumped.Genomic.SNPs, envir = e)
  
  clusterExport(cl, "e")
  parLapply(cl, 1:length(e$p.value.thresholds), calculating_scores)
  stopCluster(cl)
  
  Genes_index_for_plink <- fread("./output/Index_of_genes_and_pval_1_extended.txt")
  setnames(Genes_index_for_plink,c("Genes", "pval"))
  setkey(Genes_index_for_plink,pval)
  
  for (i in 1:length(e$p.value.thresholds)) {
    Genes_index_for_plink_p_val_thresh <- Genes_index_for_plink[pval == e$p.value.thresholds[i]]
    Genes_index_for_plink_p_val_thresh <- Genes_index_for_plink_p_val_thresh$Genes
    
    if(length(Genes_index_for_plink_p_val_thresh == 0)){
      next()
      filename <- paste0('./output/PRS_Scoring/score/whole_Genome_test_', e$Genes_index[l],'_', e$p.value.thresholds[i],".score")
      filename_extended <- paste0('./output/PRS_Scoring/score/whole_Genome_test_', e$Genes_index[l],'_', e$p.value.thresholds[i],"_extended.score")
    }else{
      for (l in 1:length(Genes_index_for_plink_p_val_thresh)) {
        filename_extended <- paste0('./output/PRS_Scoring/score/whole_Genome_test_', Genes_index_for_plink_p_val_thresh[l],'_', e$p.value.thresholds[i],"_extended.score")
        
        command <- paste0('plink --bfile ./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS --score ', filename_extended, " --out /whole_genome_testing/output/PRS_scoring/Profiles/whole_Genome_test_", Genes_index_for_plink_p_val_thresh[l],'_', e$p.value.thresholds[i], "_a_extended")
        system(command)
      } 
    } 
  }
  
}
#score
## END script!
#######################################
