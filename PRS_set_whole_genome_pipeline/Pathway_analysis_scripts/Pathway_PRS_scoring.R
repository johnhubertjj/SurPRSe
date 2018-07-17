#####################################
###########PRS PER PATHWAY###########
#####################################
# source(paste0(path_to_PRS_scripts,'Pathway_analysis_scripts/Pathway_PRS_scoring.R'))

### Start Timer
ptm <- proc.time()

### Library
library(data.table)
library(parallel)
library(base)

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
Pathway_directory <- args[5]
Gene_output_directory <- args[6]
Pathway_file_name <- args[7] # the name of the file to be accessed (must be in stationary directory)
Gene_regions <- args[8]

# This will be a problem, write in a script that defines what the arguments are (probably have to be linked with the arguments script however...stop un-needed analysis)
significance_thresholds <- as.numeric(args[c(9:length(args))])
print(significance_thresholds)

#Training_name <- "HG19_pgc.scz.full.2012-04"
#Validation_name <- "ALSPAC"
#Pathway_directory <- paste0("./", Training_name, "_", Validation_name,"_output/Pathways/")
#Pathway_file_name <- "/Users/johnhubert/Dropbox/Stationary_data/Pocklington2015_134sets_LoFi.txt"
#significance_thresholds <- c(0.05, 0.5)


##### read in pathway sets and standardise column names
pathway_sets <- fread(Pathway_file_name)
setnames(pathway_sets, c("Pathway", "Gene"))

#extract names of each pathway
factorise_column_1<- as.factor(pathway_sets$Pathway)
pathway_names <- levels(factorise_column_1)
number_of_pathways_to_analyse <- length(pathway_names)

#pathway_names <- pathway_names[1:3]
#number_of_pathways_to_analyse <- 3

# Read in summary_stats_data
Summary_stats_full_dataset <- fread(paste0("./", Training_name, "_", Validation_name,"_output/combined_", Training_name, "_table_with_CHR.POS_identifiers.txt"))

calculating_scores <- function(i,Summary_stats_full_dataset, Validation_name, Training_name, Pathway_directory, pathway_names, Gene_regions){
  
  if (Gene_regions == "both" | Gene_regions == "normal"){
    
    extra_gene_regions <- "normal"
    
    # Read in pathway_bim_file
    test_pathway <- fread(paste0(Pathway_directory, pathway_names[i],"/",Validation_name,"_",Training_name,"_",pathway_names[i],"_",extra_gene_regions,"_Clumped_whole_genome_final.bim"))
    
    
    #for (i in 1:length(Useful_pathways)) {
    #  Maindir <- paste0("~/Documents/testing_PRS_chromosome_22/test_chr5/output/", Useful_pathways[i])
    #  scoredir <- "score"
    #  profiledir <- "Profile"
    #  if (file.exists(scoredir) == FALSE){
    #    dir.create(file.path(Maindir, scoredir))
    #  }
    #  if (file.exists(profiledir) == FALSE){
    #    dir.create(file.path(Maindir, profiledir))
    #  }
    #}
    
    scoring_output_file <- paste0(Pathway_directory,pathway_names[i],"/scoring_",Training_name,"_",Validation_name,"_pathway_",pathway_names[i],"_",extra_gene_regions)
    
    #path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
    
    setnames(test_pathway,c("CHR","SNP","GD","BP","A1","A2"))  
    
    combined.test.training.clumped.Genomic.SNPs <- merge(test_pathway, Summary_stats_full_dataset, by.x="SNP", by.y="SNP", all=F, sort=F)
    combined.test.training.clumped.Genomic.SNPs$A1.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A1.y)
    combined.test.training.clumped.Genomic.SNPs$A2.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A2.y)
    
    ## check that it is merging properly here (after analysis is run)
    
    for (w in 1:length(significance_thresholds)) {
      
      a <- copy(combined.test.training.clumped.Genomic.SNPs)    
      SNPs <- a[, .I[which(P <= significance_thresholds[w])]]
      
      if (length(SNPs) != 0){
        a <- a[SNPs, .(SNP, A1.y, BETA)]
        
        filename <- paste0(scoring_output_file,'_with_', significance_thresholds[w],".score")
        
        write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
        
        write(x = paste(pathway_names[i],significance_thresholds[w],TRUE,sep = " "), file = paste0(Pathway_directory,"thresholds_to_use_",extra_gene_regions,".txt"), append = T)
        
        #      if (Useful_pathways[i] == "abnormal_learning|memory|conditioning" | Useful_pathways[i] == "abnormal_motor_capabilities|coordination|movement") {
        #        Useful_pathways[i] <- "abnormal_learning\\|memory\\|conditioning"
        #        filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],".score")
        #        path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
        #        }
        
        #      if (Useful_pathways[i] == "abnormal_motor_capabilities|coordination|movement") {
        #        Useful_pathways[i] <-"abnormal_learning_motor_capabilities\|coordination\|movement"
        #        filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],".score")
        #        path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
        #        
        #      }
        
        
      }else{
        #cat("No SNPs found for pathway ",pathway_names[i]," no score file produced for ", extra_gene_regions, " gene regions")
        write(x = paste(pathway_names[i],significance_thresholds[w],FALSE,sep = " "), file = paste0(Pathway_directory,"thresholds_to_use_",extra_gene_regions,".txt"), append = T)
        next()
      }
    }
  }
  
  if (Gene_regions == "both" | Gene_regions == "extended"){
    
    extra_gene_regions <- "extended"
    
    # Read in pathway_bim_file
    test_pathway <- fread(paste0(Pathway_directory, pathway_names[i],"/",Validation_name,"_",Training_name,"_",pathway_names[i],"_",extra_gene_regions,"_Clumped_whole_genome_final.bim"))
    
    
    #for (i in 1:length(Useful_pathways)) {
    #  Maindir <- paste0("~/Documents/testing_PRS_chromosome_22/test_chr5/output/", Useful_pathways[i])
    #  scoredir <- "score"
    #  profiledir <- "Profile"
    #  if (file.exists(scoredir) == FALSE){
    #    dir.create(file.path(Maindir, scoredir))
    #  }
    #  if (file.exists(profiledir) == FALSE){
    #    dir.create(file.path(Maindir, profiledir))
    #  }
    #}
    
    scoring_output_file <- paste0(Pathway_directory,pathway_names[i],"/scoring_",Training_name,"_",Validation_name,"_pathway_",pathway_names[i],"_",extra_gene_regions)
    
    #path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
    
    setnames(test_pathway,c("CHR","SNP","GD","BP","A1","A2"))  
    
    combined.test.training.clumped.Genomic.SNPs <- merge(test_pathway, Summary_stats_full_dataset, by.x="SNP", by.y="SNP", all=F, sort=F)
    combined.test.training.clumped.Genomic.SNPs$A1.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A1.y)
    combined.test.training.clumped.Genomic.SNPs$A2.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A2.y)
    
    ## check that it is merging properly here (after analysis is run)
    
    for (w in 1:length(significance_thresholds)) {
      
      a <- copy(combined.test.training.clumped.Genomic.SNPs)    
      SNPs <- a[, .I[which(P <= significance_thresholds[w])]]
      
      if (length(SNPs) != 0){
        a <- a[SNPs, .(SNP, A1.y, BETA)]
        
        filename_score_file <- paste0(scoring_output_file,'_with_', significance_thresholds[w],".score")
        
        write(x = paste(pathway_names[i],significance_thresholds[w],TRUE,sep = " "), file = paste0(Pathway_directory,"thresholds_to_use_",extra_gene_regions,".txt"), append = T)
        
        
        write.table(file = filename_score_file, a, row.names = F, col.names = F, quote = F, sep="\t")
        
        #      if (Useful_pathways[i] == "abnormal_learning|memory|conditioning" | Useful_pathways[i] == "abnormal_motor_capabilities|coordination|movement") {
        #        Useful_pathways[i] <- "abnormal_learning\\|memory\\|conditioning"
        #        filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],".score")
        #        path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
        #        }
        
        #      if (Useful_pathways[i] == "abnormal_motor_capabilities|coordination|movement") {
        #        Useful_pathways[i] <-"abnormal_learning_motor_capabilities\|coordination\|movement"
        #        filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],".score")
        #        path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
        #        
        #      }
        
        
      }else{
        #cat("No SNPs found for pathway ",pathway_names[i]," no score file produced for ", extra_gene_regions, " gene regions")
        write(x = paste(pathway_names[i],significance_thresholds[w],FALSE,sep = " "), file = paste0(Pathway_directory,"thresholds_to_use_",extra_gene_regions,".txt"), append = T)
        next()
      }
    }
  }
  # Bracket for function call
}

if (Gene_regions == "both" | Gene_regions == "normal"){
  
if(file.exists(paste0(Pathway_directory,"thresholds_to_use_normal.txt"))){
  file.remove(paste0(Pathway_directory,"thresholds_to_use_normal.txt"))
}
}

if (Gene_regions == "both" | Gene_regions == "extended"){
  if(file.exists(paste0(Pathway_directory,"thresholds_to_use_extended.txt"))){
    file.remove(paste0(Pathway_directory,"thresholds_to_use_extended.txt"))
  }
}


# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

# Export the environment to the cluster
clusterExport(cl, "e")
parLapply(cl, 1:number_of_pathways_to_analyse, calculating_scores, Summary_stats_full_dataset, Validation_name, Training_name, Pathway_directory, pathway_names, Gene_regions)
stopCluster(cl)