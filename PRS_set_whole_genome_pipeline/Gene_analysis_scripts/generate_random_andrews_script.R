###################################################################################################
# Generate random SNP sets from Input genes, where
# Probability of selection ~ gene_length, LD_score, MAF

# Generate score files from randomised plink sets based on significance thresholds 
# (reads in as plink format but not yet suited to it's file format)
###################################################################################################

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

# Set independent environment when needed #
e <- new.env()

# specify the different input tables #
Training_name <- args[3]
Validation_name <- args[4]
Pathway_output_directory <- args[5]
Gene_output_directory <- args[6]
gene_loc_file_name <- args[7] # The name of the file containing gene locations (magma's version)
Gene_regions <- args[8] # Whether to include/exclude the regulatory regions of a gene
rand_n = args[9]; # Number of random sets to generate for each gene-set
Pathway_file_name <- args[10] # The input file annotating genes to pathways
calculate_indep_SNPs <- args[11]
Randomised_output_directory <- args[12]
Random_scoring_directory <- args[13]
sample_replace <- args[14]
Pathways <- as.character(args[c(15:length(args))]) # All pathways

# Testing_variables
# Training_name <- "CLOZUK_PGC2noclo"
# Validation_name <- "ALSPAC"
# Gene_output_directory <- paste0(Training_name,"_",Validation_name,"_output/Genes/")
# Pathway_output_directory <- paste0(Training_name,"_",Validation_name,"_output/Pathways/")
# gene_loc_file_name <- "/home/c1020109/Stationary_data/NCBI37.3.gene.loc"
# Gene_regions <- "both"
# rand_n = 10000; # Number of random sets to generate for each gene-set
# Pathway_file_name <- "/home/c1020109/Stationary_data/Selected_Pocklington_plus_GO_pathways_SCHIZ.txt"
# Pathways <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268") 

# Read in array variables from file and rearrage to a vector#
significance_thresholds <- fread(paste0(Training_name,"_", Validation_name,"_plink_significance_thresholds_arguments_file_tmp.txt"))      
significance_thresholds <- unlist(significance_thresholds$V3)
print(significance_thresholds)


Chromosomes_to_split <- scan(file = paste0(Training_name,"_", Validation_name,"_chromosomes_to_analyse_arguments_file_tmp.txt"), sep = "\n")
print(Chromosomes_to_split)

# Specify normal or extended gene regions 

# NEEDS SEPARATE LOOPS HERE

# Read in files containing confounders #
MAF_counts <- fread(paste0(Gene_output_directory,Validation_name,"_",Training_name,"_normal_gene_regions_Clumped_whole_genome_final.frq"), colClasses = c(rep("character", 4),"numeric","integer"))
GENES_to_snps <- scan(paste0(Gene_output_directory,Training_name,"_",Validation_name,"_SNPs_normal_clumped_gene_temp.genes.annot"), what = "", sep = "\n")

if (calculate_indep_SNPs == TRUE){ 
Indep_snp_n <- fread(paste0(Gene_output_directory,Validation_name,"_",Training_name,"_normal_gene_region_information_for_randomisation_tests_ater_clumping.txt"), colClasses = rep("numeric",6))
setnames(Indep_snp_n, old = "Gene", new = "GENE_NUMBER")
setnames(Indep_snp_n, old = "Nmarkers_in_Gene", new="Nmarkers_in_Gene_independent")
}

LD_scored_snps <- fread(paste0(Gene_output_directory,Validation_name,"_",Training_name,"_normal_gene_regions_Clumped_whole_genome_final.l2.ldscore"), colClasses = c("character","character", "integer", "numeric"))
Pathway_sets <- fread(Pathway_file_name)

# Read in bim file FROM GENES DIRECTORY
# This is the clumped bim file after limiting down to gene regions (either extended or normal)
bim_file <- fread(paste0(Gene_output_directory, Validation_name, "_", Training_name, "_normal_gene_regions_Clumped_whole_genome_final.bim"))
names(bim_file) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")

# Clarify output directory and where the Summary stats file is #
Randomised_output_directory <- paste0(Pathway_output_directory,"Randomised_gene_sets_analysis/")
Summary_stats_dataset <- fread(paste0("./", Training_name, "_", Validation_name, "_output/combined_",Training_name,"_table_with_CHR.POS_identifiers.txt"))

print(Chromosomes_to_split)

###################################################################################################
#functions
adding_unread_genes <- function(MAGMA.gene.regions, clumped_SNPs, y, chromosome.number){
  
  Melted_MAGMA_list <- melt(y)
  names(Melted_MAGMA_list) <- c("SNP", "Gene")
  Melted_MAGMA_list$SNP <- as.character(Melted_MAGMA_list$SNP) 
  clumped_SNPs$SNP <- as.character(clumped_SNPs$SNP)
  merged_table_one <- merge(clumped_SNPs,Melted_MAGMA_list,by = "SNP", all = F)
  merged_table_one$Gene <- as.numeric(merged_table_one$Gene)
  Gene_clumped_SNPs <- merge(MAGMA.gene.regions,merged_table_one, by = "Gene", all.y = T)
  Gene_clumped_SNPs[,"CHR.y" :=NULL]  # remove extra_column
  setcolorder(Gene_clumped_SNPs, c("CHR.x","SNP", "GD", "BP", "A1","A2","GENE_NAME","Gene", "BP_START", "BP_END","BP_start_extended","BP_end_extended","STRAND"))
  setnames(Gene_clumped_SNPs, c("CHR.x","Gene"), c("CHR","GENE_NUMBER"))
  assign("Gene_clumped_SNPs", Gene_clumped_SNPs, envir = e)
  
}

#formula used to predict gene set membership in regression model
formula_str = 'status ~ len + L2 + MAF';

binary_membership <- function(set,element) {
  
  if (element %in% set) { 1}
  else {0}
  
}

create_random <- function(i,gene_data,set_names,set_n,formula_str,rand_n, sample_replace, Randomised_output_directory) {
 
  set_name <- set_names[i]
  set_n <- set_n[[i]]
  
  gene_data[['status']] = gene_data[[set_name]]; #glm doesn't like variables starting with a number...
  
  membership_model<- glm(as.formula(formula_str),binomial(link = 'logit'),gene_data);
  
  membership_prob = predict(membership_model, type = 'response');
 
  if (sample_replace == "TRUE"){ 
  tmp_rand = t(replicate(rand_n,sample(gene_data[['snp']],set_n,replace = TRUE, prob = membership_prob)));
  }else{
 tmp_rand = t(replicate(rand_n,sample(gene_data[['snp']],set_n,replace = FALSE, prob = membership_prob)));
}
  tmp_name = mapply(function(k) paste(set_name,k,sep='_random_'),(1:rand_n));
  
  random_mtx <-cbind(tmp_name,tmp_rand);
  random_mtx_dt <- as.data.table(random_mtx)
  
  rm(random_mtx)
  testa <- copy(random_mtx_dt)
  colnames_rndm_dt <- colnames(testa[,tmp_name := NULL])
  rm(testa)
  test2 <- melt(random_mtx_dt,measure.vars =  colnames_rndm_dt)
  final_input_list <- test2[order(test2$tmp_name),]
  final_input_list <- subset(final_input_list,  select = c(tmp_name,value))
  
  #save samples to file
  write.table(final_input_list,file = paste0(Randomised_output_directory,set_name,'_random_sets.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE);
  write.table(random_mtx_dt,file = paste0(Randomised_output_directory,set_name,'_random_sets_data_table.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE);
  
}

Calculating_scores_random <- function(i, set_names, Randomised_output_directory, significance_thresholds, rand_n, combined.test.training.clumped.Genomic.SNPs){
  
  set_name <- set_names[i]
  current_pathway <- fread(paste0(Randomised_output_directory,set_name,'_random_sets_data_table.txt'),header = F)
  current_pathway <- t(current_pathway)
  colnames(current_pathway) <- current_pathway[1,]
  current_pathway <- as.data.table(current_pathway)
  current_pathway <- current_pathway[-1]
  
  for (w in 1:length(significance_thresholds)) { 
    for (col in 1:rand_n){
      current_pathway_random_set <- current_pathway[,col, with = F]
      setnames(current_pathway_random_set,"SNP")
      a <- merge(current_pathway_random_set,combined.test.training.clumped.Genomic.SNPs, by = "SNP",all.y = F,all.x = T, sort = F)
      
      SNPs <- a[, .I[which(P <= significance_thresholds[w])]]
      
      if (length(SNPs) != 0){
        a <- a[SNPs, .(SNP, A1.y, BETA)]
        setkey(a,"SNP")
        a <- unique(a)
        
        filename <- paste0(Randomised_output_directory,'Scores/', set_name,'_random_',col,'_with_', significance_thresholds[w], ".score")
        write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
      }
      
    }
  }
  
}
###################################################################################################

## Convert from MAGMA structure to an object easily analysed in R ##
y <- strsplit(GENES_to_snps, "[[:space:]]+")
names(y) <- sapply(y, '[[', 1)
y <- lapply(y, '[', -1)

y[[1]] <- NULL
y[[1]] <- NULL

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

adding_unread_genes(MAGMA.gene.regions = current_table_name, y=y, clumped_SNPs = bim_file)

if (calculate_indep_SNPs == "TRUE"){
# Combine all datasets together with different parameters
merged1 <- merge(Indep_snp_n,e$Gene_clumped_SNPs,by = c("GENE_NUMBER","BP_START","BP_END","BP_start_extended","BP_end_extended"))
merged2 <- merge(merged1,LD_scored_snps, by = c("CHR","SNP","BP"), all = F)
merged3 <- merge(merged2,MAF_counts[, NCHROBS:= NULL], by = c("CHR","SNP"), all = T)

} else {
merged2 <- merge(e$Gene_clumped_SNPs,LD_scored_snps, by = c("CHR","SNP","BP"), all = F)
merged3 <- merge(merged2,MAF_counts[, NCHROBS:= NULL], by = c("CHR","SNP"), all = T)
}


#Indep_snp_n <- Markers_per_MB_non_independent_normal
#Indep_snp_n <- Indep_snp_n[, lapply(.SD, as.integer)]

## Read in the latest Summary stats tables after converting to one table
for (i in Chromosomes_to_split){
  assign(paste0("Training_table", i), fread(paste0(Gene_output_directory, Validation_name,"_", Training_name,"_chromosome_", i, "_normal_gene_region_information_for_randomisation_tests.txt")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_split) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}

# combine genic information on training sets into one table
combined_final_table <- rbindlist(l)
snp_n <- combined_final_table

#length of gene sets
set_len = list();

tmp <- merged3 ;
tmp[['LEN']] = (abs(tmp[['BP_END']] - tmp[['BP_START']]) + 1)/1000; #gene length (Kb)
merged_total <- tmp

if (calculate_indep_SNPs == "TRUE"){
names(merged_total) = c('chr','snp','BP','entrez_id','BP_START','BP_END','BP_start_extended','BP_end_extended','Nmarkers_in_Gene_independent','GD','A1','A2','GENE_NAME','STRAND',"L2",'A1.y','A2.y','MAF','len');

}else{
names(merged_total) = c('chr','snp','BP','GD','A1','A2','GENE_NAME','entrez_id','BP_START','BP_END','BP_start_extended','BP_end_extended','STRAND','L2','A1.y','A2.y','MAF','len');

}
print(paste('Total number of SNPs read = ',(dim(merged_total)[1])));

###################################################################################################
## Read & process gene-set data

  tmp <- Pathway_sets
  names(tmp) <- c("Pathway","entrez_id")
  
	
	for (set_name in unique(tmp[['Pathway']])) {
		
		set_genes = unique(subset(tmp, Pathway == set_name)[,2]);
		
		tmp_pathway <- tmp[Pathway == set_name]
		tmp_binary_measurements <- merge(tmp_pathway,merged_total, by = "entrez_id", all = F)
		tmp_binary_measurements <- as.data.frame(tmp_binary_measurements)  
		  
		set_len[[set_name]] = nrow(tmp_binary_measurements);
		merged_total[[set_name]] = mapply(function(x) {binary_membership(tmp_binary_measurements$snp,x)},merged_total[['snp']]);
	}

set_names = names(set_len);
set_magma_len = lapply(set_names,function(x) sum(merged_total[[x]]));
names(set_magma_len) = set_names;
number_of_pathways <- length(set_names)
write(set_names, file = paste0(Pathway_output_directory,Training_name, "_", Validation_name,"_random_pathways_to_test.txt"),ncolumns = 1)
###################################################################################################
# create random gene sets
###################################################################################################	
# all genes

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

# Export the environment to the cluster
clusterExport(cl, "e")

# Create n random sets at a time where n=number of cores on a singular node
parLapply(cl, 1:number_of_pathways, create_random, merged_total, set_names, set_magma_len, formula_str, rand_n, sample_replace, Randomised_output_directory)
stopCluster(cl)

# significance_thresholds <- 1

combined.test.training.clumped.Genomic.SNPs <- merge(bim_file, Summary_stats_dataset, by.x="SNP", by.y="SNP", all=F, sort=F)
combined.test.training.clumped.Genomic.SNPs$A1.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A1.y)
combined.test.training.clumped.Genomic.SNPs$A2.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A2.y)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

# Export the environment to the cluster
clusterExport(cl, "e")

# Create n random sets at a time where n=number of cores on a singular node
parLapply(cl, 1:number_of_pathways, Calculating_scores_random, set_names, Randomised_output_directory, significance_thresholds, rand_n, combined.test.training.clumped.Genomic.SNPs)
stopCluster(cl)

#End Timer
proc.time() - ptm

# Could do via genes instead of SNPS, correct for gene size, and number of genes, (MAF?)
# Make sure to prune the background set of SNPs beforehand, that way we don't have to worry about LD.



