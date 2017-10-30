###################################################################################################
#Generate random gene sets from all/brain-expressed genes, where
#probability of selection ~ gene_length, snp_n, indep_snp_n, snp_density, indep_snp_density (didn't use chr)
###################################################################################################
library (data.table)
e <- new.env()

## specify the different input tables #
Training_name <- "CLOZUK_PGC2noclo"
Validation_name <- "ALSPAC"
Chromosomes_to_split <- seq(1:22)
print(Chromosomes_to_split)
gene_loc_file_name <- "~/Dropbox/Stationary_data/NCBI37.3.gene.loc"

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

GENES_to_snps <- scan("~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_CLOZUK_PGC2noclo_whole_genome_clumped_annotated_by_gene.genes.annot", what = "", sep = "\n")
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
MAF_counts <- fread("~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_CLOZUK_PGC2noclo_normal_gene_regions_Clumped_whole_genome_final.frq", colClasses = c(rep("character", 4),"numeric","integer"))
output_directory_2 <- "~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/"
bim_file <- fread(paste0(output_directory_2,"ALSPAC_CLOZUK_PGC2noclo_normal_gene_regions_Clumped_whole_genome_final.bim"))
names(bim_file) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")

adding_unread_genes(MAGMA.gene.regions = current_table_name, y=y, clumped_SNPs = bim_file)

setnames(e$Gene_clumped_SNPs,old = "CHR.x", new = "CHR")

Indep_snp_n <- fread("~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_CLOZUK_PGC2noclo_normal_gene_region_information_for_randomisation_tests_ater_clumping.txt", colClasses = rep("numeric",6))
setnames(Indep_snp_n, old = "Gene", new = "GENE_NUMBER")
setnames(Indep_snp_n, old = "Nmarkers_in_Gene", new="Nmarkers_in_Gene_independent")
LD_scored_snps <- fread("~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_CLOZUK_PGC2noclo_normal_gene_regions_Clumped_whole_genome_final.l2.ldscore", colClasses = c("character","character", "integer", "numeric"))
merged1 <- merge(Indep_snp_n,e$Gene_clumped_SNPs,by = c("GENE_NUMBER","BP_START","BP_END","BP_start_extended","BP_end_extended"))
merged2 <- merge(merged1,LD_scored_snps, by = c("CHR","SNP","BP"), all = F)
merged3 <- merge(merged2,MAF_counts[, !"NCHROBS"], by = c("CHR","SNP"), all = T)

#Indep_snp_n <- Markers_per_MB_non_independent_normal
#Indep_snp_n <- Indep_snp_n[, lapply(.SD, as.integer)]





## Read in the latest Summary stats tables after converting to one table
for (i in Chromosomes_to_split){
  assign(paste0("Training_table", i), fread(paste0(output_directory_2, Validation_name,"_", Training_name,"_chromosome_", i, "_normal_gene_region_information_for_randomisation_tests.txt")),envir = .GlobalEnv)
}

l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_split) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}

combined_final_table <- rbindlist(l)
snp_n <- combined_final_table

Pocklington_sets <- fread("~/Dropbox/Stationary_data/Selected_Pocklington_plus_GO_pathways_SCHIZ.txt")
Pathways <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268") 

current_pathway <- fread("~/Documents/ALSPAC_hrc_imputed_bestguess_pathway/CLOZUK_PGC2noclo_ALSPAC_original_output/Pathways/5HT_2C/ALSPAC_original_CLOZUK_PGC2noclo_5HT_2C_Clumped_whole_genome_final.bim")

names(current_pathway) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")

#setwd("/Users/andrew/Projects (active)/2016/SNP methodology/code/random");
#source("generate_random.R");

root_dir = '/Users/andrew/Projects (active)/2016/SNP methodology';

info_dir =  paste(root_dir,'/raw/CLOZUK2',sep='');
gene_info_fname = 'CLOZUK2_PGC2.METAL.genes.txt';

set_files = c('Data Supplement - Significant Sets in MAGMA analysis.txt','GeneSets_LoFi_BrainExpSpc.txt');

out_dir = paste(root_dir,'/working/random/CLOZUK2/',sep='');

#number of random sets to generate for each gene-set
rand_n = 10000;

#length of gene sets
set_len = list();

#formula used to predict gene set membership in regression model
formula_str = 'status ~ len + L2 + MAF';
###################################################################################################
#functions

binary_membership <- function(set,element) {
	
	if (element %in% set) { 1}
	else {0}
	
}

create_random <- function(gene_data,set_name,set_n,formula_str,rand_n) {
	
	gene_data[['status']] = gene_data[[set_name]]; #glm doesn't like variables starting with a number...
	
	membership_model<- glm(as.formula(formula_str),binomial(link = 'logit'),gene_data);
	
	membership_prob = predict(membership_model, type = 'response');
	
	tmp_rand = t(replicate(rand_n,sample(gene_data[['snp']],set_n,replace = FALSE,prob = membership_prob)));
	
	tmp_name = mapply(function(k) paste(set_name,k,sep='_random_'),(1:rand_n));
	
	cbind(tmp_name,tmp_rand);
	
}
###################################################################################################
#read & process gene data
setwd(info_dir);

tmp <- merged3 ;
tmp[['LEN']] = (abs(tmp[['BP_END']] - tmp[['BP_START']]) + 1)/1000; #gene length (Kb)
merged_total <- tmp
names(merged_total) = c('chr','snp','BP','entrez_id','BP_START','BP_END','BP_start_extended','BP_end_extended','Nmarkers_in_Gene_independent','GD','A1','A2','GENE_NAME','STRAND',"L2",'A1.y','A2.y','MAF','len');

print(paste('Total number of SNPs read = ',(dim(merged_total)[1])));

###################################################################################################
## Read & process gene-set data

  tmp <- fread("~/Dropbox/Stationary_data/Selected_Pocklington_plus_GO_pathways_SCHIZ.txt")
  names(tmp) <- c("Pathway","entrez_id")
  tmp_5htc <- tmp[Pathway == "5HT_2C"]
  test_full_file <- merge(gene_data,tmp_5htc, by = "Gene")
  
  
  #Pathways <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268") 
  
	tmp = read.table(fname,header = FALSE,comment.char = '',sep='\t',as.is = TRUE);
	
	for (set_name in unique(tmp[['Pathway']])) {
		
		set_genes = unique(subset(tmp, Pathway == set_name)[,2]);
		
		tmp_pathway <- tmp[Pathway == set_name]
		tmp_binary_measurements <- merge(tmp_pathway,merged_total, by = "entrez_id", all = F)
		tmp_binary_measurements <- as.data.frame(tmp_binary_measurements)  
		  
		set_len[[set_name]] = nrow(tmp_binary_measurements);
		merged_total[[set_name]] = mapply(function(x) {binary_membership(tmp_binary_measurements$snp,x)},merged_total[['snp']]);
	}
	
set_magma_len = lapply(set_names,function(x) sum(merged_total[[x]]));
names(set_magma_len) = set_names;
set_names = names(set_len);


###################################################################################################
# create random gene sets

#start time
time1=Sys.time();

######################	
#all genes
setwd("~/Documents/testing_random_gene_sets2")

for (next_name in set_names) {
	
	random_mtx = create_random(merged_total, next_name, set_magma_len[[next_name]], formula_str, rand_n);
	random_mtx_df <- as.data.frame(random_mtx)
	random_mtx_dt <- as.data.table(random_mtx)
	
	colanames_rndm_df <- colnames(random_mtx_df[,-1])
	test2 <- melt(random_mtx_df,measure.vars =  colanames_rndm_df)
	final_input_list <- test2[order(test2$tmp_name),]
	final_input_list <- subset(final_input_list,  select = c(tmp_name,value))
	assign(paste0(next_name,"_plinkinput"), random_mtx_dt, envir = e)
	#save samples to file
    write.table(final_input_list,file = paste0(next_name,'_random_sets.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE);
    write.table(random_mtx_dt,file = paste0(next_name,'_random_sets_data_table.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE);
    
}
######################
#brain-expressed genes
setwd(paste(out_dir,'brain',sep='/'));
for (next_name in setdiff(set_names,c('fagerberg_brain_expressed'))) {
	
	random_mtx = create_random(brain_expressed_gene_data,next_name,set_magma_len[[next_name]],formula_str,rand_n);
	
	#save samples to file
    write.table(random_mtx,file = paste(next_name,'_brain_expressed.txt',sep=''),col.names = FALSE,row.names = FALSE,sep = '\t',quote = FALSE);
	
}

significance_thresholds <- 1
Summary_stats_dataset <- fread("~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/combined_CLOZUK_PGC2noclo_table_with_CHR.POS_identifiers.txt")


combined.test.training.clumped.Genomic.SNPs <- merge(bim_file, Summary_stats_dataset, by.x="SNP", by.y="SNP", all=F, sort=F)
combined.test.training.clumped.Genomic.SNPs$A1.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A1.y)
combined.test.training.clumped.Genomic.SNPs$A2.y <- toupper(combined.test.training.clumped.Genomic.SNPs$A2.y)


for (next_name in set_names){
  current_pathway <- fread(paste0("~/Documents/testing_random_gene_sets/",next_name,"_random_sets_data_table.txt"),header = F)
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
      
      filename <- paste0('./Scores/',next_name,'_random_',col,'_with_', significance_thresholds[w],".score")
      
      write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
    }
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
      
      
    }
}
}

#end time
time2=Sys.time();
tmp_time_str = paste("time taken =",difftime(time2,time1),units(difftime(time2,time1)),sep =' ');
print(tmp_time_str);


#Total number of genes read =  17619
#Time taken = 19.0934822201729 mins


# Could do via genes instead of SNPS, correct for gene size, and number of genes, (MAF?)
# Make sure to prune the background set of SNPs beforehand, that way we don't have to worry about LD.


set.seed(10)
messy <- data.frame(id = 1:4,
                    trt = sample(rep(c('control', 'treatment'), each = 2)),
                    work.T1 = runif(4),
                    home.T1 = runif(4),
                    work.T2 = runif(4),
                    home.T2 = runif(4))

