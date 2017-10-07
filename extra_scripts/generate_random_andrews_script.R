###################################################################################################
#Generate random gene sets from all/brain-expressed genes, where
#probability of selection ~ gene_length, snp_n, indep_snp_n, snp_density, indep_snp_density (didn't use chr)
###################################################################################################
library (data.table)

Indep_snp_n <- fread("~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/ALSPAC_CLOZUK_PGC2noclo_normal_gene_region_information_for_randomisation_tests_ater_clumping.txt")
setnames(Indep_snp_n,old = "Nmarkers_in_Gene",new="Nmarkers_in_Gene_independent")
output_directory_2 <- "~/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/"

## specify the different input tables #
Training_name <- "CLOZUK_PGC2noclo"
Validation_name <- "ALSPAC"
Chromosomes_to_split <- seq(1:22)
print(Chromosomes_to_split)

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

#setwd("/Users/andrew/Projects (active)/2016/SNP methodology/code/random");
#source("generate_random.R");

root_dir = '/Users/andrew/Projects (active)/2016/SNP methodology';

info_dir =  paste(root_dir,'/raw/CLOZUK2',sep='');
gene_info_fname = 'CLOZUK2_PGC2.METAL.genes.txt';

set_files = c('Data Supplement - Significant Sets in MAGMA analysis.txt','GeneSets_LoFi_BrainExpSpc.txt');

out_dir = paste(root_dir,'/working/random/CLOZUK2/',sep='');

#number of random sets to generate for each gene-set
rand_n = 1000;

#length of gene sets
set_len = list();

#formula used to predict gene set membership in regression model
formula_str = 'status ~ len + snp_n + snp_density + indep_snp_n + indep_snp_density';
###################################################################################################
#functions

binary_membership <- function(set,element) {
	
	if (element %in% set) { 1}
	else {0}
	
}

create_random <- function(gene_data,set_name,set_n,formula_str,rand_n) {
	
	gene_data[['status']] = gene_data[[set_name]]; #glm doesn't like variables starting with a number...
	
	membership_model<- glm(as.formula(formula_str),binomial(link = 'logit'),gene_data);
	
	membership_prob = predict(membership_model,type = 'response');
	
	tmp_rand = t(replicate(rand_n,sample(gene_data[['entrez_id']],set_n,replace = FALSE,prob = membership_prob)));
	
	tmp_name = mapply(function(k) paste(set_name,k,sep='_random_'),(1:rand_n));
	
	cbind(tmp_name,tmp_rand);
	
}
###################################################################################################
#read & process gene data
setwd(info_dir);

tmp <- merge(snp_n,Indep_snp_n,by = c("Gene","BP_START","BP_END","BP_start_extended","BP_end_extended"),all = F) ;
tmp[['LEN']] = (abs(tmp[['BP_END']] - tmp[['BP_START']]) + 1)/1000; #gene length (Kb)
tmp[['SNP_DENSITY']] = (1000*tmp[['Nmarkers_in_Gene']])/tmp[['LEN']];     #marker density (N per Mb)
tmp[['PARAM_DENSITY']] = (1000*tmp[['Nmarkers_in_Gene_independent']])/tmp[['LEN']];  #independent marker density (N per Mb)
tmp2 <- merge(tmp,MAGMA.gene.regions,by=c("Gene","BP_START","BP_END","BP_start_extended","BP_end_extended"))
tmp2 <- as.data.frame(tmp2)
gene_data = subset(tmp2,select = c(Gene,CHR,LEN,Nmarkers_in_Gene,SNP_DENSITY,Nmarkers_in_Gene_independent,PARAM_DENSITY));
names(gene_data) = c('entrez_id','chr','len','snp_n','snp_density','indep_snp_n','indep_snp_density');

print(paste('Total number of genes read = ',dim(gene_data)[1]));

###################################################################################################
#read & process gene-set data

  tmp <- fread("~/Dropbox/Stationary_data/Selected_Pocklington_plus_GO_pathways_SCHIZ.txt")
  tmp <- as.data.frame(tmp)
  
  #Pathways <- c("5HT_2C", "Cav2_channels", "FMRP_targets", "abnormal_behavior", "abnormal_long_term_potentiation", "abnormal_nervous_system_electrophysiology", "Calcium_ion_import_GO0070509", "Membrane_depolarization_during_action_potential_GO0086010", "Synaptic_transmission_GO0007268") 
  
	tmp = read.table(fname,header = FALSE,comment.char = '',sep='\t',as.is = TRUE);
	
	for (set_name in unique(tmp[['V1']])) {
		
		set_genes = unique(subset(tmp,V1 == set_name)[,2]);
		
		set_len[[set_name]] = length(set_genes);
		
		gene_data[[set_name]] = mapply(function(x) {binary_membership(set_genes,x)},gene_data[['entrez_id']]);
	}
	

set_names = names(set_len);

#subset of data including only brain-expressed genes
brain_expressed_gene_data = subset(gene_data,fagerberg_brain_expressed == 1);

set_magma_len = lapply(set_names,function(x) sum(gene_data[[x]]));
names(set_magma_len) = set_names;

#summary table of gene numbers
set_summary = data.frame('set' = set_names,'total' = unlist(set_len),'magma' = unlist(set_magma_len),row.names = NULL);

set_summary[['magma_brain_expressed']] = mapply(function(x) sum(brain_expressed_gene_data[[x]]),set_names);

#save summary data
setwd(out_dir);
write.table(set_summary,file = 'gene_set_summary.txt',col.names = FALSE,row.names = FALSE,sep = '\t',quote = FALSE);
###################################################################################################
# create random gene sets

#start time
time1=Sys.time();

######################	
#all genes
setwd(paste(out_dir,'all',sep='/'));
setwd("~/Documents/testing_random_gene_sets")

for (next_name in set_names) {
	
	random_mtx = create_random(gene_data,next_name,set_magma_len[[next_name]],formula_str,rand_n);
	random_mtx_df <- as.data.frame(random_mtx)
	colanames_rndm_df <- colnames(random_mtx_df[,-1])
	test2 <- melt(random_mtx_df,measure.vars =  colanames_rndm_df)
	final_input_list <- test2[order(test2$tmp_name),]
	final_input_list <- subset(final_input_list,  select = c(tmp_name,value))
	
	#save samples to file
    write.table(final_input_list,file = paste0(next_name,'_random_sets.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE);
	
}
######################
#brain-expressed genes
setwd(paste(out_dir,'brain',sep='/'));
for (next_name in setdiff(set_names,c('fagerberg_brain_expressed'))) {
	
	random_mtx = create_random(brain_expressed_gene_data,next_name,set_magma_len[[next_name]],formula_str,rand_n);
	
	#save samples to file
    write.table(random_mtx,file = paste(next_name,'_brain_expressed.txt',sep=''),col.names = FALSE,row.names = FALSE,sep = '\t',quote = FALSE);
	
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

