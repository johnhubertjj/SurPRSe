###################################################################################################
#Generate random gene sets from all/brain-expressed genes, where
#probability of selection ~ gene_length, snp_n, indep_snp_n, snp_density, indep_snp_density (didn't use chr)
###################################################################################################

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

tmp = read.table(gene_info_fname,header = TRUE,comment.char = '',sep='\t',as.is = TRUE);

tmp[['LEN']] = (abs(tmp[['STOP']] - tmp[['START']]) + 1)/1000; #gene length (Kb)
tmp[['SNP_DENSITY']] = (1000*tmp[['NSNPS']])/tmp[['LEN']];     #marker density (N per Mb)
tmp[['PARAM_DENSITY']] = (1000*tmp[['NPARAM']])/tmp[['LEN']];  #independent marker density (N per Mb)

gene_data = subset(tmp,select = c(GENE,CHR,LEN,NSNPS,SNP_DENSITY,NPARAM,PARAM_DENSITY));
names(gene_data) = c('entrez_id','chr','len','snp_n','snp_density','indep_snp_n','indep_snp_density');

print(paste('Total number of genes read = ',dim(gene_data)[1]));

###################################################################################################
#read & process gene-set data

for (fname in set_files) {
	
	tmp = read.table(fname,header = FALSE,comment.char = '',sep='\t',as.is = TRUE);
	
	for (set_name in unique(tmp[['V1']])) {
		
		set_genes = unique(subset(tmp,V1 == set_name)[,2]);
		
		set_len[[set_name]] = length(set_genes);
		
		gene_data[[set_name]] = mapply(function(x) {binary_membership(set_genes,x)},gene_data[['entrez_id']]);
	}
	
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
for (next_name in set_names) {
	
	random_mtx = create_random(gene_data,next_name,set_magma_len[[next_name]],formula_str,rand_n);
	
	#save samples to file
    write.table(random_mtx,file = paste(next_name,'.txt',sep=''),col.names = FALSE,row.names = FALSE,sep = '\t',quote = FALSE);
	
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
# 