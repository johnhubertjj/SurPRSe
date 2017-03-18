#! /bin/bash

sig=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
length_PGC=(521 1435 4787 13603 23369 40794 57657 74106 89961)

sig{[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue.txt
for i in `seq 0 8` ;
do 
	
	awk '{print $2}' ${sig[i]}CLOZUK_PGC_P_Vals_reference_file.txt > ${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue.txt
	magma --annotate filter=${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue.txt --snp-loc CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.bim --gene-loc NCBI37.3.gene.loc --out ${sig[i]}_CLOZUK_PGC_SNPs_whole_genome
	magma --bfile CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS --gene-annot ${sig[i]}_CLOZUK_PGC_SNPs_whole_genome.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ${sig[i]}gene_annotation_for_CLOZUK_without_selecting_genes
done 
