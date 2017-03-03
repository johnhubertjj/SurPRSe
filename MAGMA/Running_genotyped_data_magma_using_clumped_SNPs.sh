#!/bin/bash

cd ~/Dropbox/testing_PRS_chromosome_22/
chromosome_number=22

sig=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
for i in `seq 0 8` ;
do 
	# Create clumped plink files
	plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --extract ${sig[i]}update_Bed_files_for_MAGMA_analysis.txt --make-bed --out ${sig[i]}CLOZUK_GWAS_BGE_PRS_only_SNPS${chromosome_number}
	Rscript ~/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/MAGMA/Creation_of_filter_file_for_magma_gene.R > ${sig[i]}test.Rout ${sig[i]} ${sig[i]}CLOZUK_GWAS_BGE_PRS_only_SNPS${chromosome_number}.bim ${sig[i]}CLOZUK_PGC_P_Vals_reference_file_using_PRS.txt 	
	awk '{print $2}' ${sig[i]}CLOZUK_PGC_P_Vals_reference_file_using_PRS.txt > ${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue_using_PRS.txt
	
	magma --annotate filter=${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue_using_PRS.txt --snp-loc ${sig[i]}CLOZUK_GWAS_BGE_PRS_only_SNPS${chromosome_number}.bim --gene-loc NCBI37.3.gene.loc --out ${sig[i]}_CLOZUK_PGC_SNPs

	magma --bfile ${sig[i]}CLOZUK_GWAS_BGE_PRS_only_SNPS${chromosome_number} --gene-annot ${sig[i]}_CLOZUK_PGC_SNPs.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ${sig[i]}gene_annotation_for_CLOZUK_using_PRS_SNPS
 
done
