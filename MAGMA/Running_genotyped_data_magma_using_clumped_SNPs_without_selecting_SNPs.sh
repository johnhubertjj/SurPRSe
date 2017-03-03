#!/bin/bash

cd ~/Dropbox/testing_PRS_chromosome_22/
chromosome_number=22

sig=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)

Rscript ~/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/MAGMA/Creation_of_filter_file_for_magma_gene.R > without_selecting_genes.Rout ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}.bim CLOZUK_PGC_P_Vals_reference_file_using_PRS_without_genes_selection.txt

for i in `seq 0 8` ;
do
        awk '{print $2}' ${sig[i]}CLOZUK_PGC_P_Vals_reference_file_without_selecting_genes.txt > ${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue_without_selecting_genes.txt

        magma --annotate filter=${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue_without_selecting_genes.txt --snp-loc ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}.bim --gene-loc NCBI37.3.gene.loc --out ${sig[i]}_CLOZUK_PGC_SNPs_without_selecting_snps

        magma --bfile ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number} --gene-annot ${sig[i]}_CLOZUK_PGC_SNPs_without_selecting_snps.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ${sig[i]}gene_annotation_for_CLOZUK_without_selecting_genes

done

