#!/bin/bash

chromosome_number=$1
echo $PWD
whereami=$(uname -n)
echo "$whereami"
if [ "$whereami" == "raven13" ]; then
	/home/c1020109/magma --annotate --snp-loc ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}.bim --gene-loc NCBI37.3.gene.loc --out ./output/CLOZUK_PRS_CLUMPED_chr${chromosome_number}
else 
	magma --annotate --snp-loc ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}.bim --gene-loc NCBI37.3.gene.loc --out ./output/CLOZUK_PRS_CLUMPED_chr${chromosome_number}
