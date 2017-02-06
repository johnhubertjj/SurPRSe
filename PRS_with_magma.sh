#!/bin/bash

chromosome_number=$1
echo $PWD
magma --annotate --snp-loc ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}.bim --gene-loc NCBI37.3.gene.loc --out ./output/CLOZUK_PRS_CLUMPED_chr${chromosome_number}

