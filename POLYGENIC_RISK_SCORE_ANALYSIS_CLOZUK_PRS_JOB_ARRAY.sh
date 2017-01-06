#!/bin/bash

#PBS -q batch-long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime 10:00:00
#PBS -o PRS_CLOZUK_PGC_ANALYSIS.txt
#PBS -e PRS_CLOZUK_combine_PGC_ANALYSIS_error.txt
#PBS -j oe
#PBS -J 1-22
#PBS -r y

cd $PBS_O_WORKDIR

module load R/3.3.0
module load plink/1.9a

tar -zxvf CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX}.tar.gz
gzip -dc PGC_table${PBS_ARRAY_INDEX}.txt.gz

R CMD BATCH Clumping_CLOZUK_PGC.R 
 
plink 
 --bfile CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX}
 --extract chr${PBS_ARRAY_INDEX}PGC_CLOZUK_common_SNPs.txt
 --exclude extracted_Duplicate_snps_chr${PBS_ARRAY_INDEX}
 --clump PGC_table${PBS_ARRAY_INDEX}.txt
 --clump-kb 1000 
 --clump-p1 0.5 
 --clump-p2 0.5 
 --clump-r2 0.2 
 --out CLOZUK_PGC_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${PBS_ARRAY_INDEX}
 
tr -s ' ' '\t' < FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb.clumped_chr${PBS_ARRAY_INDEX} > CLOZUK_PGC_CLUMPED_chr${PBS_ARRAY_INDEX}.txt
cut -f 2,4,5,6 < CLOZUK_PGC_CLUMPED_chr${PBS_ARRAY_INDEX}.txt > CLOZUK_PGC_CLUMPED_FINAL_chr${PBS_ARRAY_INDEX}.txt
rm CLOZUK_PGC_CLUMPED_chr${PBS_ARRAY_INDEX}.txt

R CMD BATCH PRS_scoring_R_script.R




