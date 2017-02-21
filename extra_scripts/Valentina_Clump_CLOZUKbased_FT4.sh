#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N CLUMPFT4
#$ -t 1
#

PLINK=/share/apps/plink2

cd /home/wpcvm/SZ/Thyroid/GWAS/info0.8/CLUMPED

/share/apps/./plink2 \
  --bfile /home/wpcvm/SZ/Thyroid/GWAS/AntonioCheck/CLOZUK2016_BGE/bge_LOC_removed_A.T_C.G \
  --clump /home/wpcvm/SZ/Thyroid/GWAS/FT4_7.all_info0.8_hwe6_maf0.05.res \
  --clump-field P\
  --clump-kb 1000 \
  --clump-p1 0.5 \
  --clump-p2 0.5 \
  --clump-r2 0.2 \
  --clump-snp-field LOC \
  --out FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb

tr -s ' ' '\t' < FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb.clumped > tmp.txt
cut -f 2,4,5,6 < tmp.txt > tmp1.txt
head -n -2 tmp1.txt > FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb.clumped
cut -f 2 < FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb.clumped > FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb.clumped.snps

rm tmp.txt
rm tmp1.txt

exit;
