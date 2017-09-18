#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=10:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 21-22:2
#PBS -N c1020109_job_array

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" = *"raven"* ]]; then
 
  training_set_usually_summary="PGC_table"${chromosome_number}
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr"${chromosome_number}
  training_set_name="PGC"
  validation_set_name="CLOZUK"
  MAF_summary="NO"
  MAF_genotype="YES"
  INFO_summary="YES"
  INFO_threshold=0.6
  split_genotype_dataset_into_chromosomes="YES"
  WDPATH=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC
  
  cd $WDPATH
  path_to_scripts='/home/c1020109/PRS_scripts/'
  # Load both Plink and R
  module purge
  module load R/3.1.0
  module load plink/1.9c3
  module load python/2.7.9-mpi

  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
elif [ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  chromosome_number=22

  training_set_usually_summary="PGC_table"${chromosome_number}
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr"${chromosome_number}
  training_set_name="PGC"
  validation_set_name="CLOZUK"
  MAF_summary="NO"
  MAF_threshold=0.01
  MAF_genotype="YES"
  INFO_summary="YES"
  INFO_threshold=0.6	
  split_genotype_dataset_into_chromosomes="NO"
fi

## rewrite so that the file input is an argument for the script instead, this will work for now
shopt -s nullglob
set -- *${chromosome_number}.tar.gz
if [ "$#" -gt 0 ]; then

# unpack the CLOZUK datasets
tar -zxvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz
fi

# make directories for output and extra info
if [ ! -d "output" ]; then
   mkdir output
fi

if [ ! -d "extrainfo" ]; then
   mkdir extrainfo
fi

if [[ ${split_genotype_dataset_into_chromosomes} = "YES" ]]; then
   Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}splitting_Genotype_datasets.R ./extrainfo/splitting_${Validation_set_name}_genotype_dataset.Rout

# Run R script that removes SNPs based on INFO score and MAF
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}MAF_and_INFO_score_summary_stats_script.R ./extrainfo/PGC_remove_MAF_INFO.Rout ${training_set_usually_summary} ${training_set_name} ${MAF_summary} ${INFO_summary} ${INFO_threshold} 

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}CLOZUK_PGC_COMBINE_final.R ./extrainfo/CLOZUK_PGC_COMBINE_chr${chromosome_number}.Rout ${training_set_usually_summary} ${validation_set_usually_genotype}

if [[ ${MAF_genotype} = "YES" ]]; then
# using plink to change the names to a CHR.POS identifier and remaking the files
  plink --bfile CLOZUK_GWAS_BGE_chr${chromosome_number} --update-name ./output/CLOZUK_chr${chromosome_number}_chr.pos.txt --maf 0.01 --make-bed --out ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2
elif [[ ${MAF_genotype} = "NO" ]]; then 
  plink --bfile CLOZUK_GWAS_BGE_chr${chromosome_number} --update-name ./output/CLOZUK_chr${chromosome_number}_chr.pos.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2
fi

#re-package the original files
#tar -czvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz CLOZUK_GWAS_BGE_chr${chromosome_number}.bed CLOZUK_GWAS_BGE_chr${chromosome_number}.bim CLOZUK_GWAS_BGE_chr${chromosome_number}.fam

# create files containing duplicate SNPs
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Clumping_CLOZUK_PGC.R ./extrainfo/CLOZUK_PGC_clumpinginfo_chr${chromosome_number}.Rout ${training_set_usually_summary} ${validation_set_usually_genotype} ${training_set_name} ${validation_set_name}  

# Create files for MAGMA
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --debug --memory 20000 --exclude ./output/extracted_Duplicate_snps_chr${chromosome_number}.txt --extract ./output/chr${chromosome_number}PGC_CLOZUK_common_SNPs.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_magma_input
