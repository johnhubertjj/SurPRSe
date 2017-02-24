#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=48:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 4-5:2
#PBS -N c1020109_job_array_whole_genome


# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" == *"raven"* ]]; then
  WDPATH=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC
  cd $WDPATH
  path_to_scripts='/home/c1020109/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'

  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.9-mpi

  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  chromosome_number=22
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

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
R CMD BATCH ${path_to_scripts}CLOZUK_PGC_COMBINE_final.R ./extrainfo/CLOZUK_PGC_COMBINE_chr${chromosome_number}.Rout
 
# using plink to change the names to a CHR.POS identifier and remaking the files
plink --bfile CLOZUK_GWAS_BGE_chr${chromosome_number} --update-name ./output/CLOZUK_chr${chromosome_number}_chr.pos.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2

# re-package the original files
# tar -czvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz CLOZUK_GWAS_BGE_chr${chromosome_number}.bed CLOZUK_GWAS_BGE_chr${chromosome_number}.bim CLOZUK_GWAS_BGE_chr${chromosome_number}.fam

# create files containing duplicate SNPs
R CMD BATCH ${path_to_scripts}Clumping_CLOZUK_PGC.R ./extrainfo/CLOZUK_PGC_clumpinginfo_chr${chromosome_number}.Rout 

# Create files for MAGMA
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --exclude ./output/extracted_Duplicate_snps_chr${chromosome_number}.txt --extract ./output/chr${chromosome_number}PGC_CLOZUK_common_SNPs.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_chr22_magma_input


# Clump the datasets
# Extract the SNPs common between PGC and CLOZUK
# Remove the duplicate SNPs
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --extract ./output/chr${chromosome_number}PGC_CLOZUK_common_SNPs.txt --exclude ./output/extracted_Duplicate_snps_chr${chromosome_number}.txt --clump ./output/PGC_table${chromosome_number}_new.txt --clump-kb 1000 --clump-p1 0.5 --clump-p2 0.5 --clump-r2 0.2 --out ./output/CLOZUK_PGC_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${chromosome_number}
 
# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < ./output/CLOZUK_PGC_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${chromosome_number}.clumped > ./output/CLOZUK_PGC_CLUMPED_chr${chromosome_number}.txt
cut -f 2,4,5,6 < ./output/CLOZUK_PGC_CLUMPED_chr${chromosome_number}.txt > ./output/CLOZUK_PGC_CLUMPED_FINAL_chr${chromosome_number}.txt
rm ./output/CLOZUK_PGC_CLUMPED_chr${chromosome_number}.txt
awk '{ print $2 }' ./output/CLOZUK_PGC_CLUMPED_FINAL_chr${chromosome_number}.txt > ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt
printf "%s\n\n" "$(tail -n +2 ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt)" > ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt 

# Create clumped plink files
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --extract ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}





