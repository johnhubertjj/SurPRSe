#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=24:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 1-22
#PBS -N c1020109_job_array_whole_genome

# script requries 22 files for each validation and training set

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" == *"raven"* ]]; then
  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
  
  cd $PBS_O_WORKDIR
  path_to_scripts="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/"
  path_to_PGC_conversion="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/MAGMA/"
  number_of_files=($(find -E . -type f -regex '^./output/CLOZUK_GWAS_BGE_CLUMPED_chr[0-9]+.bed' -exec basename {} \;))
  path_to_covariate_file="~/c1020109/NCBI37.3/CLOZUK2.r7.select2PC.eigenvec.txt" 
  path_to_chromosome_length_file="~/c1020109/NCBI37.3/UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="~/c1020109/NCBI37.3/CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="~/c1020109/NCBI37.3/NCBI37.3.gene.loc"  
  
  # assign arguments here for now because there are so many
  training_set_usually_summary="PGC_table"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr"
  training_set_name="PGC"
  validation_set_name="CLOZUK"
  MAF_summary="FALSE"
  MAF_threshold=0.01
  MAF_genotype="TRUE"
  INFO_summary="TRUE"
  INFO_threshold=0.9	
  Chromosomes_to_analyse=(`seq 1 22`)
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  sig_thresholds=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
  Perform_Magma_as_well="FALSE"
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates" 
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions= "both"

  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load magma/1.06

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  chromosome_number=22

  training_set_usually_summary="PGC_table${chromosome_number}"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr${chromosome_number}"
  training_set_name="PGC"
  validation_set_name="CLOZUK"
  MAF_summary="NO"
  MAF_threshold=0.01
  MAF_genotype="YES"
  INFO_summary="YES"
  INFO_threshold=0.6	
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

# Run R script that removes SNPs based on INFO score and MAF
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}MAF_and_INFO_score_summary_stats_script.R ./extrainfo/PGC_remove_MAF_INFO${chromosome_number}.Rout ${training_set_usually_summary} ${training_set_name} ${MAF_summary} ${INFO_summary} ${INFO_threshold} 

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}CLOZUK_PGC_COMBINE_final.R ./extrainfo/CLOZUK_PGC_COMBINE_chr${chromosome_number}.Rout ${training_set_usually_summary} ${validation_set_usually_genotype} ${training_set_name} ${validation_set_name}
if [[ ${MAF_genotype} == "YES" ]]; then
   # using plink to change the names to a CHR.POS identifier and remaking the files
  plink --bfile ${validation_set_usually_genotype} --update-name ./output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt --maf 0.01 --make-bed --out ./output/${validation_set_usually_genotype}_2
elif [[ ${MAF_genotype} == "NO" ]]; then 
  plink --bfile ${validation_set_usually_genotype} --update-name ./output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt --make-bed --out ./output/${validation_set_usually_genotype}_2
fi

# re-package the original files
# tar -czvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz CLOZUK_GWAS_BGE_chr${chromosome_number}.bed CLOZUK_GWAS_BGE_chr${chromosome_number}.bim CLOZUK_GWAS_BGE_chr${chromosome_number}.fam

# create files containing duplicate SNPs
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Clumping_CLOZUK_PGC.R ./extrainfo/CLOZUK_PGC_clumpinginfo_chr${chromosome_number}.Rout ${training_set_usually_summary} ${validation_set_usually_genotype} ${training_set_name} ${validation_set_name}  


if [ ${Peform_Magma_as_well} == "True" ]; then

# Create output directory for MAGMA results
	if [ ! d "./output/MAGMA_set_analysis" ]; then
		mkdir ./output/MAGMA_set_analysis
	fi
# Create files for MAGMA
plink --bfile ./output/${validation_set_usually_genotype}_2 --exclude ./output/extracted_Duplicate_snps_chr${chromosome_number}.txt --extract ./output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt --make-bed --out ./output/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates
fi

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
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --extract ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt --make-bed --out ./output/MAGMA_set_analysis/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}


# purge all modules
module purge




