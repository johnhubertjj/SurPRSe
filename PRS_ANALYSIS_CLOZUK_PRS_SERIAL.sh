#!/bin/bash

#PBS -q workq
#PBS -P PR54
#PBS -l select=1:ncpus=10:mem=100GB
#PBS -l walltime=10:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 21-22:2
#PBS -N c1020109_job_array_whole_genome

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" == *"raven"* ]]; then
  cd $PBS_O_WORKDIR
  path_to_scripts="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/"
  path_to_PGC_conversion="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/MAGMA/"
  number_of_files=($(find -E . -type f -regex '^./output/CLOZUK_GWAS_BGE_CLUMPED_chr[0-9]+.bed' -exec basename {} \;))	
  
  # assign arguments here for now because there are so many
  training_set_usually_summary="PGC_table"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr"
  training_set_name="PGC"
  validation_set_name="CLOZUK"
  MAF_summary="FALSE"
  MAF_threshold=0.01
  MAF_genotype="TRUE"
  INFO_summary="TRUE"
  INFO_threshold=0.6	
  Chromosomes_to_analyse=(`seq 1 22`)
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  sig_thresholds=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)

# Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.9-mpi
  module load magma...something

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/whole_genome_testing/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  number_of_files=($(find -E . -type f -regex '^./output/CLOZUK_GWAS_BGE_CLUMPED_chr[0-9]+.bed' -exec basename {} \;))  
  path_to_PGC_conversion="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation"
  path_to_CLOZUK_conversion="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation"
fi


# Find the length of the array containing the names of the files
# Note double-quotes to avoid extra parsing of funny characters in filenames
echo "${number_of_files[@]}" 
length_of_array=`echo "${#number_of_files[@]}"`

# Check for make_list_file
if [ -f ./output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt ]; then
  rm ./output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
fi

# Create Make-list file by repeating the name of each file on each line of the make-list text document
# 1 line per name to prefix each .bin/.fam/.bed format
num1=1
length_of_array=`echo "$((${length_of_array} - ${num1}))"`
 
for chromosome_number in `seq 0 ${length_of_array}` ;
do
current_file=$(basename ${number_of_files[$chromosome_number]} .bed)
printf "${current_file}\n%.0s" {1} >> ./output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
done

# echo "Press CTRL+C to proceed."
# trap "pkill -f 'sleep 1h'" INT
# trap "set +x ; sleep 1h ; set -x" DEBUG

# merge all the clumped chromosomes together
cd ./output/
plink --merge-list ./${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt --make-bed --out ./${validation_set_name}_${training_set_name}_FULL_GENOME
cd ..

# merge all the PGC SNPs together into one table
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_PGC_conversion}combining_summary_stats_tables_after_conversion_to_CHR_POS.R ./extrainfo/${training_set_name}_conversion.Rout ${training_set_usually_name} ${Multiple_Training_set_tables} ${Chromosomes_to_analyse}  

# recode genotypes for input into Python
plink --bfile ./output/${validation_set_name}_${training_set_name}_FULL_GENOME --recode A --out ./output/${validation_set_name}_${training_set_name}_FULL_GENOME
plink --bfile ./output/${validation_set_name}_${training_set_name}_FULL_GENOME --freq --out ./output/${validation_set_name}_${training_set_name}_FULL_GENOME

# Get the MAF from CLOZUK and import into the PGC_summary_stats for PRS analysis
awk '{ $6=$7=$8=$10=$12=$13=$14=$15=$16=$17=""; print$0 }' ./output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt > ./output/${training_set_name}_table_for_python.txt

# match the SNP rows to the MAF rows using R script
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Prepare_both_${validation_set_name}_AND_${training_set_name}_for_PRS_MAF_Genotype.R ./extrainfo/Prepare_both_${validation_set_name}_${training_set_name}_MAF_chr.Rout ${Running_in_Serial} ${training_set_name} ${validation_set_name}

# Final removal of headings for PRS analysis
awk 'NR>1' ./output/${validation_set_name}_${training_set_name}_FULL_GENOME.raw > ./output/${validation_set_name}_${training_set_name}_FULL_GENOME_no_head.raw

# Create Final PRS outputs directory
if [ ! -d "./output/PRS_scoring" ]; then
  	mkdir ./output/PRS_scoring
fi

#########################
#### MAGMA ANALYSIS #####
#########################
if [ ! -d "./extrainfo/MAGMA_analysis" ]; then
	mkdir ./extrainfo/MAGMA_analysis
fi
### WHERE I GOT TO ####
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Creation_of_filter_file_for_magma_gene.R ./extrainfo/Creation_of_filter_file_for_magma_gene.R ${training_set_name} ${validation_set_name} ${sig_thresholds[@]}


# Run with MAGMA
sh ${path_to_scripts}PRS_with_magma.sh ${sig_thresholds[@]} ${validation_set_name} ${training_set_name}

# Run preparation for annotation file for python scripts
R CMD BATCH ${path_to_scripts}MAGMA_python_annotation_table_creator.R ./extrainfo/MAGMA_annotation_table_creator.Rout

if [ "$whereami" == "raven13" ]; then
   	python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
else
   	python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
fi
 
