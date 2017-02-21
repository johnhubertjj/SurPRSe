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

echo "Press CTRL+C to proceed."
trap "pkill -f 'sleep 1h'" INT
trap "set +x ; sleep 1h ; set -x" DEBUG

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" == *"raven"* ]]; then
  cd $PBS_O_WORKDIR
  path_to_scripts="~/c1020109/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/"
  path_to_PGC_conversion="~/c1020109/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation"
  path_to_CLOZUK_conversion="~/c1020109/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation"
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

input_file_name="./output/CLOZUK_GWAS_BGE_CLUMPED_chr"

# List the number of files matching the output from the previous Job array script "POLYGENIC_RISK_SCORE*.sh"
shopt -s nullglob
number_of_files=(*${input_file_name}*)
shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later

# Find the length of the array containing the names of the files
# Note double-quotes to avoid extra parsing of funny characters in filenames
echo "${number_of_files[@]}" 
length_of_array=`echo "${#number_of_files[@]}"`

# Check for make_list_file
if [ -f ./output/CLOZUK_PGC_FULL_GENOME_MAKE_LIST.txt ]; then
  rm ./${pathways[a]}/make_full_plink_${pathways[a]}.txt
fi

# Create Make-list file by repeating the name of each file on each line of the make-list text document
# 1 line per name to prefix each .bin/.fam/.bed format 
for chromosome_number in `seq 1 ${length_of_array}` ;
printf "${input_file_name}_${number_of_files[chromosome_number]}\n%.0s" {1} >> ./output/CLOZUK_PGC_FULL_GENOME_MAKE_LIST_INPUT.txt
done

# merge all the clumped chromosomes together
plink --merge-list ./output/CLOZUK_PGC_FULL_GENOME_MAKE_LIST.txt --make-bed --out ./output/CLOZUK_PGC_FULL_GENOME
 
# merge all the PGC SNPs together into one table
R CMD BATCH ${path_to_PGC_conversion}combining_summary_stats_tables_after_conversion_to_CHR_POS.R ./extrainfo/PGC_conversion.Rout 

# merge all the CLOZUK SNPs together into one table
R CMD BATCH ${path_to_CLOZUK_conversion}combining_clumped_CLOZUK_BGE_datasets.R ./extrainfo/PGC_conversion.Rout 

# recode genotypes for input into Python
plink --bfile ./output/CLOZUK_PGC_FULL_GENOME --recode A --out ./output/CLOZUK_PGC_FULL_GENOME
plink --bfile ./output/CLOZUK_PGC_FULL_GENOME --freq --out ./output/CLOZUK_PGC_FULL_GENOME

# Get the MAF from CLOZUK and import into the PGC_summary_stats for PRS analysis
awk '{ $6=$7=$8=$10=$12=$13=$14=$15=$16=$17=""; print$0 }' ./output/combined_PGC_table_with_CHR.POS_identifiers.txt > ./output/PGC_table_for_python.txt

# match the SNP rows to the MAF rows using R script
R CMD BATCH 'Serial' ${path_to_scripts}Prepare_both_CLOZUK_AND_PGC_for_PRS_MAF_Genotype.R ./extrainfo/Prepare_both_CLOZUK_PGC_MAF_chr.Rout

# Final removal of headings for PRS analysis
awk 'NR>1' ./output/CLOZUK_PGC_FULL_GENOME.raw > ./output/CLOZUK_PGC_FULL_GENOME_no_head.raw

# Create Final PRS outputs
if [ ! -d "./output/PRS_scoring" ]; then
  	mkdir ./output/PRS_scoring
fi

# Needs the MAGMA script
sh ${path_to_scripts}PRS_with_magma.sh ${chromosome_number}

# Run preparation for annotation file for python scripts
R CMD BATCH ${path_to_scripts}MAGMA_python_annotation_table_creator.R ./extrainfo/MAGMA_annotation_table_creator.Rout

# Run PRS python script # Loop through significance thresholds
# Change for Raven or Local
if [ "$whereami" == "raven13" ]; then
   	python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
else
   	python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
fi
 
