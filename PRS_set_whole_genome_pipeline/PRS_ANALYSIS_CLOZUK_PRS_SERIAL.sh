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

 echo "Press CTRL+C to proceed."
 trap "pkill -f 'sleep 1h'" INT
 trap "set +x ; sleep 1h ; set -x" DEBUG


if [[ "$whereami" == *"raven"* ]]; then
  
  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=NA

  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11
  module load magma/1.06

  cd $PBS_O_WORKDIR
  path_to_scripts="~/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh

  # Alter/add variables depending on what type of training dataset you have 
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments/for_${training_set_name}.txt  
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then

  cd /Volumes/PhD_storage/
  
  # Arguments
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/'
   
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh

  # Alter/add variables depending on what type of training dataset you have 
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments/for_${training_set_name}.txt  
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  chromosome_number=NA
fi 

 echo "Press CTRL+C to proceed."
 trap "pkill -f 'sleep 1h'" INT
 trap "set +x ; sleep 1h ; set -x" DEBUG

# Run Rscript to find out the important information from the previous run
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}extracting_useful_SNP_information.R ${training_set_name} ${validation_set_name} ${Raven_out_info_directory} ${INFO_summary} ${MAF_summary} ${MAF_threshold} ${INFO_threshold} ${SE_decision} ${SE_threshold} ${Chromosomes_to_analyse[@]}

   
# Find the length of the array containing the names of the files
# Note double-quotes to avoid extra parsing of funny characters in filenames
echo "${number_of_files[@]}" 
length_of_array=`echo "${#number_of_files[@]}"`
num1=1

# Check for make_list_file
if [ -f ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt ]; then
	rm ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
fi

# If running MAGMA as well
if [ ${Perform_Magma_as_well} == TRUE ]; then

	# Create output directory for MAGMA results
	if [ ! -d "./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis" ]; then
		mkdir ./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis
	fi
	
	# Calculate the number of files there are in the table and print them to the screen
	number_of_files_magma=($(find -E . -type f -regex "^./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/${validation_set_name}_GWAS_BGE_chr[0-9]+${Magma_validation_set_name}.bed" -exec basename {} \;))	
	length_of_array_magma=`echo "${#number_of_files_magma[@]}"`
	echo "${number_of_files_magma[@]}"
	length_of_array_magma=`echo "$((${length_of_array_magma} - ${num1}))"`

	# Create a make-list file specifically for magma results 
	if [ -f ./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt ]; then
		rm ./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
	fi
	
	for chromosome_number in `seq 0 ${length_of_array_magma}` ;
	do
		current_file=$(basename ${number_of_files_magma[$chromosome_number]} .bed)
		printf "${current_file}\n%.0s" {1} >> ./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
	done
	
	cd ./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/
	plink --merge-list ./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt --make-bed --out ./${validation_set_name}_${training_set_name}_MAGMA_FULL_GENOME
        cd ${PBS_O_WORKDIR}
	
	Validation_set_name_MAGMA="./${training_set_name}_${validation_set_name}_output/MAGMA_set_analysis/${validation_set_name}_${training_set_name}_MAGMA_FULL_GENOME"
	
fi
	
# Create Make-list file by repeating the name of each file on each line of the make-list text document
# 1 line per name to prefix each .bin/.fam/.bed format
length_of_array=`echo "$((${length_of_array} - ${num1}))"`
 
for chromosome_number in `seq 0 ${length_of_array}` ;
do
	current_file=$(basename ${number_of_files[$chromosome_number]} .bed)
	printf "${current_file}\n%.0s" {1} >> ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
done

# echo "Press CTRL+C to proceed."
# trap "pkill -f 'sleep 1h'" INT
# trap "set +x ; sleep 1h ; set -x" DEBUG

# merge all the clumped chromosomes together
cd ./${training_set_name}_${validation_set_name}_output/
plink --merge-list ./${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt --make-bed --out ./${validation_set_name}_${training_set_name}_FULL_GENOME
cd ..

exit 1

# merge all the PGC SNPs together into one table
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_PGC_conversion}combining_summary_stats_tables_after_conversion_to_CHR_POS.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_conversion.Rout ${training_set_usually_name} ${Multiple_Training_set_tables} ${Chromosomes_to_analyse[@]}  

# recode genotypes for input into Python
plink --bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME --recode A --out ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME
plink --bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME --freq --out ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME

# Get the MAF from CLOZUK and import into the PGC_summary_stats for PRS analysis
awk '{ $6=$7=$8=$10=$12=$13=$14=$15=$16=$17=""; print$0 }' ./${training_set_name}_${validation_set_name}_output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt > ./${training_set_name}_${validation_set_name}_output/${training_set_name}_table_for_python.txt

# match the SNP rows to the MAF rows using R script
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Prepare_both_${validation_set_name}_AND_${training_set_name}_for_PRS_MAF_Genotype.R ./${training_set_name}_${validation_set_name}_extrainfo/Prepare_both_${validation_set_name}_${training_set_name}_MAF_chr.Rout ${Running_in_Serial} ${training_set_name} ${validation_set_name}

# Final removal of headings for PRS analysis
awk 'NR>1' ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME.raw > ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_no_head.raw

# Create Final PRS outputs directory
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring" ]; then
  	mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring
fi

#########################
#### MAGMA ANALYSIS #####
#########################

if [ ${Perform_Magma_as_well} == "TRUE" ]; then
	# Create the pvalue reference file which matching SNPs to genes depending on the pvalue threshold
	Rscript ${path_to_MAGMA_scripts}RscriptEcho.R ${path_to_scripts}Creation_of_filter_file_for_magma_gene.R ./${training_set_name}_${validation_set_name}_extrainfo/Creation_of_filter_file_for_magma_gene.R ${training_set_name} ${validation_set_name} ${sig_thresholds[@]} ${validation_set_name_MAGMA}.bim ${Perform_Magma_as_well} ${INFO_summary} ${INFO_threshold}
        # Run MAGMA gene-set analysis for the whole genome
	sh ${path_to_scripts}MAGMA_analysis_whole_genome_complete.sh ${sig_thresholds[@]} ${validation_set_name} ${training_set_name} ${validation_set_name_MAGMA} ${Perform_Magma_as_well}

fi
	Perform_Magma_as_well="FALSE"
	Validataion_FULL_GENOME="./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME"
	Rscript ${path_to_MAGMA_scripts}RscriptEcho.R ${path_to_scripts}Creation_of_filter_file_for_magma_gene.R ./${training_set_name}_${validation_set_name}_extrainfo/Creation_of_filter_file_for_magma_gene.R ${training_set_name} ${validation_set_name} ${validation_FULL_GENOME}.bim ${Perform_Magma_as_well} ${sig_thresholds[@]} 

# Run with MAGMA
sh ${path_to_scripts}MAGMA_analysis_whole_genome_complete.sh ${sig_thresholds[@]} ${validation_set_name} ${training_set_name} ${validation_set_name_MAGMA} ${Perform_Magma_as_well} ${Gene_regions}

### THE FOLLOWING BELOW NEEDS TO BE IN A NEW SCRIPT ####
# Run preparation for annotation file for python scripts
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}MAGMA_python_annotation_table_creator_testing_whole_genome.R ./${training_set_name}_${validation_set_name}_extrainfo/MAGMA_annotation_table_creator_testing_whole_genome.Rout ${training_set_name} ${validation_set_name} ${path_to_chromosome_length_file} ${Validation_FULL_GENOME} ${path_to_gene_annotation_file} ${Gene_regions} ${Chromosomes_to_analyse[@]}   

# Make both directories for the profiles and the scores
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring/score" ]; then
	mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring/score
fi

if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring/score" ]; then
	mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring/Profiles
fi

# Calculate the Polygenic score using plink within an R script (for now)
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}PRS_scoring_whole_genome.R ./${training_set_name}_${validation_set_name}_extrainfo/PRS_scoring_whole_genome.Rout ${training_set_name} ${validation_set_name} ${Gene_regions} ${sig_thresholds[@]} 

if [ "$whereami" == "raven13" ]; then
   	python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
else
   	python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
fi
 
