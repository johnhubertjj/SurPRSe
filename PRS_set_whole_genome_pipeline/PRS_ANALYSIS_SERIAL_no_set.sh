#!/bin/bash

#PBS -q serial
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o /home/c1020109/merging_datasets_PRS_information
#PBS -N merging_datasets_and_calculating_PRS

if [[ "$whereami" == *"raven"* ]]; then
	cd $PBS_O_WORKDIR
fi

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

system=$3

 #echo "Press CTRL+C to proceed."
 #trap "pkill -f 'sleep 1h'" INT
 #trap "set +x ; sleep 1h ; set -x" DEBUG


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
  path_to_scripts="/home/${USER}/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh

  # Alter/add variables depending on what type of training dataset you have 
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt  
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi

if [ "$system" == "MAC" || "$system" == "LINUX" ]; then

   
  # Arguments
  Directory_to_work_from=$1
  path_to_scripts=$2
  
  #Change to correct directory (just in case)
  cd ${Directory_to_work_from}

  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh

  # Alter/add variables depending on what type of training dataset you have 
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt  
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  chromosome_number=NA
fi 

# echo "Press CTRL+C to proceed."
# trap "pkill -f 'sleep 1h'" INT
# trap "set +x ; sleep 1h ; set -x" DEBUG

if [[ "$whereami" == *"raven"* ]]; then
	echo ${Batch_job_ID}
	# move the raven.OU files to one directory
	cd /home/$USER/
	shopt -s nullglob #enable
	set -- *${Batch_job_ID}*
	if [ "$#" -gt 0 ]; then
		if [[ ! -d "${Raven_out_info_directory}" ]]; then
			mkdir ${Raven_out_info_directory}
		fi
		cp *${Batch_job_ID}* ${Raven_out_info_directory}
	else
		echo "cannot find raven output scripts, check jobid from batch run (*POLYGENIC_RISK_SCORE* script) and add to end of *new_PRS_argument* script in relevant extrainfo directory"
		exit 1
	fi
	shopt -u nullglob #disable
	cd $PBS_O_WORKDIR
fi

if [ "$system" == "MAC" || "$system" == "LINUX" ]; then
	shopt -s nullglob
	set -- *Rtestoutclusterlogfile*
	if [ "$#" -gt 0 ]; then
		if [[ ! -d ${Raven_out_info_directory} ]]; then
			mkdir ${Raven_out_info_directory} 
		fi
		cp *Rtestoutclusterlogfile* ${Raven_out_info_directory}
	else
		  
		echo "cannot find raven output scripts, check jobid from batch run (*POLYGENIC_RISK_SCORE* script) and add to end of *new_PRS_argument* script in relevant extrainfo directory"
		exit 1
	fi
	shopt -u nullglob #disable
fi

# Run Rscript to find out the important information from the previous run
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}extracting_useful_SNP_information.R ./${training_set_name}_${validation_set_name}_extrainfo/extracting_useful_SNP_information.Rout ${training_set_name} ${validation_set_name} ${Raven_out_info_directory} ${INFO_summary} ${MAF_summary} ${MAF_threshold} ${INFO_threshold} ${SE_summary} ${SE_threshold} ${Chromosomes_to_analyse[@]}

# Check for make_list_file
if [ -f ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt ]; then
	rm ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
fi

# Create Make file from list of current files in directory
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}MAKE_list_file_creator.R ./${training_set_name}_${validation_set_name}_extrainfo/MAKE_list_file_creator.Rout ${validation_set_usually_genotype_serial} ${training_set_name} ${validation_set_name} ${Chromosomes_to_analyse[@]} 
  

# merge all the clumped chromosomes together
cd ./${training_set_name}_${validation_set_name}_output/
plink --merge-list ./${validation_set_name}_${training_set_name}_FULL_GENOME_CLUMPED_MAKE_LIST_INPUT.txt --make-bed --out ./${validation_set_name}_${training_set_name}_FULL_GENOME_CLUMPED
cd ..

# merge all the PGC SNPs together into one table
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}combining_summary_stats_tables_after_conversion_to_CHR_POS.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_conversion.Rout ${training_set_name} ${validation_set_name} ${Chromosomes_to_analyse[@]}  

# Create Final PRS outputs directory
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring" ]; then
  	mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring
fi

# Create score files for each significance threshold specified 
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}PRS_scoring_whole_genome.R ./${training_set_name}_${validation_set_name}_extrainfo/PRS_scoring_whole_genome.Rout ${training_set_name} ${validation_set_name} ${sig_thresholds[@]} 

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
	filename="./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_whole_genome_significance_threshold_at_${i}.score"
	plink --bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_FULL_GENOME_CLUMPED --score ${filename} --out ./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_whole_genome_significance_threshold_at_${i}
done

#Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}PRS_whole_genome_calc_log_regres_with_sig_thresh.R ./${training_set_name}_${validation_set_name}_extrainfo/PRS_genome_calc_log_regres.Rout ${training_set_name} ${validation_set_name} ${sig_thresholds[@]}

