#!/bin/bash

#PBS -q workq
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=24:00:00
#PBS -o /home/c1020109/parse_summary_stats.log
#PBS -j oe
#PBS -N c1020109_Parse_summary_stats_whole_genome

# Script attempts to resolve what type of summary stats you have and how it influences the arguments further down the pipeline

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" == *"raven"* ]]; then
  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
  
  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11
  module load magma/1.06

  cd $PBS_O_WORKDIR
  path_to_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  
  # Assign the shell variables
  source ${path_to_scripts}PRS_arguments_script.sh 
  cat ${path_to_scripts}PRS_arguments_script.sh 

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_cross_disorder/
  
  # Arguments
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/'
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 
  
fi

#  echo "Press CTRL+C to proceed."
#  trap "pkill -f 'sleep 1h'" INT
#  trap "set +x ; sleep 1h ; set -x" DEBUG

# Find out what type the input file is...quickly
shopt -s nullglob #enable
set -- *${training_set_usually_summary}*
if [ "$#" -eq 1 ]; then
	file_extension=$(for f in *${training_set_usually_summary}*; do printf "%s\n" "${f##*.}"; done | sort -u) 
	if [[ ${file_extension} == "txt" ]]; then
		training_set_full_path=`echo *${training_set_usually_summary}*`
		echo "input file is .txt file"
	elif [[ ${file_extension} != "txt" ]]; then
		summary_input_files=$(file *${training_set_usually_summary}* | cut -f2 -d' ')
		if [[ ${summary_input_files} == "ASCII" ]]; then
			training_set_full_path=`echo *${training_set_usually_summary}*`
			echo "input file is ASCII file"
		elif [[ ${summary_input_files} != "ACSII" ]]; then
			summary_input_files2=$(file *${training_set_usually_summary}* | cut -f3 -d' ')
			if [[ ${summary_input_files2} == "compressed" || ${summary_input_files2} == "archive" ]]; then
			        training_set_full_path=`echo *${training_set_usually_summary}*`
				${path_to_scripts}unarc ${training_set_full_path}
			        echo "hi"
			fi
		fi
	fi
elif [ "$#" -gt 1 ]; then
echo "${training_set_usually_summary} is either mispelt or has more than one input file. Please use only one input file"
exit 1
fi
shopt -u nullglob # disable

  
