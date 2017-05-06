#!/bin/bash

#PBS -q serial
#PBS -P PR54
#PBS -l select=1:ncpus=1
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -o /home/c1020109/Summary_stats_info
#PBS -N c1020109_separate_summary_stats_whole_genome

cd $PBS_O_WORKDIR

# Script attempts to resolve what type of summary stats you have and how it influences the arguments further down the pipeline
echo "I'm running"
# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" == *"raven"* ]]; then
  # assign a new variable for the PBS_ARRAY_variable
  
  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11
  module load magma/1.06

  path_to_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  
  # make directories for output and extra info
  if [ ! -d "output" ]; then
     mkdir output
  fi
  
  if [ ! -d "extrainfo" ]; then
     mkdir extrainfo
  fi
  
  # Assign the shell variables
  source ${path_to_scripts}PRS_arguments_script.sh 
  cat ${path_to_scripts}PRS_arguments_script.sh 

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_cross_disorder/  
  
  # Arguments
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/'
   
  # make directories for output and extra info
  if [ ! -d "output" ]; then
     mkdir output
  fi
  
  if [ ! -d "extrainfo" ]; then
     mkdir extrainfo
  fi
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 

fi


Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Taking_in_summary_stat_data_for_chromosome_conversion.R ./extrainfo/${training_set_name}_summary_stats_for_chromosome_conversion.Rout ${training_set_original_filename} ${training_set_name} ${MAF_summary} ${INFO_summary} ${SE_summary} ${SE_threshold} ${Chromosomes_to_analyse[@]}  

source ./extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt 

