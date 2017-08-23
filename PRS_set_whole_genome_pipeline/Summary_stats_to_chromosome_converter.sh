#!/bin/bash

#PBS -q serial
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=8GB
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

system=$3

if [[ "$whereami" == *"raven"* ]]; then
  # assign a new variable for the PBS_ARRAY_variable
  
  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11
  module load magma/1.06

  path_to_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  
  
  # Assign the shell variables
  source ${path_to_scripts}PRS_arguments_script.sh 
  cat ${path_to_scripts}PRS_arguments_script.sh 

  # make directories for output and extra info
  if [[ ! -d "${training_set_name}_${validation_set_name}_output" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_output
  fi
  
  if [[ ! -d "${training_set_name}_${validation_set_name}_extrainfo" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_extrainfo
  fi
fi

if [[ "$system" == "MAC" || "$system" == "LINUX"]]; then
  Directory_to_work_from=$1
  cd ${Directory_to_work_from}  
  
  # Arguments
  path_to_scripts=$2
   
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 

  if [[ ! -d "${training_set_name}_${validation_set_name}_output" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_output
  fi
  
  if [[ ! -d "${training_set_name}_${validation_set_name}_extrainfo" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_extrainfo
  fi
   
fi

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Taking_in_summary_stat_data_for_chromosome_conversion.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_summary_stats_for_chromosome_conversion.Rout ${training_set_original_filename} ${training_set_name} ${validation_set_name} ${MAF_summary} ${INFO_summary} ${SE_summary} ${SE_threshold} ${Chromosomes_to_analyse[@]}  

source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt 

if [ "${split_by_chromosome_required}" == "TRUE" ]; then
   for i in "${Chromosomes_to_analyse[@]}" ;
   do
	plink --bfile ${validation_set_full_name_without_chromosome} --chr ${i} --make-bed --out ${validation_set_usually_genotype_serial}${i}
	tar -zcvf ${validation_set_usually_genotype_serial}${i}.tar.gz ${validation_set_usually_genotype_serial}${i}.{bim,bed,fam,log}
	rm ${validation_set_usually_genotype_serial}${i}.{bed,bim,fam,log} 
   done
fi


if [[ "$whereami" == *"raven"* ]]; then
#Purge all modules
module purge
fi

