#!/bin/bash
 cd $PBS_O_WORKDIR
	
 
 # Script attempts to resolve what type of summary stats you have and how it influences the arguments further down the pipeline
 echo "I'm running"
 # Run locally or on ARCCA
 
 system=$3

if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then
  
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

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Taking_in_summary_stat_data_for_chromosome_conversion.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_summary_stats_for_chromosome_conversion.Rout\
 ${path_to_training_dataset}\
 ${training_set_name}\
 ${validation_set_name}\
 ${MAF_summary}\
 ${INFO_summary}\
 ${SE_summary}\
 ${SE_threshold}\
 ${Chromosomes_to_analyse[@]}  

source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt 

if [ "${split_by_chromosome_required}" = "TRUE" ]; then
   for i in "${Chromosomes_to_analyse[@]}" ;
   do
	plink --bfile ${path_to_validation_dataset} --chr ${i} --make-bed --out ${validation_set_usually_genotype_serial}${i}
	tar -zcvf ${validation_set_usually_genotype_serial}${i}.tar.gz ${validation_set_usually_genotype_serial}${i}.{bim,bed,fam,log}
	rm ${validation_set_usually_genotype_serial}${i}.{bed,bim,fam,log} 
   done
else
   
   for i in "${Chromosomes_to_analyse[@]}" ;
   do
	cp -v "${path_to_validation_dataset}_chr${i}.tar.gz" ${PBS_O_WORKDIR}
   done
fi 




if [[ "${Using_raven}" = "TRUE" ]]; then
#Purge all modules
module purge
fi

