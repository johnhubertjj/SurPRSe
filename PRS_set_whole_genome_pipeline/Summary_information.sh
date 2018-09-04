#! /bin/bash

Summary_analysis=FALSE

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Summary_table_for_analysis.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Summary_table_for_analysis.Rout\
 ${validation_set_full_name_without_chromosome}\
 ${training_set_name}\
 ${validation_set_name}\
 ${Summary_analysis}\
 ${Chromosomes_to_analyse[@]}


if [[ "${Using_raven}" = "FALSE" ]]; then

 # Quick fix for permissions on local and on the cluster
sudo parallel ${path_to_scripts}Summary_Collate_scripts.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system} ::: ${training_set_name} ::: ${validation_set_name}

else

# For use on local
parallel ${path_to_scripts}Summary_Collate_scripts.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system} ::: ${training_set_name} ::: ${validation_set_name}

fi

Summary_analysis=TRUE
 
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Summary_table_for_analysis.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Summary_table_for_analysis.Rout\
 ${validation_set_full_name_without_chromosome}\
 ${training_set_name}\
 ${validation_set_name}\
 ${Summary_analysis}\
 ${Chromosomes_to_analyse[@]}
