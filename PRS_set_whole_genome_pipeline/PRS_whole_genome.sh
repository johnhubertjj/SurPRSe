
echo "Whole Genome Full PRS commencing...."
 

if [[ "${Using_raven}" = "FALSE" ]]; then

# Quick fix for permissions on local and on the cluster
sudo parallel ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_WHOLE_GENOME_PRS_JOB_ARRAY.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system} ::: ${training_set_name} ::: ${validation_set_name}
else
# For use on local
parallel ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_WHOLE_GENOME_PRS_JOB_ARRAY.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system} ::: ${training_set_name} ::: ${validation_set_name} 
fi


source ${path_to_scripts}PRS_ANALYSIS_SERIAL_no_set.sh


