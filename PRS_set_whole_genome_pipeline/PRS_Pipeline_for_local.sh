#!/bin/bash

path_to_scripts="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_pipeline_scripts/PRS_set_whole_genome_pipeline/"
path_to_pathway_scripts="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts/"
Directory_to_work_from="/Users/johnhubert/CLOZUK_ALSPAC_DATA_PRS"
chromosomes=(`seq 1 22`)
DATE=`date +%Y-%m-%d`

sh ${path_to_scripts}Summary_stats_to_chromosome_converter ~/Documents/CLOZUK_ALSPAC_DATA_PRS

parallel ./${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh ::: ${chromosomes[@]}

# Use sed to extract the argument Extra_analysis from your output scripts in question #output the arguments to a temporary file in order to read it
Extra_analysis=(`sed -n 's/Extra_analyses\=//p' ${path_to_scripts}PRS_arguments_script.sh`)
Name_of_extra_analysis=(`sed -n 's/Name_of_extra_analysis\=//p' ${path_to_scripts}PRS_arguments_script.sh`)

if [[ ${testing[2]} == "TRUE" ]]; then

	if [[ ${Name_of_extra_analysis[2]} == "Pathways" ]]; then
		sh ${path_to_pathway_scripts}Pathway_analysis.sh 
	fi
else

# Need an alternative to Raven's log files to extract locally on the computer probably output important information to one file
sh ${path_to_scripts}PRS_ANALYSIS_SERIAL_no_set.sh

