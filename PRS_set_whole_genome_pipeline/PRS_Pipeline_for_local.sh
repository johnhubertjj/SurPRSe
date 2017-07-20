#!/bin/bash
#exec &> logfile.txt
path_to_scripts="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
path_to_pathway_scripts="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts/"
Directory_to_work_from=`pwd`
chromosomes=(`seq 1 22`)

sh ${path_to_scripts}Summary_stats_to_chromosome_converter.sh ${Directory_to_work_from}


# Use sed to extract the argument Extra_analysis from your output scripts in question #output the arguments to a temporary file in order to read it
Extra_analyses=(`sed -n 's/Extra_analyses\=//p' ${path_to_scripts}PRS_arguments_script.sh`)
Name_of_extra_analysis=(`sed -n 's/Name_of_extra_analysis\=//p' ${path_to_scripts}PRS_arguments_script.sh`)


if [ "${Extra_analyses[1]}" = "TRUE" ]; then
	echo hi
	if [[ "${Name_of_extra_analysis[1]}" == "Pathways" ]]; then
		echo hi
		sh ${path_to_pathway_scripts}Pathway_analysis.sh ${Directory_to_work_from} 
	fi
else
	# calculate polygenic scores for the whole genome across different chromosomes	
#	sudo parallel ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh ::: ${chromosomes[@]} ::: ${Directory_to_work_from}
	
	# Need an alternative to Raven's log files to extract locally on the computer probably output important information to one file
	sh ${path_to_scripts}PRS_ANALYSIS_SERIAL_no_set.sh ${Directory_to_work_from}
fi
