#! /bin/bash

#PBS -q serial_long
#PBS -P PR54
#PBS -l select=1:ncpus=12:mem=40GB
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -o /home/c102019/Summary_stats_info
#PBS -N PRS_pipeline_test

log_file_name=${1}
exec &> "${log_file_name}"_logfile.txt

whereami=$(uname -n)
echo "$whereami"

if [ "$whereami" = "v1711-0ab8c3db.mobile.cf.ac.uk" ]; then
home_OS="/Users"
extra_path="/johnhubert/Documents/PhD_scripts"
system=MAC

elif [ "$whereami" = "johnhubert-ThinkPad-P50" ]; then
home_OS="/home"
extra_path="/johnhubert/Documents"
system=LINUX # oh god programmers are going to hate me for using this argument

elif [[ "$whereami" = *"raven"* ]]; then 
home_OS=${HOME}
extra_path=/PhD_scripts"
system=LINUX #Too late to change now...its official, Raven runs on Linux because my scripts says so.
fi

path_to_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
path_to_pathway_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts/"
path_to_gene_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Gene_analysis_scripts/"

Directory_to_work_from=`pwd`

source ${path_to_scripts}PRS_arguments_script.sh
cat ${path_to_scripts}PRS_arguments_script.sh

# Create arguments files for easier input into Rscripts for parallelisation # 
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Chromosome_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_Chromosome_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Chromosomes_to_analyse[@]}\

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}sig_thresholds_lower_bounds_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_sig_thresholds_lower_bounds_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${sig_thresholds_lower_bounds[@]}\

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}sig_thresholds_plink_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_sig_thresholds_plink_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${sig_thresholds[@]}\

${path_to_scripts}Summary_stats_to_chromosome_converter.sh ${Directory_to_work_from} ${path_to_scripts} ${system}

# calculate polygenic scores for the whole genome across different chromosomes	

sudo parallel ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system}

# exit 0

# Use sed to extract the argument Extra_analysis from your output scripts in question #output the arguments to a temporary file in order to read it
Extra_analyses=(`sed -n 's/Extra_analyses\=//p' ${path_to_scripts}PRS_arguments_script.sh`)
Name_of_extra_analysis=(`sed -n 's/Name_of_extra_analysis\=//p' ${path_to_scripts}PRS_arguments_script.sh`)


if [ "${Extra_analyses[1]}" = "TRUE" ]; then
	echo hi
length_of_extra_analysis_array=`echo ${#Name_of_extra_analysis[@]}`
	if [ "${length_of_extra_analysis_array}" -gt  "2" ]; then
		echo hi all
		Name_of_extra_analysis_specific=(Genes Pathways)
		#${path_to_pathway_scripts}Pathway_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_pathway_scripts} ${Name_of_extra_analysis_specific[1]} 
		 ${path_to_gene_scripts}Gene_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_gene_scripts} ${Name_of_extra_analysis_specific[0]}
	
	elif  [[ "${Name_of_extra_analysis[1]}" = "Pathways" ]]; then
		echo hi Pat
		echo ${Name_of_extra_analysis[1]}
 ${path_to_pathway_scripts}Pathway_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_pathway_scripts} ${Name_of_extra_analysis[1]} 
	
	elif [[ "${Name_of_extra_analysis[1]}" = "Genes" ]]; then
		echo hi Gene
		 ${path_to_gene_scripts}Gene_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_gene_scripts} ${Name_of_extra_analysis[1]}
	fi
else
	# Need an alternative to Raven's log files to extract locally on the computer probably output important information to one file
	${path_to_scripts}PRS_ANALYSIS_SERIAL_no_set.sh ${Directory_to_work_from} ${path_to_scripts} ${system}
fi

if [[ "${Using_raven}" = "TRUE" ]]; then
#Purge all modules
module purge
fi

