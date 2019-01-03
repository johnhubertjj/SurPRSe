#! /bin/bash

#SBATCH -p htc
#SBATCH --account=scw1143
#SBATCH --tasks-per-node=22
#SBATCH --mem-per-cpu=10G
#SBATCH -t 1-00:00:00
#SBATCH --job-name=PRS_tutorial_test
#SBATCH -o /home/c.c1020109/Summary_stats_info_Biobank_PRS

echo "hi"

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

else

home_OS=${HOME}

# extra_path must be NULL or a path
extra_path=NULL
system=LINUX #Too late to change now...its official, Raven runs on Linux because my scripts says so.
fi

if [[ "${extra_path}" != "NULL" ]]; then
path_to_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
path_to_pathway_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts/"
path_to_gene_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Gene_analysis_scripts/"

else
path_to_scripts="${home_OS}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
path_to_pathway_scripts="${home_OS}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts/"
path_to_gene_scripts="${home_OS}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Gene_analysis_scripts/"
fi


# Preparation script; only alter if stated within the wiki/tutorial of this software
source ${path_to_scripts}Pipeline_Preparation.sh

log_file_name="${validation_set_name}_${training_set_name}_PRS_analysis"
 exec &> "${log_file_name}"_logfile.txt

# Create argument text files as alternative arguments to other scripts when defining chromosomes or p-value significance thresholds
source ${path_to_scripts}Plink_arguments_files_R_scripts.sh

# Convert input files into a chromosome format
source ${path_to_scripts}Summary_stats_to_chromosome_converter.sh

# calculate polygenic scores for the whole genome across different chromosomes

if [[ "${Using_raven}" = "FALSE" ]]; then

# Quick fix for permissions on local and on the cluster
sudo parallel ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system} ::: ${training_set_name} ::: ${validation_set_name}
else
# For use on local
parallel ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${system} ::: ${training_set_name} ::: ${validation_set_name}
fi

# exit 0

# Use sed to extract the argument Extra_analysis from your output scripts in question #output the arguments to a temporary file in order to read it

# Extra_analyses=(`sed -n 's/Extra_analyses\=//p' ${path_to_scripts}PRS_arguments_script.sh`)
# Name_of_extra_analysis=(`sed -n 's/Name_of_extra_analysis\=//p' ${path_to_scripts}PRS_arguments_script.sh`)


if [ "${Extra_analyses}" = "TRUE" ]; then
	echo hi
length_of_extra_analysis_array=`echo ${#Name_of_extra_analysis[@]}`
	if [ "${length_of_extra_analysis_array}" -eq  "2" ]; then
		echo hi all
		Name_of_extra_analysis_specific=(Pathways Genes)
		${path_to_gene_scripts}Gene_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_gene_scripts} ${Name_of_extra_analysis_specific[1]}
		${path_to_pathway_scripts}Pathway_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_pathway_scripts} ${path_to_gene_scripts} ${Name_of_extra_analysis_specific[0]}

	elif  [[ "${Name_of_extra_analysis[0]}" = "Pathways" ]]; then
		echo hi Pat
		echo ${Name_of_extra_analysis[0]}
 		${path_to_pathway_scripts}Pathway_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_pathway_scripts} ${Name_of_extra_analysis[0]}

	elif [[ "${Name_of_extra_analysis[1]}" = "Genes" ]]; then
		echo hi Gene
		${path_to_gene_scripts}Gene_analysis.sh ${Directory_to_work_from} ${path_to_scripts} ${system} ${path_to_gene_scripts} ${Name_of_extra_analysis[1]}
	fi


else
	# Need an alternative to Raven's log files to extract locally on the computer probably output important information to one file

# Source the altered arguments from earlier, unsure if it would have already been performed

 source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
 cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

# Source the final stages of the whole genome analysis...output files could use some structure though
 source	${path_to_scripts}PRS_ANALYSIS_SERIAL_no_set.sh

fi

if [[ "${Full_genome_PRS_extra_analysis}" = "TRUE" ]]; then

source ${path_to_scripts}PRS_whole_genome.sh

fi

if [[ "${MAGMA_gene_set_analysis}" = "TRUE" ]]; then

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}MAGMA_extract_SNP_list.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_MAGMA_extract_SNP_list.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${validation_set_usually_genotype_serial}\
 ${Chromosomes_to_analyse[@]}

source ${path_to_scripts}MAGMA_gene_set_analysis.sh

fi

if [ "${Extra_analyses}" = "TRUE" ]; then
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Collate_all_PRS_files_together.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_combine_PRS_results.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Extra_analyses}\
 ${Full_genome_PRS_extra_analysis}\
 ${Gene_regions}\
 ${whole_genome_genic}\
 ${Gene_specific_PRS}\
 ${sig_thresholds[@]}
fi

 
if [[ "${Using_raven}" = "TRUE" ]]; then
#Purge all modules
module purge
# rsync -avz . /neurocluster/filesync/c1020109/.
fi
