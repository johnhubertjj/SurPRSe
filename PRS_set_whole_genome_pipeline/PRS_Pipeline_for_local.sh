#! /bin/bash

#PBS -q serial
#PBS -P PR54
#PBS -l ncpus=16
#PBS -l mem=40gb
#PBS -l walltime=24:00:00
#PBS -o /home/c1020109/Summary_stats_info_Biobank_PRS


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

elif [[ "$whereami" = *"raven"* ]]; then 
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


source ${path_to_scripts}PRS_arguments_script.sh

if [ ${Using_raven} = "TRUE" ]; then
echo ${PBS_O_WORKDIR}
cd $PBS_O_WORKDIR 
fi

Directory_to_work_from=`pwd`

log_file_name="${validation_set_name}_${training_set_name}_PRS_analysis"
 exec &> "${log_file_name}"_logfile.txt

cat ${path_to_scripts}PRS_arguments_script.sh


  if [[ ! -d "${training_set_name}_${validation_set_name}_output" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_output
  fi

  if [[ ! -d "${training_set_name}_${validation_set_name}_extrainfo" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_extrainfo
  fi

# use the require function to make sure that the packages are up to date with the software

# Create arguments files for easier input into Rscripts for parallelisation # 
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Chromosome_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_Chromosome_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Chromosomes_to_analyse[@]}

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}sig_thresholds_lower_bounds_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_sig_thresholds_lower_bounds_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${sig_thresholds_lower_bounds[@]}

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}sig_thresholds_plink_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_sig_thresholds_plink_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${sig_thresholds[@]}

${path_to_scripts}Summary_stats_to_chromosome_converter.sh ${Directory_to_work_from} ${path_to_scripts} ${system}

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
	${path_to_scripts}PRS_ANALYSIS_SERIAL_no_set.sh ${Directory_to_work_from} ${path_to_scripts} ${system}
fi

if [[ "${Using_raven}" = "TRUE" ]]; then
#Purge all modules
module purge
# rsync -avz . /neurocluster/filesync/c1020109/. 
fi

