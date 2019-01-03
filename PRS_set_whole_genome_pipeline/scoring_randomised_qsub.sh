#! /bin/bash

#PBS -q serial
#PBS -P PR54
#PBS -l ncpus=16
#PBS -l mem=45gb
#PBS -l walltime=24:00:00
#PBS -o /home/c1020109/Summary_stats_info


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
extra_path="/PhD_scripts"
system=LINUX #Too late to change now...its official, Raven runs on Linux because my scripts says so.
fi

path_to_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
path_to_pathway_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts/"
path_to_gene_scripts="${home_OS}${extra_path}/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Gene_analysis_scripts/"


source ${path_to_scripts}PRS_arguments_script.sh

if [ ${Using_raven} = "TRUE" ]; then
echo ${SLURM_SUBMIT_DIR}
cd $SLURM_SUBMIT_DIR
fi

Directory_to_work_from=`pwd`

log_file_name="${validation_set_name}_${training_set_name}_PRS_analysis"
# exec &> "${log_file_name}"_logfile.txt

cat ${path_to_scripts}PRS_arguments_script.sh

Name_of_extra_analysis=Pathways
Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/"
Pathway_file_name=${Pathway_output_directory}Pathway_names.txt
while IFS='' read -r line; do pathways+=("$line"); done <$Pathway_file_name
Gene_output_directory="./${training_set_name}_${validation_set_name}_output/Genes/"
Random_directory=${Pathway_output_directory}Randomised_gene_sets_analysis/
Random_scoring_directory=${Pathway_output_directory}Randomised_gene_sets_analysis/Scores/
pathways_for_randomisation=(`awk '{ print $1 }' ${Pathway_output_directory}${training_set_name}_${validation_set_name}_random_pathways_to_test.txt`)

# Run the thing you want to do over the servers
parallel ${path_to_pathway_scripts}test_script_randomised_plink.sh ::: ${pathways_for_randomisation[@]} ::: ${path_to_scripts} ::: ${Name_of_extra_analysis} ::: ${Gene_output_directory} ::: ${Random_scoring_directory} ::: ${permutations}

