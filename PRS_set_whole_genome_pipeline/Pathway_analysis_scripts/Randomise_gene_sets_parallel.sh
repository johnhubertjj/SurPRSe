#! /bin/bash

# Assign arguments
pathways=$1
path_to_scripts=$2
path_to_pathway_scripts=$3
path_to_Gene_scripts=$4
system=$5
Name_of_extra_analysis=$6
Gene_output_directory=$7

# Assign the shell variables
source ${path_to_scripts}/PRS_arguments_script.sh
cat ${path_to_scripts}/PRS_arguments_script.sh

# Alter/add variables depending on what type of Training dataset you have
source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}"

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_Gene_scripts}generate_random_andrews_script.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_generate_random_andrews_script.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${Gene_output_directory}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${Gene_regions}\
 ${permutations}\
 

