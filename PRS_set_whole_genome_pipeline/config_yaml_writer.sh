path_to_scripts=$1
system=$2

#if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then

#source ${path_to_scripts}/PRS_arguments_script.sh
#source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
#Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/"
#Pathway_file_name=${Pathway_output_directory}Pathway_names.txt
#while IFS='' read -r line; do pathways+=("$line"); done <$Pathway_file_name
#Gene_output_directory="./${training_set_name}_${validation_set_name}_output/Genes/"
#Random_directory=${Pathway_output_directory}Randomised_gene_sets_analysis/
#Random_scoring_directory=${Pathway_output_directory}Randomised_gene_sets_analysis/Scores/

#fi
training_set_name=CLOZUK1
sed -i -e "/Training_name:/ s/: .*/: ${training_set_name}/" "$3"


