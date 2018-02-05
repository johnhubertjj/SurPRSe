path_to_scripts=$1
system=$2

Pathway_output_directory=$4
Pathway_output_directory=$(echo $Pathway_output_directory | sed 's_/_\\/_g')
Gene_output_directory=$5
Gene_output_directory=$(echo $Gene_output_directory | sed 's_/_\\/_g')

Random_directory=$6
Random_directory=$(echo $Random_directory | sed 's_/_\\/_g')

Random_scoring_directory=$7
Random_scoring_directory=$(echo $Random_scoring_directory | sed 's_/_\\/_g')

Pathway_names=${Pathway_output_directory}Pathway_names.txt
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
source ${path_to_scripts}PRS_arguments_script.sh

# training_set_name=CLOZUK1
sed -i -e "/Training_name:/ s/: .*/: ${training_set_name}/" "$3"
sed -i -e "/Validation_name:/ s/: .*/: ${validation_set_name}/" "$3"
sed -i -e "/Pathway_output_directory:/ s/: .*/: "${Pathway_output_directory}"/" "$3"
sed -i -e "/Gene_output_directory:/ s/: .*/: ${Gene_output_directory}/" "$3"
sed -i -e "/Gene_location_filename:/ s/: .*/: ${Gene_location_filename}/" "$3"
sed -i -e "/Gene_regions:/ s/: .*/: ${Gene_regions}/" "$3"
sed -i -e "/permutations:/ s/: .*/: ${permutations}/" "$3"
sed -i -e "/calculate_indep_SNPs:/ s/: .*/: ${calculate_indep_SNPs}/" "$3"
sed -i -e "/Random_directory:/ s/: .*/: ${Random_directory}/" "$3"
sed -i -e "/Random_scoring_directory:/ s/: .*/: ${Random_scoring_directory}/" "$3"
sed -i -e "/sample_replace:/ s/: .*/: ${sample_replace}/" "$3"
sed -i -e "/Pathway_names:/ s/: .*/: ${Pathway_names}/" "$3"
