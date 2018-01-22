#! /bin/bash

# Read in directory information
pathway=$1
path_to_scripts=$2
Gene_output_directory=$4
Random_scoring_directory=$5

# Assign the shell variables
source ${path_to_scripts}/PRS_arguments_script.sh
cat ${path_to_scripts}/PRS_arguments_script.sh

# Alter/add variables depending on what type of Training dataset you have
source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

# Overwrite arguments where needed
Name_of_extra_analysis=$3
permutations=$6

sig_thresholds=(`awk '{ print $3 }' ${training_set_name}_${validation_set_name}_plink_significance_thresholds_arguments_file_tmp.txt`)


# Create PRS profiles for each significance threshold specified
	
for i in "${sig_thresholds[@]}" ;
do

for l in `seq 1 ${permutations}` ;
do
	echo ${l}
	filename="${Random_scoring_directory}${pathway}_random_${l}_with_${i}"
	echo ${filename}

	plink --bfile ${Gene_output_directory}${validation_set_name}_${training_set_name}_normal_gene_regions_Clumped_whole_genome_final --score ${Random_scoring_directory}${pathway}_random_${l}_with_${i}.score --out ${filename}


done
done


