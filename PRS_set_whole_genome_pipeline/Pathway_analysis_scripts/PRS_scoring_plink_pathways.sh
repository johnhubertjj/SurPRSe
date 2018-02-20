#! /bin/bash

pathways=$1
path_to_scripts=$2
path_to_pathway_scripts=$3
system=$4

# Assign the shell variables
source ${path_to_scripts}/PRS_arguments_script.sh
cat ${path_to_scripts}/PRS_arguments_script.sh

# Alter/add variables depending on what type of Training dataset you have
source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

# Overwrite arguments if argument is in array format of (Pathways Genes)
Name_of_extra_analysis=$5
Gene_regions=$6

Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}"

if [[ ${Gene_regions} = "both" || ${Gene_regions} = "normal" ]]; then

extra_gene_regions=normal

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
	filename="${Pathway_output_directory}/${pathways}/scoring_${training_set_name}_${validation_set_name}_pathway_${pathways}_${extra_gene_regions}_with_${i}"
	plink --bfile ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_${pathways}_${extra_gene_regions}_Clumped_whole_genome_final --score ${filename}.score --out ${filename}
done

fi

if [[ ${Gene_regions} = "both" || ${Gene_regions} = "extended" ]]; then

extra_gene_regions=extended

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
	filename="${Pathway_output_directory}/${pathways}/scoring_${training_set_name}_${validation_set_name}_pathway_${pathways}_${extra_gene_regions}_with_${i}"
	plink --bfile ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_${pathways}_${extra_gene_regions}_Clumped_whole_genome_final --score ${filename}.score --out ${filename}
done

fi

chmod -R a+rwx ${Pathway_output_directory}/${pathways} 

