#! /bin/bash

pathways=$1
path_to_scripts=$2
path_to_pathway_scripts=$3
system=$4
Name_of_extra_analysis=$5

# Assign the shell variables
source ${path_to_scripts}/PRS_arguments_script.sh
cat ${path_to_scripts}/PRS_arguments_script.sh

# Alter/add variables depending on what type of Training dataset you have
source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}"

for a in ${Chromosomes_to_analyse[@]};
do

	if [ ! -d $Pathway_output_directory/${pathways} ]; then
		mkdir $Pathway_output_directory/${pathways}
	fi

	if [ -f $Pathway_output_directory/${pathways}/make_full_plink_${pathways}.txt ]; then
		rm $Pathway_output_directory/${pathways}/make_full_plink_${pathways}.txt
	fi

	PRSACPJA_data_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}${a}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"
	
	plink --bfile ${PRSACPJA_data_output} --extract ${Pathway_output_directory}/chromosome_${a}_${pathways}_SNPs_for_clumping_normal_gene_regions.txt --make-bed --out ${Pathway_output_directory}/${pathways}/chromosome_${a}_${validation_set_name}_${training_set_name}_${pathways}

done

# Create Make file from list of current files in directory
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_pathway_scripts}MAKE_list_file_creator_for_sets.R ${Pathway_output_directory}/${pathways}/MAKE_list_file_creator.Rout ${validation_set_usually_genotype_serial} ${training_set_name} ${validation_set_name} ${Pathway_output_directory} ${pathways} ${Chromosomes_to_analyse[@]}

plink --merge-list ${Pathway_output_directory}/${pathways}/make_full_plink_${pathways}.txt --make-bed --out ${Pathway_output_directory}/${pathways}/${pathways}_Clumping_input

plink --bfile ${Pathway_output_directory}/${pathways}/${pathways}_Clumping_input --clump ${training_set_name}_${validation_set_name}_output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt --clump-p1 $p1 --clump-p2 $p2 --clump-r2 $r2 --clump-kb $window --out ${Pathway_output_directory}/${pathways}/${pathways}_${validation_set_name}_${training_set_name}_pathways_r0.2_removed_AT_CG_1000kb
 
# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < ${Pathway_output_directory}/${pathways}/${pathways}_${validation_set_name}_${training_set_name}_pathways_r0.2_removed_AT_CG_1000kb.clumped > ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_CLUMPED_${pathways}.txt
cut -f 2,4,5,6 < ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_CLUMPED_${pathways}.txt > ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_${pathways}.txt
rm ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_CLUMPED_${pathways}.txt
awk '{ print $2 }' ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_${pathways}.txt > ${Pathway_output_directory}/${pathways}/CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt
printf "%s\n\n" "$(tail -n +2 ${Pathway_output_directory}/${pathways}/CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt)" > ${Pathway_output_directory}/${pathways}/CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt
 
# Create clumped plink files
plink --bfile ${Pathway_output_directory}/${pathways}/${pathways}_Clumping_input --extract ${Pathway_output_directory}/${pathways}/CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt --make-bed --out ${Pathway_output_directory}/${pathways}/${validation_set_name}_${training_set_name}_${pathways}_Clumped_whole_genome_final

chmod -R a+rwx ${Pathway_output_directory}/${pathways} 

