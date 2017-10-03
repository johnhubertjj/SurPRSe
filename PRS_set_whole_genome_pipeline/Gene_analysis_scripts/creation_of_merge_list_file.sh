#! /bin/bash

Directory_to_work_from=$1
path_to_scripts=$2
path_to_gene_scripts=$3
system=$5

# Assign the shell variables
source ${path_to_scripts}/PRS_arguments_script.sh
cat ${path_to_scripts}/PRS_arguments_script.sh

# Alter/add variables depending on what type of Training dataset you have
source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

Gene_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}"
Name_of_extra_analysis=Genes

for a in ${Chromosomes_to_analyse[@]};
do

PRS_gene_bim_file_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"

if [[ ${Gene_regions} = "both" ]]; then
	Gene_regions=normal
	if [ -f $Gene_output_directory/make_full_plink_genes_normal.txt ]; then
		rm $Gene_output_directory/make_full_plink_genes_normal.txt
	fi
		
	plink --bfile ${PRS_gene_bim_file_output} --extract ${Gene_output_directory}/chromosome_${a}_SNPs_for_clumping_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/chromosome_${a}_${validation_set_name}_${training_set_name}_normal_gene_regions

# Create Make file from list of current files in directory
	Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}MAKE_list_file_creator_for_sets.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/MAKE_list_file_creator_gene_specific.Rout\
 ${validation_set_usually_genotype_serial}\
 ${training_set_name}\
 ${validation_set_name}\
 ${Gene_output_directory}\
 ${Gene_regions}\
 ${Chromosomes_to_analyse[@]}

	plink --merge-list ${Gene_output_directory}/make_full_plink_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions

	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --clump ${training_set_name}_${validation_set_name}_output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt --clump-p1 $p1 --clump-p2 $p2 --clump-r2 $r2 --clump-kb $window --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb
 
	# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
	tr -s ' ' '\t' < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb.clumped > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	cut -f 2,4,5,6 < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt
	rm ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	awk '{ print $2 }' ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
	printf "%s\n\n" "$(tail -n +2 ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt)" > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
 
	# Create clumped plink files
	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --extract ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_Clumped_whole_genome_final


	# add extended stuff here
	Gene_regions=extended
	if [ -f $Gene_output_directory/make_full_plink_genes_extended.txt ]; then
		rm $Gene_output_directory/make_full_plink_genes_extended.txt
	fi



	plink --bfile ${PRS_gene_bim_file_output} --extract ${Gene_output_directory}/chromosome_${a}_SNPs_for_clumping_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/chromosome_${a}_${validation_set_name}_${training_set_name}_normal_gene_regions

# Create Make file from list of current files in directory
	Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}MAKE_list_file_creator_for_sets.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/MAKE_list_file_creator_gene_specific.Rout\
 ${validation_set_usually_genotype_serial}\
 ${training_set_name}\
 ${validation_set_name}\
 ${Gene_output_directory}\
 ${Gene_regions}\
 ${Chromosomes_to_analyse[@]}

	plink --merge-list ${Gene_output_directory}/make_full_plink_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions

	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --clump ${training_set_name}_${validation_set_name}_output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt --clump-p1 $p1 --clump-p2 $p2 --clump-r2 $r2 --clump-kb $window --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb
 
	# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
	tr -s ' ' '\t' < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb.clumped > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	cut -f 2,4,5,6 < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt
	rm ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	awk '{ print $2 }' ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
	printf "%s\n\n" "$(tail -n +2 ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt)" > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
 
	# Create clumped plink files
	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --extract ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_Clumped_whole_genome_final

	Gene_regions=both
	chmod -R a+rwx ${Gene_output_directory}/

elif [[ ${Gene_regions} = "normal" ]]; then

 
	if [ -f $Gene_output_directory/make_full_plink_genes_normal.txt ]; then
		rm $Gene_output_directory/make_full_plink_genes_normal.txt
	fi
		
	plink --bfile ${PRS_gene_bim_file_output} --extract ${Gene_output_directory}/chromosome_${a}_SNPs_for_clumping_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/chromosome_${a}_${validation_set_name}_${training_set_name}_normal_gene_regions

# Create Make file from list of current files in directory
	Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}MAKE_list_file_creator_for_sets.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/MAKE_list_file_creator_gene_specific.Rout\
 ${validation_set_usually_genotype_serial}\
 ${training_set_name}\
 ${validation_set_name}\
 ${Gene_output_directory}\
 ${Gene_regions}\
 ${Chromosomes_to_analyse[@]}

	plink --merge-list ${Gene_output_directory}/make_full_plink_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions

	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --clump ${training_set_name}_${validation_set_name}_output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt --clump-p1 $p1 --clump-p2 $p2 --clump-r2 $r2 --clump-kb $window --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb
 
	# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
	tr -s ' ' '\t' < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb.clumped > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	cut -f 2,4,5,6 < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt
	rm ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	awk '{ print $2 }' ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
	printf "%s\n\n" "$(tail -n +2 ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt)" > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
 
	# Create clumped plink files
	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --extract ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_Clumped_whole_genome_final

elif [[ ${Gene_regions} = "extended" ]]; then

	if [ -f $Gene_output_directory/make_full_plink_genes_extended.txt ]; then
		rm $Gene_output_directory/make_full_plink_genes_extended.txt
	fi


	plink --bfile ${PRS_gene_bim_file_output} --extract ${Gene_output_directory}/chromosome_${a}_SNPs_for_clumping_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/chromosome_${a}_${validation_set_name}_${training_set_name}_normal_gene_regions

# Create Make file from list of current files in directory
	Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}MAKE_list_file_creator_for_sets.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/MAKE_list_file_creator_gene_specific.Rout\
 ${validation_set_usually_genotype_serial}\
 ${training_set_name}\
 ${validation_set_name}\
 ${Gene_output_directory}\
 ${Gene_regions}\
 ${Chromosomes_to_analyse[@]}

	plink --merge-list ${Gene_output_directory}/make_full_plink_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions

	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --clump ${training_set_name}_${validation_set_name}_output/combined_${training_set_name}_table_with_CHR.POS_identifiers.txt --clump-p1 $p1 --clump-p2 $p2 --clump-r2 $r2 --clump-kb $window --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb
 
	# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
	tr -s ' ' '\t' < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_r${r2}_p1_${p1}_p2_${p2}_removed_AT_CG_window_${window}kb.clumped > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	cut -f 2,4,5,6 < ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt > ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt
	rm ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_normal_gene_regions.txt
	awk '{ print $2 }' ${Gene_output_directory}/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_normal_gene_regions.txt > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
	printf "%s\n\n" "$(tail -n +2 ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt)" > ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt
 
	# Create clumped plink files
	plink --bfile ${Gene_output_directory}/${validation_set_name}_${training_set_name}_Clumping_input_normal_gene_regions --extract ${Gene_output_directory}/CLUMPED_EXTRACT_${validation_set_name}_normal_gene_regions.txt --make-bed --out ${Gene_output_directory}/${validation_set_name}_${training_set_name}_normal_gene_regions_Clumped_whole_genome_final

	chmod -R a+rwx ${Gene_output_directory}/
fi

