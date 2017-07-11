#! /bin/bash

cd ~/Documents/testing_PRS_chromosome_22/test_chr5/output
pathways=$1
path_to_scripts=$2

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

	if [ -f ./$Pathway_output_directory/${pathways}/make_full_plink_${pathways}.txt ]; then
		rm ./$Pathway_output_directory/${pathways}/make_full_plink_${pathways}.txt
	fi

	plink --bfile ${validation_set_usually_genotype_serial}${a}_magma_input --extract chromosome_${a}_${pathways}SNPs_for_clumping_regular_gene_regions.txt --make-bed --out chromosome_${a}_${validation_set_name}_${training_set_name}_${pathways}

#printf "chromosome_${number_of_files[i]}_CLOZUK_PGC_${pathways[a]} %.0s" {1..2} >> ./${pathways[a]}/make_full_plink_${pathways[a]}.txt
printf "chromosome_${a}_${validation_set_name}_${training_set_name}_${pathways}\n%.0s" {1} >> ./${Pathway_output_directory}/${pathways}/make_full_plink_${pathways}.txt  

done

plink --merge-list ./${Pathway_output_directory}/${pathways}/make_full_plink_${pathways}.txt --make-bed --out ./${Pathway_output_directory}/${pathways}/${pathways}_Clumping_input

plink --bfile ./${Pathway_output_directory}/${pathways}_Clumping_input --clump ${PATH_for_PGC}combined_${training_set_name}_table_with_CHR.POS_identifiers.txt --clump-p1 $p1 --clump-p2 $p2 --clump-r2 $r2 --clump-kb $window --out ${Pathway_output_directory}/${pathways}_${validation_set_name}_${training_set_name}_pathways_r0.2_removed_AT_CG_1000kb
 
# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < ${Pathway_output_directory}/${pathways}_${validation_set_name}_${training_set_name}_pathways_r0.2_removed_AT_CG_1000kb.clumped > ${Pathway_output_directory}/pathways_${validation_set_name}_${training_set_name}_CLUMPED_${pathways}.txt
cut -f 2,4,5,6 < ${Pathway_output_directory}/pathways_${validation_set_name}_${training_set_name}_CLUMPED_${pathways}.txt > pathways_${validation_set_name}_${training_set_name}_CLUMPED_FINAL_${pathways}.txt
rm ${Pathway_output_directory}/pathways_${validation_set_name}_${training_set_name}_CLUMPED_${pathways}.txt
awk '{ print $2 }' ${Pathway_output_directory}/pathways_${validation_set_name}_${training_set_name}_CLUMPED_FINAL_${pathways}.txt > ${Pathway_output_directory}/pathways_CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt
printf "%s\n\n" "$(tail -n +2 ${Pathway_output_directory}/pathways_CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt)" > ${Pathway_output_directory}/pathways_CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt
 
# Create clumped plink files
plink --bfile ${pathways}_Clumping_input --extract pathways_CLUMPED_EXTRACT_${validation_set_name}_${pathways}.txt --make-bed --out ${Pathway_output_directory}/pathways_${validation_set_name}_GWAS_BGE_CLUMPED_${pathways}

