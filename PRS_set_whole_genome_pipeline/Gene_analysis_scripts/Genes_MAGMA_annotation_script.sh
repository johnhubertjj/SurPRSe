#!/bin/bash

chromosome_number=$1
Directory_to_work_from=$2
path_to_scripts=$3
path_to_gene_scripts=$4
system=$5

# Assign the shell variables
source ${path_to_scripts}/PRS_arguments_script.sh
cat ${path_to_scripts}/PRS_arguments_script.sh

# Alter/add variables depending on what type of Training dataset you have
source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt

Name_of_extra_analysis=Genes

gene_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}"
gene_bim_file=${gene_output_directory}/${validation_set_usually_genotype_serial}${chromosome_number}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates.bim
# input the right format here

if [[ ${Gene_regions} = "both" ]]; then
        magma --annotate window=35,10 --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${gene_output_directory}/${chromosome_number}_${training_set_name}_${validation_set_name}_SNPs_extended_gene_temp
	
        magma --annotate --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${gene_output_directory}/${chromosome_number}_${training_set_name}_${validation_set_name}_SNPs_gene_temp

elif [[ ${Gene_regions} = "extended" ]]; then
        magma --annotate window=35,10 --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${gene_output_directory}/${chromosome_number}_${training_set_name}_${validation_set_name}_SNPs_extended_gene_temp

elif [[ ${Gene_regions} = "normal" ]]; then
        magma --annotate --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${gene_output_directory}/${chromosome_number}_${training_set_name}_${validation_set_name}_SNPs_gene_temp

fi

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}Assign_SNPs_to_genes_whole_genome.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${chromosome_number}_assiging_SNPs_to_genes_whole_genome.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${validation_set_usually_genotype_serial}${chromosome_number}\
 ${gene_bim_file}\
 ${Name_of_extra_analysis}\
 ${path_to_stationary_data}${Gene_location_filename}\
 ${Gene_regions}\
 ${chromosome_number}

# Add annotation script (aka select out SNPs for each chromosome and then run plink to remake files only including genci SNPS (with or without regulatory regions) 
# END script        
