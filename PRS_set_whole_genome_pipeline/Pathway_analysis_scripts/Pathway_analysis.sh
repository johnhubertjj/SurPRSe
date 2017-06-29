#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=45GB
#PBS -l walltime=24:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 1-22
#PBS -N c1020109_job_array_whole_genome

# script requries 22 files for each validation and training set

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" == *"raven"* ]]; then
  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
  
  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11
  module load magma/1.06

  cd $PBS_O_WORKDIR
  path_to_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  path_to_pathway_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts"
  # Assign the shell variables
  source ${path_to_scripts}PRS_arguments_script.sh 
  cat ${path_to_scripts}PRS_arguments_script.sh 
   
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
 
elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_cross_disorder/
  
  # Arguments
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/'
  path_to_pathway_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts"
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 
  
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi 

if [[ ${Name_of_extra_analysis} == "Pathways" ]]; then
 
# Create output directory for Extra_analysis results
	if [ ! -d "./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}" ]; then
        	mkdir ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
	fi
 
Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}PATHWAYS_PRS_COLLECTING_MAGMA_INFO.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/PATHWAYS_PRS_COLLECTING_MAGMA_INFO.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${valudation_set_usually_genotype_serial}\
 ${Name_of_extra_analysis}\
 ${Pathway_file_name}\
 ${Gene_location_filename}\
 ${Chromosomes_to_analyse[@]}
 
 

# From the above script, identify the number of pathways you want to analyse (probably safest to write to a file, port to a variable and then delete the file)
# Also a text-delimited file with each line specifying a pathway name to be used
# The seperate gene_loc files belonging to previously specified analysis
# Then use magma to create annotation files:

PRSACPJA_data_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}${Chromosomes_to_analyse[@]}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"
Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/" 
# separate script below:

if [[ ${Gene_regions} == "Extended" ]]; then

	while IFS='' read -r line || [[ -n "$line" ]]; 
	do
		for i in ${Chromosomes_to_analyse[@]};
		do		    
			magma --annotate window=35,10 --snp-loc ${PRSACPJA_data_output} --gene-loc ${Pathway_output_directory}${training_set_name}_${validation_set_name}_${line}_chromosome_${Chromosomes_to_analyse[i]}_extended_temp.gene.loc --out ${Pathway_output_directory}${Chromosomes_to_analyse[i]}_${training_set_name}_${validation_set_name}_SNPs_${line}_extended_pathway_temp
		done

	done < "$1"
fi

if [[ ${Gene_regions} == "Regular" ]]; then

	while IFS='' read -r line || [[ -n "$line" ]]; 
	do
		for i in ${Chromosomes_to_analyse[@]};
		do		    
			magma --annotate --snp-loc ${PRSACPJA_data_output} --gene-loc ${Pathway_output_directory}${training_set_name}_${validation_set_name}_${line}_chromosome_${Chromosomes_to_analyse[i]}_temp.gene.loc --out ${Pathway_output_directory}${Chromosomes_to_analyse[i]}_${training_set_name}_${validation_set_name}_SNPs_${line}_pathway_temp
		done

	done < "$1"
fi


Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Assign_SNPS_to_genes_from_pathways.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/assiging_SNPs_to_genes_from_pathways.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${valudation_set_usually_genotype_serial}\
 ${Name_of_extra_analysis}\
 ${Pathway_file_name}\
 ${Gene_location_filename}\
 ${Gene_regions}\
 ${Chromosomes_to_analyse[@]}

# magma --bfile CLOZUK_GWAS_BGE_chr22_magma_input_2 --gene-annot ${chr[i]}_CLOZUK_PGC_SNPs_pathway.genes.annot --out ${chr[i]}gene_annotation_for_CLOZUK_test

# Then write a new R script with previous details and finish, NOTE MAKE IT EASY TO DELETE USELESS FILES 

 
