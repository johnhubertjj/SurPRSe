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

#Am I Local or on a server? 
system=$3

if [[ "$whereami" = *"raven"* ]]; then
  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
  
  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11
  module load magma/1.06

  cd $SLURM_SUBMIT_DIR
  path_to_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
  path_to_pathway_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/Pathway_analysis_scripts"
  # Assign the shell variables
  source ${path_to_scripts}PRS_arguments_script.sh 
  cat ${path_to_scripts}PRS_arguments_script.sh 
   
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi
 
if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then
  
  Directory_to_work_from=$1
  cd ${Directory_to_work_from}
  
  # Arguments
  path_to_scripts=$2
  path_to_pathway_scripts=$4
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 
  
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi 
##########################################################################################################################
###################                  ###############                                 ###############                     #
################### PATHWAY ANALYSIS ###############      PATHWAY ANALYSIS           ############### PATHWAY ANALYSIS    #
##########################################################################################################################
if [[ ${Name_of_extra_analysis} = "Pathways" ]]; then
 
Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/" 

# Create output directory for Extra_analysis results
	if [ ! -d "./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}" ]; then
        	mkdir ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
	fi

sudo chmod  g+rwx ${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
sudo chmod  g+rwx ${training_set_name}_${validation_set_name}_output/
 
if [ -e "${Pathway_output_directory}Pathways_analysis_empty_pathways_info_file.txt" ]; then
	rm "${Pathway_output_directory}Pathways_analysis_empty_pathways_info_file.txt"
fi

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}PATHWAYS_PRS_COLLECTING_MAGMA_INFO.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/PATHWAYS_PRS_COLLECTING_MAGMA_INFO.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${validation_set_usually_genotype_serial}\
 ${Name_of_extra_analysis}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${path_to_stationary_data}${Gene_location_filename}\
 ${Chromosomes_to_analyse[@]}
 
 
# From the above script, identify the number of pathways you want to analyse (probably safest to write to a file, port to a variable and then delete the file)
# Also a text-delimited file with each line specifying a pathway name to be used
# The seperate gene_loc files belonging to previously specified analysis
# Then use magma to create annotation files:

Pathway_file_name=${Pathway_output_directory}Pathway_names.txt
# separate script below:

# Annotate using MAGMA here, separate out between extended and normal gene regions.
# Essentially, this will read an input file with the unique pathway names created previously and the gene-loc file\
# as an input from PATHWAYS_PRS_COLLECTING_MAGMA_INFO.R to annotate.
# The pathways used will then be read into an array variable so that we don't have to keep reading this file and it is within the BASH environment 

if [[ ${Gene_regions} = "Extended" ]]; then

	while IFS='' read -r line || [[ -n "$line" ]]; 
	do
		for i in ${Chromosomes_to_analyse[@]};
		do		    
			
			PRSACPJA_data_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}${i}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates.bim"
			magma --annotate window=35,10 --snp-loc ${PRSACPJA_data_output} --gene-loc ${Pathway_output_directory}${training_set_name}_${validation_set_name}_${line}_chromosome_${Chromosomes_to_analyse[i]}_extended_temp.gene.loc --out ${Pathway_output_directory}${Chromosomes_to_analyse[i]}_${training_set_name}_${validation_set_name}_SNPs_${line}_extended_pathway_temp
		done
	pathways+=("$line")
	done < "$Pathway_file_name"
fi

if [[ "${Gene_regions}" = "normal" ]]; then
echo normal gene regions
echo ${Pathway_file_name}
	while IFS='' read -r line || [[ -n "$line" ]]; 
	do
		for i in ${Chromosomes_to_analyse[@]};
		do		    
			PRSACPJA_data_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}${i}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates.bim"
			magma --annotate --snp-loc ${PRSACPJA_data_output} --gene-loc ${Pathway_output_directory}${training_set_name}_${validation_set_name}_${line}_chromosome_${i}_temp.gene.loc --out ${Pathway_output_directory}${i}_${training_set_name}_${validation_set_name}_SNPs_${line}_pathway_temp
		done

	done <"$Pathway_file_name"
fi

while IFS='' read -r line; do pathways+=("$line"); done <$Pathway_file_name

if [ -e "${Pathway_output_directory}Pathways_analysis_empty_pathways_info_file_run2.txt" ]; then
	rm "${Pathway_output_directory}Pathways_analysis_empty_pathways_info_file_run2.txt"
fi

if [ -e "${Pathway_output_directory}MAGMA_empty_files_after_analysis.txt" ]; then
	rm "${Pathway_output_directory}MAGMA_empty_files_after_analysis.txt"
fi


Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Assign_SNPs_to_genes_from_pathways.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${pathways}_assiging_SNPs_to_genes_from_pathways.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${validation_set_usually_genotype_serial}\
 ${Name_of_extra_analysis}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${path_to_stationary_data}${Gene_location_filename}\
 ${Gene_regions}\
 ${Chromosomes_to_analyse[@]}


# merge all the PGC SNPs together into one table
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}combining_summary_stats_tables_after_conversion_to_CHR_POS.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_conversion.Rout ${training_set_name} ${validation_set_name} ${Chromosomes_to_analyse[@]}

chmod -R g+rwx ${path_to_scripts} 
echo ${pathways[@]}

if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then

sudo parallel ${path_to_pathway_scripts}creation_of_merge_list_file.sh ::: ${pathways[@]} ::: ${path_to_scripts} ::: ${path_to_pathway_scripts} ::: ${system}

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Pathway_PRS_scoring.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Pathway_PRS_scoring.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${sig_thresholds[@]}

sudo parallel ${path_to_pathway_scripts}PRS_scoring_plink_pathways.sh ::: ${pathways[@]} ::: ${path_to_scripts} ::: ${path_to_pathway_scripts} ::: ${system}

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Collate_all_pathways.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Collate_all_pathways.Rout\ 
 ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${sig_thresholds}

# now just require the collate all paths script here... (can do at home)
else
	exit 1
fi
fi
# ADD ALEX's PIPELINE SCRIPTS HERE FOR PRS ANALYSIS AND TO GET ALL THE PATHWAYS IN THE RIGHT FORMAT IN PERL (but only the collate_all_paths.pl)

# magma --bfile CLOZUK_GWAS_BGE_chr22_magma_input_2 --gene-annot ${chr[i]}_CLOZUK_PGC_SNPs_pathway.genes.annot --out ${chr[i]}gene_annotation_for_CLOZUK_test

# Then write a new R script with previous details and finish, NOTE MAKE IT EASY TO DELETE USELESS FILES 



