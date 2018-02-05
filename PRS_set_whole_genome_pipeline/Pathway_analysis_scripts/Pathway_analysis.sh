#/bin/bash

# script requries 22 files for each validation and training set

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

#Am I Local or on a server? 
system=$3
 
if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then
  
  Directory_to_work_from=$1
  cd ${Directory_to_work_from}
  
  # Arguments
  path_to_scripts=$2
  path_to_pathway_scripts=$4
  path_to_gene_scripts=$5

  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 
  
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi

Name_of_extra_analysis=$6 
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

if [ ${Using_raven} = "FALSE" ]; then
  sudo chmod  g+rwx ${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
  sudo chmod  g+rwx ${training_set_name}_${validation_set_name}_output/
 else 
  chmod  u+rwx ${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
  chmod  u+rwx ${training_set_name}_${validation_set_name}_output/
fi

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

# it stops here, run interactive session here to solve...

if [[ ${Gene_regions} = "extended" || ${Gene_regions} = "both" ]]; then
echo "extended gene regions"
	while IFS='' read -r line || [[ -n "$line" ]]; 
	do
		for i in ${Chromosomes_to_analyse[@]};
		do		    
			
			PRSACPJA_data_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}${i}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates.bim"
			magma --annotate window=35,10 --snp-loc ${PRSACPJA_data_output} --gene-loc ${Pathway_output_directory}${training_set_name}_${validation_set_name}_${line}_chromosome_${i}_temp.gene.loc --out ${Pathway_output_directory}${i}_${training_set_name}_${validation_set_name}_SNPs_${line}_extended_pathway_temp
		done
	done <"$Pathway_file_name"
fi

if [[ "${Gene_regions}" = "normal" || ${Gene_regions} = "both" ]]; then
echo "normal gene regions"
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

if [[ "${Gene_regions}" = "normal" || ${Gene_regions} = "both" ]]; then

	if [ -e "${Pathway_output_directory}MAGMA_empty_files_after_analysis_normal.txt" ]; then
		rm "${Pathway_output_directory}MAGMA_empty_files_after_analysis_normal.txt"
	fi
fi

if [[ "${Gene_regions}" = "extended" || ${Gene_regions} = "both" ]]; then

	if [ -e "${Pathway_output_directory}MAGMA_empty_files_after_analysis_extended.txt" ]; then
		rm "${Pathway_output_directory}MAGMA_empty_files_after_analysis_extended.txt"
	fi
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

if [ ${Using_raven} = "FALSE" ]; then
  chmod -R g+rwx ${path_to_scripts} 

else
  chmod -R u+rwx ${path_to_scripts}

fi

echo ${pathways[@]}

if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then

if [ ${Using_raven} = "FALSE" ]; then
sudo parallel ${path_to_pathway_scripts}creation_of_merge_list_file.sh ::: ${pathways[@]} ::: ${path_to_scripts} ::: ${path_to_pathway_scripts} ::: ${system} ::: ${Name_of_extra_analysis}

else
parallel ${path_to_pathway_scripts}creation_of_merge_list_file.sh ::: ${pathways[@]} ::: ${path_to_scripts} ::: ${path_to_pathway_scripts} ::: ${system} ::: ${Name_of_extra_analysis}

fi

# If randomising Gene_sets (aka ran Gene specific PRS first) then run this script:

if [[ "$randomise" = TRUE ]]; then
# start a snakemake file from here: should link everything up together properly #
# Set up extra arguments from Genes directory and for PRS scoring for randomised sets
Gene_output_directory="./${training_set_name}_${validation_set_name}_output/Genes/" 

mkdir ${Pathway_output_directory}Randomised_gene_sets_analysis/
Random_directory=${Pathway_output_directory}Randomised_gene_sets_analysis/

mkdir ${Pathway_output_directory}Randomised_gene_sets_analysis/Scores/
Random_scoring_directory=${Pathway_output_directory}Randomised_gene_sets_analysis/Scores/

# Got to here; need to figure out file paths...again...also need to un gunzip LD score files
Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}generate_random_andrews_script.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_generate_random_andrews_script.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${Gene_output_directory}\
 ${path_to_stationary_data}${Gene_location_filename}\
 ${Gene_regions}\
 ${permutations}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${calculate_indep_SNPs}\
 ${Random_directory}\
 ${Random_scoring_directory}\
 ${sample_replace}\
 ${pathways[@]}

pathways_for_randomisation=(`awk '{ print $1 }' ${Pathway_output_directory}${training_set_name}_${validation_set_name}_random_pathways_to_test.txt`)
path_to_config_file=${path_to_scripts}scoring_config_yaml.yaml

${path_to_scripts}config_yaml_writer.sh ${path_to_scripts} ${system} ${path_to_config_file} ${Pathway_output_directory} ${Gene_output_directory} ${Random_directory} ${Random_scoring_directory} 

cp ${path_to_config_file} ${PBS_O_WORKDIR}

if [[ ${Using_raven} = FALSE ]]; then
sudo parallel ${path_to_pathway_scripts}test_script_randomised_plink.sh ::: ${pathways_for_randomisation[@]} ::: ${path_to_scripts} ::: ${Name_of_extra_analysis} ::: ${Gene_output_directory} ::: ${Random_scoring_directory} ::: ${permutations} 

else
parallel ${path_to_pathway_scripts}test_script_randomised_plink.sh ::: ${pathways_for_randomisation[@]} ::: ${path_to_scripts} ::: ${Name_of_extra_analysis} ::: ${Gene_output_directory} ::: ${Random_scoring_directory} ::: ${permutations} 

fi

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_gene_scripts}Collate_all_pathways_random.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_collate_all_pathways_random.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Gene_output_directory}\
 ${Pathway_output_directory}\
 ${Random_scoring_directory}\
 ${path_to_stationary_data}${Gene_location_filename}\
 ${Gene_regions}\
 ${permutations}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${pathways[@]}
 
Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Pathway_PRS_scoring.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Pathway_PRS_scoring.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${Gene_output_directory}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${sig_thresholds[@]}

if [ ${Using_raven} = FALSE ]; then
sudo parallel ${path_to_pathway_scripts}PRS_scoring_plink_pathways.sh ::: ${pathways[@]} ::: ${path_to_scripts} ::: ${path_to_pathway_scripts} ::: ${system} ::: ${Name_of_extra_analysis}
else
parallel ${path_to_pathway_scripts}PRS_scoring_plink_pathways.sh ::: ${pathways[@]} ::: ${path_to_scripts} ::: ${path_to_pathway_scripts} ::: ${system} ::: ${Name_of_extra_analysis}
fi

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Collate_all_pathways.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Collate_all_pathways.Rout ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${sig_thresholds[@]}

# now just require the collate all paths script here... (can do at home)
else
	exit 1
fi
fi
fi

find "${Random_scoring_directory}" -type f -name '*.score' -maxdepth 1 >> "${Random_directory}score_files_to_compress_or_delete.txt"
tar -zcv -T "${Random_directory}score_files_to_compress_or_delete.txt" -f "${Random_directory}${training_set_name}_${validation_set_name}_tarball_for_scores.tar.gz"


# ADD ALEX's PIPELINE SCRIPTS HERE FOR PRS ANALYSIS AND TO GET ALL THE PATHWAYS IN THE RIGHT FORMAT IN PERL (but only the collate_all_paths.pl)

# magma --bfile CLOZUK_GWAS_BGE_chr22_magma_input_2 --gene-annot ${chr[i]}_CLOZUK_PGC_SNPs_pathway.genes.annot --out ${chr[i]}gene_annotation_for_CLOZUK_test

# Then write a new R script with previous details and finish, NOTE MAKE IT EASY TO DELETE USELESS FILES 

