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

  cd $PBS_O_WORKDIR
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
  path_to_gene_scripts=$4
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 
  
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi 

##########################################################################################################################
###################                  ###############                                 ###############                     #
################### GENE ANALYSIS    ###############      GENE ANALYSIS              ############### GENE    ANALYSIS    #
##########################################################################################################################

# Read in arguments to check what type of analysis you will need to do #
Name_of_extra_analysis=$5

# Check that the analysis is indeed for Genes
if [[ ${Name_of_extra_analysis} = "Genes" ]]; then
 
Gene_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/" 

# Create output directory for Extra_analysis results
	if [ ! -d "./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}" ]; then
        	mkdir ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
	fi

sudo chmod  g+rwx ${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}
sudo chmod  g+rwx ${training_set_name}_${validation_set_name}_output/
sudo chmod  g+rwx ${path_to_gene_scripts}Genes_MAGMA_annotation_script.sh 

if [ -e "${Pathway_output_directory}Pathways_analysis_empty_pathways_info_file.txt" ]; then
	rm "${Pathway_output_directory}Pathways_analysis_empty_pathways_info_file.txt"
fi

#run annotations of SNPs to genes, include useful information files
sudo parallel ${path_to_gene_scripts}Genes_MAGMA_annotation_script.sh ::: ${Chromosomes_to_analyse[@]} ::: ${Directory_to_work_from} ::: ${path_to_scripts} ::: ${path_to_gene_scripts} ::: ${system} 


# From the above script, identify the number of pathways you want to analyse (probably safest to write to a file, port to a variable and then delete the file)
# Also a text-delimited file with each line specifying a pathway name to be used
# The seperate gene_loc files belonging to previously specified analysis
# Then use magma to create annotation files:

# separate script below:

# Annotate using MAGMA here, separate out between extended and normal gene regions.
# Essentially, this will read an input file with the unique pathway names created previously and the gene-loc file\
# as an input from PATHWAYS_PRS_COLLECTING_MAGMA_INFO.R to annotate.
# The pathways used will then be read into an array variable so that we don't have to keep reading this file and it is within the BASH environment 

# merge all the PGC SNPs together into one table
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}combining_summary_stats_tables_after_conversion_to_CHR_POS.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_conversion.Rout ${training_set_name} ${validation_set_name} ${Chromosomes_to_analyse[@]}

chmod -R g+rwx ${path_to_scripts} 

if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then

${path_to_gene_scripts}creation_of_merge_list_file.sh ${Directory_to_work_from} ${path_to_scripts} ${path_to_gene_scripts} ${system}

# re-annotate the large clumped file again (unfortunately, might be unneeded, but might need to be quick

if [[ ${Gene_regions} = "both" ]]; then

Gene_regions=normal

	gene_bim_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final.bim
        gene_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final
	ending_name=${Gene_regions}_gene_regions_Clumped_whole_genome_final
        
	magma --annotate --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${Gene_output_directory}${training_set_name}_${validation_set_name}_SNPs_normal_clumped_gene_temp
	
ldsc.py --bfile ${gene_file}\ 
 --l2\
 --ld-wind-kb ${window}\
 --out ${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final

plink --bfile ${gene_file}\
 --freq\
 --out ${gene_file}

if [[ ${whole_genome_genic} = "TRUE" ]]; then

# Create Final PRS outputs directory
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring" ]; then
        mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring
fi

# Create score files for each significance threshold specified 
Rscript ${path_to_scripts}RscriptEcho.R\
       	${path_to_scripts}PRS_scoring_whole_genome.R\
       	./${training_set_name}_${validation_set_name}_extrainfo/PRS_scoring_whole_genome.Rout\
       	${training_set_name}\
       	${validation_set_name}\
       	${gene_bim_file}\
        ${ending_name}\
	${sig_thresholds[@]}

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
        filename="./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}.score"
        plink --bfile ${gene_file} --score ${filename} --out ./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}
done

fi

Gene_regions=extended

gene_bim_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final.bim
gene_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final
ending_name=${Gene_regions}_gene_regions_Clumped_whole_genome_final
        
	magma --annotate window=35,10 --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${Gene_output_directory}${training_set_name}_${validation_set_name}_SNPs_extended_clumped_gene_temp

ldsc.py --bfile ${gene_file}\ 
 --l2\
 --ld-wind-kb ${window}\
 --out ${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final

plink --bfile ${gene_file}\
 --freq\
 --out ${gene_file}

if [[ ${whole_genome_genic} = "TRUE" ]]; then

# Create Final PRS outputs directory
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring" ]; then
        mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring
fi

# Create score files for each significance threshold specified 
Rscript ${path_to_scripts}RscriptEcho.R\
       	${path_to_scripts}PRS_scoring_whole_genome.R\
       	./${training_set_name}_${validation_set_name}_extrainfo/PRS_scoring_whole_genome.Rout\
       	${training_set_name}\
       	${validation_set_name}\
       	${gene_bim_file}\
        ${ending_name}\
	${sig_thresholds[@]}

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
        filename="./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}.score"
        plink --bfile ${gene_file} --score ${filename} --out ./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}
done

fi
Gene_regions=both

elif [[ ${Gene_regions} = "extended" ]]; then
	
gene_bim_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final.bim
gene_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final
ending_name=${Gene_regions}_gene_regions_Clumped_whole_genome_final
	
	magma --annotate window=35,10 --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${Gene_output_directory}${training_set_name}_${validation_set_name}_SNPs_extended_clumped_gene_temp

ldsc.py --bfile ${gene_file}\ 
 --l2\
 --ld-wind-kb ${window}\
 --out ${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final

plink --bfile ${gene_file}\
 --freq\
 --out ${gene_file}

if [[ ${whole_genome_genic} = "TRUE" ]]; then

# Create Final PRS outputs directory
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring" ]; then
        mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring
fi

# Create score files for each significance threshold specified 
Rscript ${path_to_scripts}RscriptEcho.R\
       	${path_to_scripts}PRS_scoring_whole_genome.R\
       	./${training_set_name}_${validation_set_name}_extrainfo/PRS_scoring_whole_genome.Rout\
       	${training_set_name}\
       	${validation_set_name}\
       	${gene_bim_file}\
        ${ending_name}\
	${sig_thresholds[@]}

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
        filename="./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}.score"
        plink --bfile ${gene_file} --score ${filename} --out ./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}
done

fi
elif [[ ${Gene_regions} = "normal" ]]; then

gene_bim_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final.bim
gene_file=${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final
ending_name=${Gene_regions}_gene_regions_Clumped_whole_genome_final

	 magma --annotate --snp-loc ${gene_bim_file} --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${Gene_output_directory}${training_set_name}_${validation_set_name}_SNPs_normal_clumped_gene_temp

ldsc.py --bfile ${gene_file}\ 
 --l2\
 --ld-wind-kb ${window}\
 --out ${Gene_output_directory}${validation_set_name}_${training_set_name}_${Gene_regions}_gene_regions_Clumped_whole_genome_final

plink --bfile ${gene_file}\
 --freq\
 --out ${gene_file}

if [[ ${whole_genome_genic} = "TRUE" ]]; then

# Create Final PRS outputs directory
if [ ! -d "./${training_set_name}_${validation_set_name}_output/PRS_scoring" ]; then
        mkdir ./${training_set_name}_${validation_set_name}_output/PRS_scoring
fi

# Create score files for each significance threshold specified 
Rscript ${path_to_scripts}RscriptEcho.R\
       	${path_to_scripts}PRS_scoring_whole_genome.R\
       	./${training_set_name}_${validation_set_name}_extrainfo/PRS_scoring_whole_genome.Rout\
       	${training_set_name}\
       	${validation_set_name}\
       	${gene_bim_file}\
        ${ending_name}\
	${sig_thresholds[@]}

# Create PRS profiles for each significance threshold specified
for i in "${sig_thresholds[@]}" ;
do
        filename="./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}.score"
        plink --bfile ${gene_file} --score ${filename} --out ./${training_set_name}_${validation_set_name}_output/PRS_scoring/${training_set_name}_${validation_set_name}_${ending_name}_significance_threshold_at_${i}
done

fi

fi

if [[ ${Gene_specific_PRS} = "FALSE" ]]; then

exit 1
fi
 

mkdir ${Gene_output_directory}Genes_PRS/ 

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_Gene_scripts}Gene_specific_polygenic_risk_score_and_clumped_information_for_randomised_sets.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Pathway_PRS_scoring.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${validation_set_usually_genotype_serial}\
 ${Gene_output_directory}\
 ${path_to_stationary_data}${Gene_location_filename}\
 ${Gene_regions}\
 ${sig_thresholds[@]}


## Unfinished, below is similar procedure for pathways ## 
exit 1

Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_pathway_scripts}Collate_all_pathways.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_${validation_set_name}_Collate_all_pathways.Rout\ 
 ${training_set_name}\
 ${validation_set_name}\
 ${Pathway_output_directory}\
 ${path_to_stationary_data}${Pathway_filename}\
 ${sig_thresholds[@]}

# Run Magma Gene-set analysis for comparison to PRS if required
 
# magma --bfile ~/Desktop/ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr20_consensus_with_CLOZUK_PGC2noclo_more_sig_thresh_flipped_alleles_no_duplicates --gene-annot ~/Desktop/test.genes.annot --out ~/Desktop/testing_magma_output

else
	exit 1
fi
fi
# ADD ALEX's PIPELINE SCRIPTS HERE FOR PRS ANALYSIS AND TO GET ALL THE PATHWAYS IN THE RIGHT FORMAT IN PERL (but only the collate_all_paths.pl)

# magma --bfile CLOZUK_GWAS_BGE_chr22_magma_input_2 --gene-annot ${chr[i]}_CLOZUK_PGC_SNPs_pathway.genes.annot --out ${chr[i]}gene_annotation_for_CLOZUK_test

# Then write a new R script with previous details and finish, NOTE MAKE IT EASY TO DELETE USELESS FILES 
