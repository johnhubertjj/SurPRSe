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
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_remove_MAF_INFO${chromosome_number}.Rout\

# From the above script, identify the number of pathways you want to analyse (probably safest to write to a file, port to a variable and then delete the file)
# Also a text-delimited file with each line specifying a pathway name to be used
# The seperate gene_loc files belonging to previously specified analysis
# Then use magma to create annotation files:
 
# separate script below:
#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Text read from file: $line"
done < "$1"

chmod +x rr.sh
./rr.sh filename.txt

#! /bin/bash

chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

pathways=("FMRP_targets" "abnormal_behavior" "abnormal_nervous_system_electrophysiology" "abnormal_learning|memory|conditioning"  "abnormal_CNS_synaptic_transmission" "Cav2_channels" "abnormal_synaptic_transmission" "5HT_2C" "abnormal_long_term_potentiation" "abnormal_motor_capabilities|coordination|movement" "abnormal_behavioral_response_to_xenobiotic" "abnormal_associative_learning" "Lek2015_LoFintolerant_90" "BGS_top2_mean" "BGS_top2_max")


for i in `seq 0 21` ;
do 
	for l in `seq 0 14`;
	do
		magma --annotate window=35,10 --snp-loc CLOZUK_GWAS_BGE_chr${chr[i]}_magma_input.bim --gene-loc ${pathways[l]}_chromosome_${chr[i]}.gene.loc --out ${chr[i]}_CLOZUK_PGC_SNPs_${pathways[l]}pathway
	done
done 

# magma --bfile CLOZUK_GWAS_BGE_chr22_magma_input_2 --gene-annot ${chr[i]}_CLOZUK_PGC_SNPs_pathway.genes.annot --out ${chr[i]}gene_annotation_for_CLOZUK_test

# Then write a new R script with previous details and finish, NOTE MAKE IT EASY TO DELETE USELESS FILES 

 
