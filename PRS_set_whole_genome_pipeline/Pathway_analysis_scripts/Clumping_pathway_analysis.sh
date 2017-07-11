#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=10:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 8-9:2
#PBS -N c1020109_job_array_pathways

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" == *"raven"* ]]; then
  path_to_scripts='/home/c1020109/PRS_scripts/'
  # Load both Plink and R
  module purge
  module load plink/1.9c3

  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
 
  num1=1
  current_pathway_output=`echo "$((${chromosome_number} - ${num1}))"`

  WDPATH=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC/pathways_PRS_Antonio_paper/${pathways[${current_pathway_output}]}/
  PATH_for_PGC=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC/
  
  cd $WDPATH

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
fi

PRSACPJA_data_output="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/${validation_set_usually_genotype_serial}${Chromosomes_to_analyse[@]}_consensus_with_    ${training_set_name}_flipped_alleles_no_duplicates"
Pathway_output_directory="./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis}/" 
Pathway_file_name=${Pathway_output_directory}Pathway_analysis.txt

while IFS='' read -r line || [[ -n "$line" ]];
do
	pathways+=("$line")
done < "$Pathway_file_name"


plink --bfile ${pathways[${current_pathway_output}]}_Clumping_input --clump ${PATH_for_PGC}combined_PGC_table_with_CHR.POS_identifiers.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000 --out ${pathways[${current_pathway_output}]}_CLOZUK_PGC_pathways_r0.2_removed_AT_CG_1000kb

# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < ${pathways[${current_pathway_output}]}_CLOZUK_PGC_pathways_r0.2_removed_AT_CG_1000kb.clumped > pathways_CLOZUK_PGC_CLUMPED_${pathways[${current_pathway_output}]}.txt
cut -f 2,4,5,6 < pathways_CLOZUK_PGC_CLUMPED_${pathways[${current_pathway_output}]}.txt > pathways_CLOZUK_PGC_CLUMPED_FINAL_${pathways[${current_pathway_output}]}.txt
rm pathways_CLOZUK_PGC_CLUMPED_${pathways[${current_pathway_output}]}.txt
awk '{ print $2 }' pathways_CLOZUK_PGC_CLUMPED_FINAL_${pathways[${current_pathway_output}]}.txt > pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt
printf "%s\n\n" "$(tail -n +2 pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt)" > pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt 

# Create clumped plink files
plink --bfile ${pathways[${current_pathway_output}]}_Clumping_input --extract pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt --make-bed --out pathways_CLOZUK_GWAS_BGE_CLUMPED_${pathways[${current_pathway_output}]}


