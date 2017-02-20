#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=10:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
<<<<<<< HEAD
#PBS -J 1-15
=======
#PBS -J 21-22:2
>>>>>>> f1453e0ee59bc2ed7d1165fcb816ea27845efcb2
#PBS -N c1020109_job_array_pathways

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" == *"raven"* ]]; then
  path_to_scripts='/home/c1020109/PRS_scripts/'
  # Load both Plink and R
  module purge
<<<<<<< HEAD
  module load plink/1.9c3
=======
  module load R/3.1.0
  module load plink/1.9c3
  module load python/2.7.9-mpi
>>>>>>> f1453e0ee59bc2ed7d1165fcb816ea27845efcb2

  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}
  
  pathways=("FMRP_targets" "abnormal_behavior" "abnormal_nervous_system_electrophysiology" "abnormal_learning|memory|conditioning"  "abnormal_CNS_synaptic_transmission" "Cav2_channels" "abnormal_synaptic_transmission" "5HT_2C" "abnormal_long_term_potentiation" "abnormal_motor_capabilities|coordination|movement" "abnormal_behavioral_response_to_xenobiotic" "abnormal_associative_learning" "Lek2015_LoFintolerant_90" "BGS_top2_mean" "BGS_top2_max")

  num1=1
  current_pathway_output=`echo "$((${chromosome_number} - ${num1}))"`

<<<<<<< HEAD
  WDPATH=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC/pathways_PRS_Antonio_paper/${pathways[${current_pathway_output}]}/
  PATH_for_PGC=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC/
=======
  WDPATH=/scratch/$USER/PR54/PGC_CLOZUK_PRS/PRS_CLOZUK_PGC/pathways_PRS_Antonio_paper/output/${pathways[${current_pathway_output}]}/
  
>>>>>>> f1453e0ee59bc2ed7d1165fcb816ea27845efcb2
  cd $WDPATH

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  chromosome_number=22
fi



<<<<<<< HEAD
plink --bfile ${pathways[${current_pathway_output}]}_Clumping_input --clump ${PATH_for_PGC}combined_PGC_table_with_CHR.POS_identifiers.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000 --out ${pathways[${current_pathway_output}]}_CLOZUK_PGC_pathways_r0.2_removed_AT_CG_1000kb
=======
plink --bfile ${pathways[${current_pathway_output}]}_Clumping_input --clump combined_PGC_table_with_CHR.POS_identifiers.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000 --out ${pathways[${current_pathway_output}]}_CLOZUK_PGC_pathways_r0.2_removed_AT_CG_1000kb
>>>>>>> f1453e0ee59bc2ed7d1165fcb816ea27845efcb2

# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < ${pathways[${current_pathway_output}]}_CLOZUK_PGC_pathways_r0.2_removed_AT_CG_1000kb.clumped > pathways_CLOZUK_PGC_CLUMPED_${pathways[${current_pathway_output}]}.txt
cut -f 2,4,5,6 < pathways_CLOZUK_PGC_CLUMPED_${pathways[${current_pathway_output}]}.txt > pathways_CLOZUK_PGC_CLUMPED_FINAL_${pathways[${current_pathway_output}]}.txt
rm pathways_CLOZUK_PGC_CLUMPED_${pathways[${current_pathway_output}]}.txt
awk '{ print $2 }' pathways_CLOZUK_PGC_CLUMPED_FINAL_${pathways[${current_pathway_output}]}.txt > pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt
printf "%s\n\n" "$(tail -n +2 pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt)" > pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt 

# Create clumped plink files
plink --bfile ${pathways[${current_pathway_output}]}_Clumping_input --extract pathways_CLUMPED_EXTRACT_CLOZUK_${pathways[${current_pathway_output}]}.txt --make-bed --out pathways_CLOZUK_GWAS_BGE_CLUMPED_${pathways[${current_pathway_output}]}


