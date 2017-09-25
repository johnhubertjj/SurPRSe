#!/bin/bash

exec &> Rtestoutclusterlogfile$1.txt
#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=24:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 1-22
#PBS -N PRS_whole_genome

# script requries 22 files for each validation and training set

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
system=$4

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
  
  # Assign the shell variables
  source ${path_to_scripts}PRS_arguments_script.sh 
  cat ${path_to_scripts}PRS_arguments_script.sh 
   
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi
 
if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then
  # add directory to work from at the top of the script (just in case)
  Directory_to_work_from=$2
  cd ${Directory_to_work_from} 
  
  #add the chromosome this copy of the script will analyse
  chromosome_number=$1
  echo $chromosome_number
  
  # path to where the scripts are located
  path_to_scripts=$3
  
  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh 
  
  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  echo ${validation_set_usually_genotype}
fi  

## rewrite so that the file input is an argument for the script instead, this will work for now
## -gt is greater than
## set essentially sets the arguments, so by nullglobbing everything, set will record the instances in which the pattern exists in the directory

shopt -s nullglob #enable
set -- *${chromosome_number}.tar.gz
if [ "$#" -gt 0 ]; then

# unpack the CLOZUK datasets
# remember to delete after unpacking at the end of the analysis
tar -zxvf ${validation_set_usually_genotype}.tar.gz
shopt -u nullglob # disable
fi

if [[ ${remove_IDs} = "TRUE" ]]; then

plink --bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype} --remove CLOZUK2_IDs_remove_plink_file_chr_${chromosome_number}.txt --make-bed --out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}

fi

# make directories for output and extra info
if [ ! -d "${training_set_name}_${validation_set_name}_output" ]; then
   mkdir ${training_set_name}_${validation_set_name}_output
fi

if [ ! -d "${training_set_name}_${validation_set_name}_extrainfo" ]; then
   mkdir ${training_set_name}_${validation_set_name}_extrainfo
fi
 
# Run R script that removes SNPs based on INFO score and MAF
Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_scripts}MAF_and_INFO_score_summary_stats_script.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_remove_MAF_INFO${chromosome_number}.Rout\
 ${training_set_usually_summary}\
 ${training_set_name}\
 ${validation_set_name}\
 ${MAF_summary}\
 ${MAF_threshold}\
 ${INFO_summary}\
 ${INFO_threshold}\
 ${SE_summary}\
 ${SE_threshold}\
 ${chromosome_number} 

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_scripts}CLOZUK_PGC_COMBINE_final.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${validation_set_name}_${training_set_name}_COMBINE_chr${chromosome_number}.Rout\
 ${training_set_usually_summary}\
 ${validation_set_usually_genotype}\
 ${training_set_name}\
 ${validation_set_name}\
 ${chromosome_number}\
 ${CHR_name}\
 ${SNP_name}\
 ${BP_name}\
 ${A1_name}\
 ${A2_name}\
 ${OR_name}\
 ${BETA_name}\
 ${system}\
 ${Number_of_frequency_columns}  

if [[ ${MAF_genotype} = "TRUE" ]]; then
if [[ ${Missing_geno} = "TRUE" ]]; then
     	# using plink to change the names to a CHR.POS identifier and remaking the files
  plink \
--bfile ${validation_set_usually_genotype} \
--update-name ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt \
--maf ${MAF_threshold} \
--geno ${genotype_missingness_check} \
--make-bed  \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2
else 

   # using plink to change the names to a CHR.POS identifier and remaking the files
  plink \
--bfile ${validation_set_usually_genotype} \
--update-name ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt \
--maf ${MAF_threshold} \
--make-bed  \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2
fi

elif [[ ${MAF_genotype} = "FALSE" ]]; then 
if [[ ${Missing_geno} = "TRUE" ]]; then

plink \
--bfile ${validation_set_usually_genotype} \
--update-name ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt \
--geno ${genotype_missingness_check} \
--make-bed \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2

else

   # using plink to change the names to a CHR.POS identifier and remaking the files
  plink \
--bfile ${validation_set_usually_genotype} \
--update-name ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt \
--make-bed  \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2
fi
fi

# re-package the original files
# tar -czvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz\
# CLOZUK_GWAS_BGE_chr${chromosome_number}.bed\
# CLOZUK_GWAS_BGE_chr${chromosome_number}.bim\
# CLOZUK_GWAS_BGE_chr${chromosome_number}.fam

if [[ ${Extra_analyses} = "TRUE" ]]; then
# NEED A CHECK ON WHAT EXTRA ANALYSES YOU ARE RUNNING HERE! 
Name_of_extra_analysis=(`sed -n 's/Name_of_extra_analysis\=//p' ${path_to_scripts}PRS_arguments_script.sh`)
length_of_extra_analysis_array=`echo ${#Name_of_extra_analysis[@]}`

if [ "${length_of_extra_analysis_array}" -gt  "2" ]; then
Name_of_extra_analysis_specific=(Pathways Genes)
# Create output directory for Extra_analysis results
	if [ ! -d "./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis_specific[0]}" ]; then
		mkdir ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis_specific[0]}
	fi

# Create files for Clumping
	plink \
--bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2 \
--exclude ./${training_set_name}_${validation_set_name}_output/extracted_Duplicate_snps_${validation_set_name}_${training_set_name}_chr${chromosome_number}.txt \
--extract ./${training_set_name}_${validation_set_name}_output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt \
--make-bed \
--out ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis_specific[0]}/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates

echo "Arguments script stated that this analysis is for ${Name_of_extra_analysis_specific[0]}\
, therefore clumping will be not be performed for Polygenic risk scores"

# Create output directory for Extra_analysis results
	if [ ! -d "./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis_specific[1]}" ]; then
		mkdir ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis_specific[1]}
	fi

# Create files for Clumping
	plink \
--bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2 \
--exclude ./${training_set_name}_${validation_set_name}_output/extracted_Duplicate_snps_${validation_set_name}_${training_set_name}_chr${chromosome_number}.txt \
--extract ./${training_set_name}_${validation_set_name}_output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt \
--make-bed \
--out ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis_specific[1]}/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates

echo "Arguments script stated that this analysis is for ${Name_of_extra_analysis_specific[1]}\
, therefore clumping will be not be performed for Polygenic risk scores"


else

	if [ ! -d "./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis[1]}" ]; then
		mkdir ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis[1]}
	fi

# Create files for Clumping
	plink \
--bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2 \
--exclude ./${training_set_name}_${validation_set_name}_output/extracted_Duplicate_snps_${validation_set_name}_${training_set_name}_chr${chromosome_number}.txt \
--extract ./${training_set_name}_${validation_set_name}_output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt \
--make-bed \
--out ./${training_set_name}_${validation_set_name}_output/${Name_of_extra_analysis[1]}/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates

echo "Arguments script stated that this analysis is for ${Name_of_extra_analysis[1]}\
, therefore clumping will be not be performed for Polygenic risk scores"


fi

# Exit so that clumping is not performed on the whole dataset
exit 0

else
	
	plink \
--bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2 \
--exclude ./${training_set_name}_${validation_set_name}_output/extracted_Duplicate_snps_${validation_set_name}_${training_set_name}_chr${chromosome_number}.txt \
--extract ./${training_set_name}_${validation_set_name}_output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt \
--make-bed \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates

fi

## Limit to genic SNPs here 
##########################
##########################

# Clump the datasets
# Extract the SNPs common between PGC and CLOZUK
# Remove the duplicate SNPs
plink \
--bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates \
--extract ./${training_set_name}_${validation_set_name}_output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt \
--exclude ./${training_set_name}_${validation_set_name}_output/extracted_Duplicate_snps_${validation_set_name}_${training_set_name}_chr${chromosome_number}.txt \
--clump ./${training_set_name}_${validation_set_name}_output/${training_set_name}_table${chromosome_number}_new.txt \
--clump-kb ${window} \
--clump-p1 ${p1} \
--clump-p2 ${p2} \
--clump-r2 ${r2} \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${chromosome_number}
 
# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' \
< ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${chromosome_number}.clumped \
> ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_CLUMPED_chr${chromosome_number}.txt

cut -f 2,4,5,6 \
< ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_CLUMPED_chr${chromosome_number}.txt \
> ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_chr${chromosome_number}.txt

rm ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_CLUMPED_chr${chromosome_number}.txt

awk '{ print $2 }' \
./${training_set_name}_${validation_set_name}_output/${validation_set_name}_${training_set_name}_CLUMPED_FINAL_chr${chromosome_number}.txt \
> ./${training_set_name}_${validation_set_name}_output/CLUMPED_EXTRACT_${validation_set_name}_chr${chromosome_number}.txt

printf "%s\n\n" \
"$(tail -n +2 ./${training_set_name}_${validation_set_name}_output/CLUMPED_EXTRACT_${validation_set_name}_chr${chromosome_number}.txt)" \
> ./${training_set_name}_${validation_set_name}_output/CLUMPED_EXTRACT_${validation_set_name}_chr${chromosome_number}.txt 

	
# Create clumped plink files
plink \
--bfile ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates \
--extract ./${training_set_name}_${validation_set_name}_output/CLUMPED_EXTRACT_${validation_set_name}_chr${chromosome_number}.txt \
--make-bed \
--out ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_${training_set_name}_CLUMPED

# Clean up original datasets to only leave the tar.gz file 
shopt -s nullglob
set -- *${validation_set_usually_genotype}*
if [ "$#" -gt 3 ]; then
	rm -rf ${validation_set_usually_genotype}.{bim,bed,fam,log} 
	shopt -u nullglob
fi

# Clean up other unneeded files
rm ./${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype}_2.{bim,bed,fam,log}
rm ./${training_set_name}_${validation_set_name}_output/${validation_set_name}_chr${chromosome_number}_chr.pos.txt
rm ./${training_set_name}_${validation_set_name}_output/chr${chromosome_number}${training_set_name}_${validation_set_name}_common_SNPs.txt
rm ./${training_set_name}_${validation_set_name}_output/extracted_Duplicate_snps_${validation_set_name}_${training_set_name}_chr${chromosome_number}.txt 
rm ./${training_set_name}_${validation_set_name}_output/CLUMPED_EXTRACT_${validation_set_name}_chr${chromosome_number}.txt 

# Append job ID to PRS_arguments_script
if [[ "$PBS_ARRAYID" = 1 ]]; then
	echo "Batch_job_ID=${PBS_JOBID}" >> ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
fi

#purge all modules 
if [[ "$whereami" = *"raven"* ]]; then
	module purge
fi

