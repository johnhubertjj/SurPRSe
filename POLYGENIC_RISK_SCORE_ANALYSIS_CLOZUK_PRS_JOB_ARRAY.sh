#!/bin/bash

#PBS -q batch_long
#PBS -P PR54
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime 10:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 21-22:2
#PBS -N c1020109_job_array_whole_genome

echo "Press CTRL+C to proceed."
trap "pkill -f 'sleep 1h'" INT
trap "set +x ; sleep 1h ; set -x" DEBUG

# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
if [[ "$whereami" == *"raven"* ]]; then
  cd $PBS_O_WORKDIR
  path_to_scripts='~/c1020109/PRS_scripts/'

  # Load both Plink and R
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.9-mpi

  # assign a new variable for the PBS_ARRAY_variable
  chromosome_number=${PBS_ARRAY_INDEX}

elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
  cd ~/Documents/testing_PRS_chromosome_22/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  chromosome_number=22
fi

## rewrite so that the file input is an argument for the script instead, this will work for now
if [ -f *.tar.gz ] && [ ! -f *.bim  ] && [ ! -f *.bed ] && [ !  -f *.fam ]; then
# unpack the CLOZUK datasets
tar -zxvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz
fi

# make directories for output and extra info
if [ ! -d "output" ]; then
   mkdir output
fi

if [ ! -d "extrainfo" ]; then
   mkdir extrainfo
fi

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
R CMD BATCH ${path_to_scripts}CLOZUK_PGC_COMBINE_final.R ./extrainfo/CLOZUK_PGC_COMBINE_chr${chromosome_number}.Rout
 
# using plink to change the names to a CHR.POS identifier and remaking the files
plink --bfile CLOZUK_GWAS_BGE_chr${chromosome_number} --update-name CLOZUK_chr${chromosome_number}_chr.pos.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2

#re-package the original files
#tar -czvf CLOZUK_GWAS_BGE_chr${chromosome_number}.tar.gz CLOZUK_GWAS_BGE_chr${chromosome_number}.bed CLOZUK_GWAS_BGE_chr${chromosome_number}.bim CLOZUK_GWAS_BGE_chr${chromosome_number}.fam

# create files containing duplicate SNPs
R CMD BATCH ${path_to_scripts}Clumping_CLOZUK_PGC.R ./extrainfo/CLOZUK_PGC_clumpinginfo_chr${chromosome_number}.Rout 

# Create files for MAGMA
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --exclude ./output/extracted_Duplicate_snps_chr${chromosome_number}.txt --extract ./output/chr${chromosome_number}PGC_CLOZUK_common_SNPs.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_chr22_magma_input


# Clump the datasets
# Extract the SNPs common between PGC and CLOZUK
# Remove the duplicate SNPs
# plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --extract ./output/chr${chromosome_number}PGC_CLOZUK_common_SNPs.txt --exclude ./output/extracted_Duplicate_snps_chr${chromosome_number}.txt --clump ./output/PGC_table${chromosome_number}_new.txt --clump-kb 1000 --clump-p1 0.5 --clump-p2 0.5 --clump-r2 0.2 --out ./output/CLOZUK_PGC_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${chromosome_number}
 
# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < ./output/CLOZUK_PGC_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${chromosome_number}.clumped > ./output/CLOZUK_PGC_CLUMPED_chr${chromosome_number}.txt
cut -f 2,4,5,6 < ./output/CLOZUK_PGC_CLUMPED_chr${chromosome_number}.txt > ./output/CLOZUK_PGC_CLUMPED_FINAL_chr${chromosome_number}.txt
rm ./output/CLOZUK_PGC_CLUMPED_chr${chromosome_number}.txt
awk '{ print $2 }' ./output/CLOZUK_PGC_CLUMPED_FINAL_chr${chromosome_number}.txt > ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt
printf "%s\n\n" "$(tail -n +2 ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt)" > ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt 

# Create clumped plink files
plink --bfile ./output/CLOZUK_GWAS_BGE_chr${chromosome_number}_2 --extract ./output/CLUMPED_EXTRACT_CLOZUK_chr${chromosome_number}.txt --make-bed --out ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}

# recode genotypes for input into Python
plink --bfile ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number} --recode A --out ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}
plink --bfile ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number} --freq --out ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}

# Get the MAF from CLOZUK and import into the PGC_summary_stats for PRS analysis
awk '{ $6=$7=$8=$10=$12=$13=$14=$15=$16=$17=""; print$0}' ./output/PGC_table${chromosome_number}_new.txt > ./output/PGC_table_for_python${chromosome_number}.txt

# match the SNP rows to the MAF rows using R script
R CMD BATCH ${path_to_scripts}Prepare_both_CLOZUK_AND_PGC_for_PRS_MAF_Genotype.R ./extrainfo/Prepare_both_CLOZUK_PGC_MAF_chr${chromosome_number}.Rout

# Final removal of headings for PRS analysis
awk 'NR>1' ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}.raw > ./output/CLOZUK_GWAS_BGE_CLUMPED_chr${chromosome_number}_no_head.raw 

# Create Final PRS outputs
if [ ! -d "./output/PRS_scoring" ]; then
   mkdir ./output/PRS_scoring
fi

# Needs the MAGMA script
sh ${path_to_scripts}PRS_with_magma.sh ${chromosome_number} 

# Run preparation for annotation file for python scripts
R CMD BATCH ${path_to_scripts}MAGMA_python_annotation_table_creator.R ./extrainfo/MAGMA_annotation_table_creator.Rout

# Run PRS python script # Loop through significance thresholds
# Change for Raven or Local
if [ "$whereami" == "raven13" ]; then
   python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py
else
   python ${path_to_scripts}PRS_scoring_parallel_clump_maf_JJ.py 
fi

#### alter below ####
# Calculate the PRS and profiles
# R CMD BATCH PRS_scoring_R_script.R ./extrainfo/PRS_SCORING_INFO_CLOZUK_PGC_chr${chromosome_number}.Rout

#significant=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)

#for i in `seq 1 ${#significant[@]}` ;
#do
#	plink --silent --bfile CLOZUK_GWAS_BGE_chr${chromosome_number} --exclude extracted_Duplicate_snps_chr${chromosome_number}.txt --score score/scoring_PGC_CLOZUK_chromosome_${chromosome_number}_${significant[i]}.score --out profiles/chr_${chromosome_number}_clump_r0.2_1000kb_${significant[i]}
	
#done





