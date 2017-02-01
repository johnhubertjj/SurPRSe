#!/bin/bash

#PBS -q batch-long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime 10:00:00
#PBS -o /home/c1020109/
#PBS -e /home/c1020109/
#PBS -j oe
#PBS -J 21-22:2
#PBS -N c1020109_job_array

cd $PBS_O_WORKDIR

# Load both Plink and R
module load R/3.3.0
module load plink/1.9a

# unpack the CLOZUK datasets
tar -zxvf CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX}.tar.gz

# make directories for output and extra info
if [ ! -d "output" ]; then
   mkdir output
fi

if [ ! -d "extrainfo" ]; then
   mkdir extrainfo
fi

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
R CMD BATCH CLOZUK_PGC_COMBINE_final.R ./extrainfo/Rout_files/CLOZUK_PGC_COMBINE_chr${PBS_ARRAY_INDEX}.Rout
 
# using plink to change the names to a CHR.POS identifier and remaking the files
plink 
 --bfile CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX} 
 --update-name CLOZUK_chr${PBS_ARRAY_INDEX}_chr.pos.txt 
 --make-bed 
 --out CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX}_2

# create files containing duplicate SNPs
R CMD BATCH Clumping_CLOZUK_PGC.R ./extrainfo/Rout_files/CLOZUK_PGC_clumpinginfo_chr${PBS_ARRAY_INDEX}.Rout 

# Clump the datasets
# Extract the SNPs common between PGC and CLOZUK
# Remove the duplicate SNPs
plink 
 --bfile CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX}_2
 --extract chr${PBS_ARRAY_INDEX}PGC_CLOZUK_common_SNPs.txt
 --exclude extracted_Duplicate_snps_chr${PBS_ARRAY_INDEX}
 --clump PGC_table${PBS_ARRAY_INDEX}.txt
 --clump-kb 1000 
 --clump-p1 0.5 
 --clump-p2 0.5 
 --clump-r2 0.2 
 --out CLOZUK_PGC_bge_removed_A.T_C.G.target_r0.2_1000kb_non_verbose_chr${PBS_ARRAY_INDEX}
 
# Clean up the files to leave a dataset that can be read into R/Python as well as a list of SNPs to extract for the CLUMPED plink files
tr -s ' ' '\t' < FT4.training_bge_LOC_removed_A.T_C.G.target_r0.2_1000kb.clumped_chr${PBS_ARRAY_INDEX} > CLOZUK_PGC_CLUMPED_chr${PBS_ARRAY_INDEX}.txt
cut -f 2,4,5,6 < CLOZUK_PGC_CLUMPED_chr${PBS_ARRAY_INDEX}.txt > CLOZUK_PGC_CLUMPED_FINAL_chr${PBS_ARRAY_INDEX}.txt
rm CLOZUK_PGC_CLUMPED_chr${PBS_ARRAY_INDEX}.txt
awk '{ print $2}' CLOZUK_PGC_CLUMPED_FINAL_chr${PBS_ARRAY_INDEX} > CLUMPED_EXTRACT_CLOZUK_chr${PBS_ARRAY_INDEX}.txt
printf "%s\n\n" "$(tail -n +2 CLUMPED_EXTRACT_CLOZUK_chr${PBS_ARRAY_INDEX}.txt)" > CLUMPED_EXTRACT_CLOZUK_chr${PBS_ARRAY_INDEX}.txt 

# Create clumped plink files
plink 
 --bfile CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX}_2 
 --extract CLUMPED_EXTRACT_CLOZUK_chr${PBS_ARRAY_INDEX}.txt 
 --make-bed 
 --out CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX}

# recode genotypes for input into Python
plink --bfile CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX} --recodeA --out CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX}
plink --bfile CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX} --freq --out CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX}

# Get the MAF from CLOZUK and import into the PGC_summary_stats for PRS analysis
awk '{ $6=$7=$8=$10=$12=$13=$14=$15=$16=$17=""; print$0}' PGC_table${PBS_ARRAY_INDEX}.txt > PGC_table_for_python${PBS_ARRAY_INDEX}.txt

# match the SNP rows to the MAF rows using R script
R CMD BATCH Prepare_both_CLOZUK_AND_PGC_for_PRS_MAF_Genotype.R

awk 'NR>1' CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX}.raw > CLOZUK_GWAS_BGE_CLUMPED_chr${PBS_ARRAY_INDEX}_no_head.raw 


#### alter below ####
# Calculate the PRS and profiles
R CMD BATCH PRS_scoring_R_script.R ./extrainfo/Rout_files/PRS_SCORING_INFO_CLOZUK_PGC_chr${PBS_ARRAY_INDEX}.Rout

significant=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)

for i in `seq 1 ${#significant[@]}` ;
do
	plink --silent --bfile CLOZUK_GWAS_BGE_chr${PBS_ARRAY_INDEX} --exclude extracted_Duplicate_snps_chr${PBS_ARRAY_INDEX}.txt --score score/scoring_PGC_CLOZUK_chromosome_${PBS_ARRAY_INDEX}_${significant[i]}.score --out profiles/chr_${PBS_ARRAY_INDEX}_clump_r0.2_1000kb_${significant[i]}
	
done





