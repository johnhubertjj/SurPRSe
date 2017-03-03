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

