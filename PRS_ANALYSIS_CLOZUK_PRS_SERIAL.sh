
input_file_name="./output/CLOZUK_GWAS_BGE_CLUMPED_chr"

# List the number of files matching the output from the previous Job array script "POLYGENIC_RISK_SCORE*.sh"
shopt -s nullglob
number_of_files=(*${input_file_name}*)
shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later

# Find the length of the array containing the names of the files
echo "${number_of_files[@]}"  # Note double-quotes to avoid extra parsing of funny characters in filenames
length_of_array=`echo "${#number_of_files[@]}"`

# Check for make_list_file
if [ -f ./output/CLOZUK_PGC_FULL_GENOME_MAKE_LIST.txt ]; then
rm ./${pathways[a]}/make_full_plink_${pathways[a]}.txt
fi

for chromosome_number in `seq 1 ${length_of_array}` ;
printf "${input_file_name}_${number_of_files[chromosome_number]}\n%.0s" {1} >> ./output/CLOZUK_PGC_FULL_GENOME_MAKE_LIST_INPUT.txt
done

# merge all the clumped chromosomes together
plink --merge-list ./output/CLOZUK_PGC_FULL_GENOME_MAKE_LIST.txt --make-bed --out ./output/CLOZUK_PGC_FULL_GENOME
 
# recode genotypes for input into Python
plink --bfile ./output/CLOZUK_PGC_FULL_GENOME --recode A --out ./output/CLOZUK_PGC_FULL_GENOME
plink --bfile ./output/CLOZUK_PGC_FULL_GENOME --freq --out ./output/CLOZUK_PGC_FULL_GENOME

# Get the MAF from CLOZUK and import into the PGC_summary_stats for PRS analysis
awk '{ $6=$7=$8=$10=$12=$13=$14=$15=$16=$17=""; print$0}' ./output/PGC_table${chromosome_number}_new.txt > ./output/PGC_table_for_python${chromosome_number}.txt

# match the SNP rows to the MAF rows using R script
R CMD BATCH ${path_to_scripts}Prepare_both_CLOZUK_AND_PGC_for_PRS_MAF_Genotype.R ./extrainfo/Prepare_both_CLOZUK_PGC_MAF_chr${chromosome_number}.Rout

# Final removal of headings for PRS analysis
awk 'NR>1' ./output/CLOZUK_PGC_FULL_GENOME.raw > ./output/CLOZUK_PGC_FULL_GENOME_no_head.raw

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
 
