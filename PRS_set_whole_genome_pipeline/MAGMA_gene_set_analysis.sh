# Create a SNPLOC file and convert the gene-set input into a GMT file (un-needed but worth a check for compatibility in any case)
path_to_geneset_file=${path_to_stationary_data}${Pathway_filename}

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}SNP_loc_creator_and_gmt_formatter.R  ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_MAGMA_snploc_creator.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${path_to_stationary_data}${Pathway_filename} 

Pathway_filename_minus_extension="${Pathway_filename%.*}"

# Create Make file from list of current files in directory
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}MAKE_LIST_file_creator_unclumped_files.R ./${training_set_name}_${validation_set_name}_extrainfo/MAKE_list_file_creator_before_clumping.Rout ${training_set_name} ${validation_set_name} ${validation_set_usually_genotype_serial} ${Chromosomes_to_analyse[@]}

plink --merge-list ${training_set_name}_${validation_set_name}_output/Make_file_${validation_set_name}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates_full_genome_before_clumping.txt --make-bed --out ${validation_set_name}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates_full_genome_before_clumping

# Begin MAGMA analyses
magma --annotate window=35,10 --snp-loc combined_${training_set_name}_${validation_set_name}_table_with_CHR.POS_identifiers_flipped_alleles_no_duplicates_MAF_INFOsnploc.txt --gene-loc ${path_to_stationary_data}${Gene_location_filename} --out ${training_set_name}_${validation_set_name}_annotation_MAGMA

magma --bfile ${validation_set_name}_consensus_with_${training_set_name}_flipped_alleles_no_duplicates_full_genome_before_clumping --pval combined_${training_set_name}_${validation_set_name}_table_with_CHR.POS_identifiers_flipped_alleles_no_duplicates_MAF_INFO.txt N=${number_of_samples_training} --gene-annot ${training_set_name}_${validation_set_name}_annotation_MAGMA.genes.annot --out ${training_set_name}_${validation_set_name}_MAGMA

magma --gene-results ${training_set_name}_${validation_set_name}_MAGMA.genes.raw --set-annot ${Pathway_filename_minus_extension}.gmt --model fwer=100000 alpha=0.05 --out ${training_set_name}_${validation_set_name}_gene_set_analysis
