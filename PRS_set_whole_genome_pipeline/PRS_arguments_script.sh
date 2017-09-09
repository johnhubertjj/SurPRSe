# Arguments to be read into each analysis
# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" = *"raven"* ]]; then
  # State paths to the relevant stationary folders required for the analysis
  path_to_PGC_conversion="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/MAGMA/"
  path_to_covariate_file="/home/c1020109/NCBI37.3/CLOZUK2.r7.select2PC.eigenvec.txt" 
  path_to_chromosome_length_file="/home/c1020109/NCBI37.3/UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="/home/c1020109/NCBI37.3/CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="/home/c1020109/NCBI37.3/NCBI37.3.gene.loc"  
  
  # assign arguments here for now because there are so many
  # Datasets
  training_set_usually_summary="CLOZUK_PGC2noclo_table${chromosome_number}"
  training_set_original_filename="daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt"
  validation_set_usually_genotype="ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr${chromosome_number}"
  validation_set_usually_genotype_serial="ALSPAC_hrc_imputed_step3_mri_brain_measurements_only_chr"
  validation_set_full_name_without_chromosome="ALSPAC_hrc_imputed_step3_mri_brain_measurements_only"
  training_set_name="CLOZUK_PGC2noclo"
  validation_set_name="PGCnoCLOZUK"  
  Pathway_filename="Pocklington2015_134sets_LoFi.txt"
  Gene_location_filename="NCBI37.3.gene.loc"
  # Split_by_chromosome for genotype?
  split_by_chromosome_required="FALSE" 
  # MAF, INFO and SE
  MAF_summary="FALSE"
  MAF_threshold=0.1
  MAF_genotype="TRUE"
  INFO_summary="TRUE"
  INFO_threshold=0.9
  SE_summary="FALSE"
  SE_threshold=5	
  Raven_out_info_directory="/home/c1020109/${training_set_name}_${validation_set_name}_clumping_log_files/"
  # The number of chromosomes you wish to analyse (PRS_serial)
  Chromosomes_to_analyse=(`seq 1 22`) 
  # Clumping Arguments
  p1=1
  p2=1
  r2=0.1
  window=500
  # PRS_serial arguments
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  sig_thresholds=(5e-08 1e-06 1e-04 0.01 0.05 0.1 0.2 0.5 1)
  Extra_analyses=TRUE
  Name_of_extra_analysis=Pathways
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates" 
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions="both" #either ( "extended" "normal" "both" )
  external_harddrive="FALSE"


elif [ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
home_OS="/Users"
elif [ "$whereami" = 'johnhubert-ThinkPad-P50' ]; then
home_OS="/home"
fi

if [[ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' || "$whereami" = 'johnhubert-ThinkPad-P50' ]]; then

  path_to_PGC_conversion="Summary_stat_manipulation"
  path_to_CLOZUK_conversion="Genotype_dataset_manipulation"
  path_to_covariate_file="$home_OS/johnhubert/Dropbox/whole_genome_testing/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt"
  path_to_chromosome_length_file="$home_OS/johnhubert/Dropbox/whole_genome_testing/Stationary_data/UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="$home_OS/johnhubert/Dropbox/whole_genome_testing/output/Stationary_data/CLOZUK.r7.GWAS_IDs.fam"
  path_to_stationary_data="$home_OS/johnhubert/Dropbox/Stationary_data/" 
  
  # Datasets
  training_set_usually_summary="IQ_GWAS_2017_table${chromosome_number}"
  training_set_original_filename="IQ_GWAS_2017_why_have_you_used_stupid_headings_sniekersetal.txt"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr${chromosome_number}"
  validation_set_usually_genotype_serial="CLOZUK_GWAS_BGE_chr"
  validation_set_full_name_without_chromosome="CLOZUK_GWAS_BGE"
  training_set_name="IQ_GWAS_2017"
  validation_set_name="CLOZUK" 
  
  # Pathway datasets
  Pathway_filename="Pocklington2015_134sets_LoFi_2sets_morphology_notmorphology_deduplicated.txt"
  Gene_location_filename="NCBI37.3.gene.loc"
  
  # Split_by_chromosome for genotype?
  split_by_chromosome_required="FALSE" 
  # Do a missingness check? 
  Missing_geno=FALSE
  genotype_missingness_check=FALSE
  # MAF, INFO and SE
  MAF_summary="FALSE"
  MAF_threshold=0.01
  MAF_genotype=TRUE
  INFO_summary="TRUE"
  INFO_threshold=0.9
  SE_summary="FALSE"
  SE_threshold=5
  Raven_out_info_directory="${training_set_name}_${validation_set_name}_clumping_log_files/"
  # The number of chromosomes you wish to analyse (PRS_serial)
  Chromosomes_to_analyse=(`seq 1 22`)
  
  # Clumping Arguments
  p1=0.5
  p2=0.5
  r2=0.2
  window=1000
  
  # Arguments for PRS_serial script
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  sig_thresholds=(1e-04 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
  Extra_analyses=FALSE
  Name_of_extra_analysis=Pathways
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions=normal #either ( "extended" "normal" "both" )	
  external_harddrive="FALSE"
fi

if [ "$external_harddrive" = "TRUE" ]; then
  path_to_harddrive=/Volumes/HD-PCU2
  cp $path_to_harddrive/CLOZUK_data/${validation_set_usually_genotype}.tar.gz .
  cp $path_to_harddrive/PGC_noCLOZUK_data/${training_set_usually_summary}.txt .
fi
