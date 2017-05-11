# Arguments to be read into each analysis
# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" == *"raven"* ]]; then
  # State paths to the relevant stationary folders required for the analysis
  path_to_PGC_conversion="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="/home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/MAGMA/"
  number_of_files=($(find -E . -type f -regex '^./${training_set_name}_${validation_set_name}_output/CLOZUK_GWAS_BGE_CLUMPED_chr[0-9]+.bed' -exec basename {} \;))
  path_to_covariate_file="/home/c1020109/NCBI37.3/CLOZUK2.r7.select2PC.eigenvec.txt" 
  path_to_chromosome_length_file="/home/c1020109/NCBI37.3/UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="/home/c1020109/NCBI37.3/CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="/home/c1020109/NCBI37.3/NCBI37.3.gene.loc"  
  
  # assign arguments here for now because there are so many
  # Datasets
  training_set_usually_summary="PGC_BIP_table${chromosome_number}"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr${chromosome_number}"
  training_set_original_filename="daner_PGC_BIP32b_mds7a"
  training_set_name="PGC_BIP"
  validation_set_name="CLOZUK"
  # MAF, INFO and SE
  MAF_summary="FALSE"
  MAF_threshold=0.01
  MAF_genotype="TRUE"
  INFO_summary="TRUE"
  INFO_threshold=0.9
  SE_summary="FALSE"
  SE_threshold=5	
  Raven_out_info_directory="/Users/johnhubert/Documents/John_CLOZUK_PGC_PRS/CLOZUK_PRS_clumping_log_files"
  # The number of chromosomes you wish to analyse (PRS_serial)
  Chromosomes_to_analyse=(`seq 1 22`) 
  # Clumping Arguments
  p1=0.5
  p2=0.5
  r2=0.2
  window=1000
  # PRS_serial arguments
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  sig_thresholds=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
  Perform_Magma_as_well="TRUE"
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates" 
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions= "both"
  external_harddrive="FALSE"


elif [ "$whereami" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then

  path_to_PGC_conversion="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation"
  path_to_CLOZUK_conversion="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation"
  path_to_covariate_file="/Users/Dropbox/whole_genome_testing/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt"
  path_to_chromosome_length_file="/Users/Dropbox/whole_genome_testing/Stationary_data/UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="/Users/Dropbox/whole_genome_testing/output/Stationary_data/CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="/Users/Dropbox/whole_genome_testing/output/Stationary_data/NCBI37.3.gene.loc"

  chromosome_number=14
  # Datasets
  training_set_usually_summary="Neurot_Assoc_Biobank_table${chromosome_number}"
  training_set_original_filename="Neurot_Assoc_Biobank_PCA_Imputed_91370_8cov_array_noDupl_miss0.05_info0.4_12.01.16.res"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr${chromosome_number}"
  validation_set_usually_genotype_serial="CLOZUK_GWAS_BGE_chr"
  training_set_name="Neurot_Assoc_Biobank"
  validation_set_name="CLOZUK" 
  # MAF, INFO and SE
  MAF_summary="FALSE"
  MAF_threshold=0.01
  MAF_genotype="TRUE"
  INFO_summary="TRUE"
  INFO_threshold=0.9
  SE_summary="FALSE"
  SE_threshold=5
  Raven_out_info_directory="/Volumes/PhD_storage/${training_set_name}_${validation_set_name}_clumping_log_files/"
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
  sig_thresholds=(0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5)
  Perform_Magma_as_well="FALSE"
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions="both"	
  external_harddrive="FALSE"
fi

if [ "$external_harddrive" == "TRUE" ]; then
  path_to_harddrive=/Volumes/HD-PCU2
  cp $path_to_harddrive/CLOZUK_data/${validation_set_usually_genotype}.tar.gz .
  cp $path_to_harddrive/PGC_noCLOZUK_data/${training_set_usually_summary}.txt .
fi
