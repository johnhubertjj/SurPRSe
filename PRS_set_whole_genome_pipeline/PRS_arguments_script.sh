# Arguments to be read into each analysis
# Run locally or on ARCCA

whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" = *"raven"* ]]; then
  
  #Path_to_locations_on_the_server
  path_to_PRS_scripts="${HOME}/PhD_scripts/Schizophrenia_pipeline_scipts/"
  # Re-assign to the training_set_usually_genotype and validation_full_name_without_chromosome using the sed command
  path_to_test_dataset="/scratch/${USER}/PR54/ALSPAC_test_pipeline/ALSPAC_hrc_imputed_step3_mri_brain_measurements_only"
  path_to_training_dataset="/scratch/${USER}/PR54/ALSPAC_test_pipeline/CLOZUK_PGC2noclo.METAL.assoc.dosage"

  # State paths to the relevant stationary folders required for the analysis
  path_to_PGC_conversion="${path_to_PRS_scripts}/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="${path_to_PRS_scripts}/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="${path_to_PRS_scripts}/MAGMA/"
  
  # Path to where stationary file is kept
  path_to_stationary_file="${HOME}/Stationary_data/"
  path_to_covariate_file="${path_to_stationary_data}CLOZUK2.r7.select2PC.eigenvec.txt" 
  path_to_chromosome_length_file="${path_to_stationary_data}UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="${path_to_stationary_data}CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="${path_to_stationary_file}NCBI37.3.gene.loc"  
  
  # assign arguments here for now because there are so many
  # Datasets
  training_set_name="CLOZUK_PGC2noclo"
  validation_set_name="ALSPAC_brain_measurements_only" 
  training_set_original_filename="CLOZUK_PGC2noclo.METAL.assoc.dosage"
  validation_set_full_name_without_chromosome="ALSPAC_hrc_imputed_step3_mri_brain_measurements_only"
  
  # DO NOT ALTER!!!
  validation_set_usually_genotype="${validation_set_full_name_without_chromosome}_chr${chromosome_number}"
  validation_set_usually_genotype_serial="${validation_set_full_name_without_chromosome}_chr"
  training_set_usually_summary="${training_set_name}_table${chromosome_number}"
  
  #Pathway datasets
  Pathway_filename="Pocklington2015_134sets_LoFi.txt"
  Gene_location_filename="NCBI37.3.gene.loc"
  
  # Split_by_chromosome for genotype?
  split_by_chromosome_required=FALSE 
  
  # Do a missingness check?
  Missing_geno=FALSE
  genotype_missingness_check=FALSE
  
  # MAF, INFO and SE
  MAF_summary=FALSE
  MAF_threshold=0.01
  MAF_genotype=TRUE
  INFO_summary=TRUE
  INFO_threshold=0.9
  SE_summary=FALSE
  SE_threshold=5	
  Raven_out_info_directory="${training_set_name}_${validation_set_name}_clumping_log_files/"
  
  # The number of chromosomes you wish to analyse (PRS_serial)
  Chromosomes_to_analyse=(`seq 1 22`) 
  
  # Clumping Arguments
  p1=1
  p2=1
  r2=0.1
  window=500
  
  # Arguments for PRS_serial script
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  
  # Significance thresholds specifications
  sig_thresholds=(5e-08 1e-06 1e-04 0.01 0.05 0.1 0.2 0.5 1)
  sig_thresholds_lower_bounds=(0 0 0 0 0 0 0 0 0)
 
  #Arguments specific to PRS set analysis
  Extra_analyses=FALSE
  Name_of_extra_analysis=NULL
  randomise=FALSE
  permutations=10000
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates" 
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions=both #either ( "extended" "normal" "both" )
  whole_genome_genic=FALSE
  Gene_specific_PRS=FALSE


elif [ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
home_OS="/Users"
elif [ "$whereami" = 'johnhubert-ThinkPad-P50' ]; then
home_OS="/home"
fi

if [[ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' || "$whereami" = 'johnhubert-ThinkPad-P50' ]]; then

  #path_to_CLOZUK="/mnt/databank/CLOZUK/GWAS/BGE/RSupdate"
  #path_to_Biobank="/c8000xd3/big-wpcvm/UKBB/"	
  path_to_PGC_conversion="Summary_stat_manipulation"
  path_to_CLOZUK_conversion="Genotype_dataset_manipulation"
  path_to_stationary_data="${home_OS}/johnhubert/Dropbox/Stationary_data/" 
  path_to_covariate_file="${path_to_stationary_data}CLOZUK2.r7.select2PC.eigenvec.txt"
  path_to_chromosome_length_file="${path_to_stationary_data}UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="${path_to_stationary_data}CLOZUK.r7.GWAS_IDs.fam"
 
   
  # Datasets

  # Dataset naming is unnessarily confusing here, essentially:\
  # - training set name and validation set name are the names that THE USER wants \
  # to call each dataset, they DO NOT need to match the file names
  # - the beginning (before '_table') of training_set_usually_summary needs to match to training_set_name
  # - all other arguments in dataset arguments need to match file names with relevant add on to names, ONLY change these parts of the argument
  # - eventually, I will just have four arguments of dataset names and file names only (I'm sorry for this Katherine!)
 
  training_set_name="BIP_PGC2"
  validation_set_name="COGSv2016_IMPUTE2" 
  training_set_original_filename="daner_PGC_BIP32b_mds7a"
  validation_set_full_name_without_chromosome="COGSv2016_IMPUTE2_Missing_hwe"
  
  # DO NOT ALTER!!!
  validation_set_usually_genotype="${validation_set_full_name_without_chromosome}_chr${chromosome_number}"
  validation_set_usually_genotype_serial="${validation_set_full_name_without_chromosome}_chr"
  training_set_usually_summary="${training_set_name}_table${chromosome_number}"
  
  # Pathway datasets
  Pathway_filename="Selected_Pocklington_plus_GO_pathways_SCHIZ.txt"
  Gene_location_filename="NCBI37.3.gene.loc"
  
  # These arguments below change from being in quotes to not, essentially, DON'T use quotes because they are not recognised in some POSIX style bash scripts.\
  # However, some arguments are read by R and I haven't tested whether ALL arguments can be without quotes yet
  # Again I don't like it but was short on time!

  # Split_by_chromosome for genotype?
  split_by_chromosome_required=FALSE 
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
  p1=1
  p2=1
  r2=0.1	
  window=500
  
  # Arguments for PRS_serial script
  # DEFUNCT ARGUMENTS?#
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  
  # Significance thresholds specifications
  sig_thresholds=(5e-08 1e-06 1e-04 0.01 0.05 0.1 0.2 0.5 1)
  sig_thresholds_lower_bounds=(0 0 0 0 0 0 0 0 0)
  # Lower bounds as optional format to match up with plink lower bounds
  
  # arguments specific to PRS set analyses
  Extra_analyses=FALSE
  Name_of_extra_analysis=NULL
  randomise=FALSE
  permutations=10000
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions=both #either ( extended normal both )	
  whole_genome_genic=FALSE
  Gene_specific_PRS=FALSE
  external_harddrive="FALSE"
fi

# This is defunct code as far as I'm aware, but I need to check back #
if [ "$external_harddrive" = "TRUE" ]; then
  path_to_harddrive=/Volumes/HD-PCU2
  cp $path_to_harddrive/CLOZUK_data/${validation_set_usually_genotype}.tar.gz .
  cp $path_to_harddrive/PGC_noCLOZUK_data/${training_set_usually_summary}.txt .
fi
