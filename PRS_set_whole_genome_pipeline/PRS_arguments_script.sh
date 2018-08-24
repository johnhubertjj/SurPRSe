# Arguments to be read into each analysis
# Run locally or on ARCCA

whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" = *"raven"* ]]; then

  # Assign a variable to show that you are on raven
  Using_raven=TRUE

  # Load modules for analysis
  module purge
  module load R/3.3.0
  module load plink/1.9c3
  module load python/2.7.11-genomics
  module load magma/1.06
  module load parallel/20170322

  # Path_to_locations_on_the_server
  path_to_PRS_scripts="${HOME}/Schizophrenia_pipeline_scripts/"

  # Re-assign to the training_set_usually_genotype and validation_full_name_without_chromosome using the sed command
  path_to_validation_dataset="/scratch/${USER}/PRS_tutorial/eur-1000g-phase3-chrall-mac5"
  path_to_training_dataset="/scratch/$USER/PRS_tutorial/iPSYCH-PGC_ASD_Nov2017.tsv.gz"

  # State paths to the relevant stationary folders required for the analysis
  path_to_PGC_conversion="${path_to_PRS_scripts}/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="${path_to_PRS_scripts}/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="${path_to_PRS_scripts}/MAGMA/"

  # Path to where stationary file is kept
  path_to_stationary_data="${HOME}/Stationary_data/"
  path_to_covariate_file="${path_to_stationary_data}COGSv2016.noREL.txt"
  path_to_chromosome_length_file="${path_to_stationary_data}UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="${path_to_stationary_data}CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="${path_to_stationary_data}NCBI37.3.gene.loc"

  # Datasets
  training_set_name="iPsych_GWAS"
  validation_set_name="thousand_genomes"

  # Sample sizes
  number_of_samples_training=78308

  # DO NOT ALTER!!!
  training_set_original_filename=`echo "${path_to_training_dataset}" | sed 's:.*/::'`
  validation_set_full_name_without_chromosome=`echo "${path_to_validation_dataset}" | sed 's:.*/::'`
  validation_set_usually_genotype_serial="${validation_set_full_name_without_chromosome}_chr"
  training_set_usually_summary="${training_set_name}_table"

  # Path to LDscore regression python script
  ldsc="/home/c1020109/ldsc/ldsc.py"

  # Pathway datasets
  Pathway_filename="Pocklington2015_134sets_LoFi.txt"
  Gene_location_filename="NCBI37.3.gene.loc"
  calculate_indep_SNPs=FALSE

  # Split_by_chromosome for genotype?
  split_by_chromosome_required=TRUE

  # Do a missingness check?
  Missing_geno=TRUE
  genotype_missingness_check=0.1
  HWE_thresh=1e-06
  hwe_p_test="midp"
  include_noncntrl="include-nonctrl"

  # MAF, INFO and SE
  MAF_summary=FALSE
  MAF_threshold=0.1
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

  # Arguments specific to PRS set analysis
  Extra_analyses=FALSE
  Full_genome_PRS_extra_analysis=FALSE
  Name_of_extra_analysis=(Pathways Genes)
  randomise=FALSE
  sample_replace=FALSE
  permutations=1000
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"

  # Pathway specific arguments
  Use_files_for_pattern_matching_make_list=TRUE

  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions=both #either ( "extended" "normal" "both" )
  whole_genome_genic=FALSE
  Gene_specific_PRS=FALSE

  # MAGMA analysis
  MAGMA_gene_set_analysis=TRUE

elif [ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then
home_OS="/Users"
elif [ "$whereami" = 'johnhubert-ThinkPad-P50' ]; then
home_OS="/home"
fi

if [[ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' || "$whereami" = 'johnhubert-ThinkPad-P50' ]]; then

  # Assign a variable to show that you are on raven
  Using_raven=FALSE

  # Path_to_locations_on_the_server
  path_to_PRS_scripts="${HOME}/PhD_scripts/Schizophrenia_pipeline_scripts/"

  # Re-assign to the training_set_usually_genotype and validation_full_name_without_chromosome using the sed command
  path_to_validation_dataset="/scratch/${USER}/PR54/PGC_CLOZUK_PRS/COGS_CLOZUK_datasets/COGSv2016_IMPUTE2"
  path_to_training_dataset="/scratch/${USER}/PR54/Summary_stats_files/IQ_2017_GWAS.txt"

  # State paths to the relevant stationary folders required for the analysis
  path_to_PGC_conversion="${path_to_PRS_scripts}/Summary_stat_manipulation/"
  path_to_CLOZUK_conversion="${path_to_PRS_scripts}/Genotype_dataset_manipulation/"
  path_to_MAGMA_scripts="${path_to_PRS_scripts}/MAGMA/"

  # Path to where stationary file is kept
  path_to_stationary_data="${HOME}/Stationary_data/"
  path_to_covariate_file="${path_to_stationary_data}COGSv2016.noREL.txt"
  path_to_chromosome_length_file="${path_to_stationary_data}UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
  path_to_new_fam_file="${path_to_stationary_data}CLOZUK.r7.GWAS_IDs.fam"
  path_to_gene_annotation_file="${path_to_stationary_data}NCBI37.3.gene.loc"

  # Datasets
  training_set_name="IQ"
  validation_set_name="COGS"

  # Sample sizes
  number_of_samples_training=78000

  # DO NOT ALTER!!!
  training_set_original_filename=`echo "${path_to_training_dataset}" | sed 's:.*/::'`
  validation_set_full_name_without_chromosome=`echo "${path_to_validation_dataset}" | sed 's:.*/::'`
  validation_set_usually_genotype_serial="${validation_set_full_name_without_chromosome}_chr"
  training_set_usually_summary="${training_set_name}_table"

  # Path to LDscore regression python script
  ldsc="/home/c1020109/ldsc/ldsc.py"

  # Pathway datasets
  Pathway_filename="Pocklington2015_134sets.txt"
  Gene_location_filename="NCBI37.3.gene.loc"
  calculate_indep_SNPs=FALSE

  # Split_by_chromosome for genotype?
  split_by_chromosome_required=TRUE

  # Do a missingness check?
  Missing_geno=TRUE
  genotype_missingness_check=0.1
  HWE_thresh=1e-06
  hwe_p_test="midp"
  include_noncntrl="include-nonctrl"

  # MAF, INFO and SE
  MAF_summary=FALSE
  MAF_threshold=0.1
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

  # Arguments specific to PRS set analysis
  Extra_analyses=TRUE
  Full_genome_PRS_extra_analysis=TRUE
  Name_of_extra_analysis=(Pathways Genes)
  randomise=FALSE
  sample_replace=FALSE
  permutations=1000
  Magma_validation_set_name="_consensus_with_${training_set_name}_flipped_alleles_no_duplicates"

  # Pathway specific arguments
  Use_files_for_pattern_matching_make_list=TRUE

  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions=both #either ( "extended" "normal" "both" )
  whole_genome_genic=TRUE
  Gene_specific_PRS=FALSE

  # MAGMA analysis
  MAGMA_gene_set_analysis=FALSE
fi
