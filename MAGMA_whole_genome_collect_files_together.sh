# assign arguments here for now because there are so many
  training_set_usually_summary="PGC_table"
  validation_set_usually_genotype="CLOZUK_GWAS_BGE_chr"
  training_set_name="PGC"
  validation_set_name="CLOZUK"
  MAF_summary="FALSE"
  MAF_threshold=0.01
  MAF_genotype="TRUE"
  INFO_summary="TRUE"
  INFO_threshold=0.9
  Chromosomes_to_analyse=(`seq 1 22`)
  Multiple_Training_set_tables="TRUE"
  Running_in_Serial="TRUE"
  sig_thresholds=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
  Perform_Magma_as_well="True"
  Magma_validation_set_name="_magma_input"
  # either "extended" "normal" or "both" : change to a numerical input in the future
  Gene_regions="both"

# Load both Plink and R
#  module purge
#  module load R/3.3.0
#  module load plink/1.9c3
#  module load python/2.7.11
#  module load magma/1.06

  cd ~/Dropbox/whole_genome_testing/
  path_to_scripts='/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/'
  number_of_files=($(find -E . -type f -regex "^./output/CLOZUK_GWAS_BGE_CLUMPED_chr[0-9]+.bed" -exec basename {} \;))
  path_to_PGC_conversion="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Summary_stat_manipulation"
  path_to_CLOZUK_conversion="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Genotype_dataset_manipulation"





# If running MAGMA as well
# Create output directory for MAGMA results
if [ ! -d "./output/MAGMA_set_analysis" ]; then
        mkdir ./output/MAGMA_set_analysis
fi

num1=1
        
# Calculate the number of files there are in the table and print them to the screen
number_of_files_magma=($(find -E . -type f -regex "^./output/${validation_set_name}_GWAS_BGE_chr[0-9]+${Magma_validation_set_name}.bed" -exec basename {} \;))       
length_of_array_magma=`echo "${#number_of_files_magma[@]}"`
echo "${number_of_files_magma[@]}"
length_of_array_magma=`echo "$((${length_of_array_magma} - ${num1}))"`

# Create a make-list file specifically for magma results 
if [ -f ./output/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt ]; then
       rm ./output/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
fi
        
for chromosome_number in `seq 0 ${length_of_array_magma}` ;
    do
        current_file=$(basename ${number_of_files_magma[$chromosome_number]} .bed)
        printf "${current_file}\n%.0s" {1} >> ./output/${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt
    done
        
cd ./output/
plink --merge-list ./${validation_set_name}_GWAS_BGE_${Magma_validation_set_name}_FULL_GENOME_MAKE_LIST_INPUT.txt --make-bed --out ./${validation_set_name}_${training_set_name}_MAGMA_FULL_GENOME
        
        

