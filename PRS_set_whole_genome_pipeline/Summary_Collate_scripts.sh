
# Run locally or on ARCCA
whereami=$(uname -n)
echo "$whereami"
system=$4


if [[ "$system" = "MAC" || "$system" = "LINUX" ]]; then
  # add directory to work from at the top of the script (just in case)
  Directory_to_work_from=$2
  cd ${Directory_to_work_from}

  #add the chromosome this copy of the script will analyse
  chromosome_number=$1
  echo $chromosome_number

  # path to where the scripts are located
  path_to_scripts=$3

  # Assign the shell variables
  source ${path_to_scripts}/PRS_arguments_script.sh
  cat ${path_to_scripts}/PRS_arguments_script.sh

  # Alter/add variables depending on what type of Training dataset you have
  source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
  echo ${validation_set_usually_genotype_serial}${chromosome_number}
fi

Summary_analysis=TRUE

# Run R script that will combine PGC and CLOZUK to an individual table
# Output is in PGC_CLOZUK_SNP_table.txt
Rscript ${path_to_scripts}RscriptEcho.R\
 ${path_to_scripts}CLOZUK_PGC_COMBINE_final.R\
 ./${training_set_name}_${validation_set_name}_extrainfo/${validation_set_name}_${training_set_name}_COMBINE_chr${chromosome_number}.Rout\
 ${training_set_usually_summary}${chromosome_number}\
 ${training_set_name}_${validation_set_name}_output/${validation_set_usually_genotype_serial}${chromosome_number}_2\
 ${training_set_name}\
 ${validation_set_name}\
 ${chromosome_number}\
 ${CHR_name}\
 ${SNP_name}\
 ${BP_name}\
 ${A1_name}\
 ${A2_name}\
 ${OR_name}\
 ${BETA_name}\
 ${system}\
 ${Number_of_frequency_columns}\
 ${path_to_scripts}\
 ${Summary_analysis}

