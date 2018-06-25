# use the require function to make sure that the packages are up to date with the software

# Create arguments files for easier input into Rscripts for parallelisation # 
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Chromosome_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_Chromosome_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${Chromosomes_to_analyse[@]}

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}sig_thresholds_lower_bounds_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_sig_thresholds_lower_bounds_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${sig_thresholds_lower_bounds[@]}

Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}sig_thresholds_plink_arguments_text_file.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_sig_thresholds_plink_arguments_text_file.Rout\
 ${training_set_name}\
 ${validation_set_name}\
 ${sig_thresholds[@]}

