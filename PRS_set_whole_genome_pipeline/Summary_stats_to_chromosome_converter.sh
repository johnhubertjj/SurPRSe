
Rscript ${path_to_scripts}RscriptEcho.R ${path_to_scripts}Taking_in_summary_stat_data_for_chromosome_conversion.R ./${training_set_name}_${validation_set_name}_extrainfo/${training_set_name}_summary_stats_for_chromosome_conversion.Rout\
 ${path_to_training_dataset}\
 ${training_set_name}\
 ${validation_set_name}\
 ${MAF_summary}\
 ${INFO_summary}\
 ${SE_summary}\
 ${SE_threshold}\
 ${Chromosomes_to_analyse[@]}  

source ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt
cat ./${training_set_name}_${validation_set_name}_extrainfo/new_PRS_set_arguments_for_${training_set_name}.txt 

if [ "${split_by_chromosome_required}" = "TRUE" ]; then
   for i in "${Chromosomes_to_analyse[@]}" ;
   do
	plink --bfile ${path_to_validation_dataset} --chr ${i} --make-bed --out ${validation_set_usually_genotype_serial}${i}
	tar -zcvf ${validation_set_usually_genotype_serial}${i}.tar.gz ${validation_set_usually_genotype_serial}${i}.{bim,bed,fam,log}
	rm ${validation_set_usually_genotype_serial}${i}.{bed,bim,fam,log} 
   done
else
   
   for i in "${Chromosomes_to_analyse[@]}" ;
   do
	cp -v "${path_to_validation_dataset}_chr${i}.tar.gz" ${SLURM_SUBMIT_DIR}
   done
fi 


