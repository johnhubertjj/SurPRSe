# Sources arguments to the main script, ensures directory is local, and creates extra folders where needed

# When deciding better structure for results in pathway analysis and combining whole genome with gene-set analysis, put the directory structure here

#Add arguments to the main body
source ${path_to_scripts}PRS_arguments_script.sh

if [ ${Using_raven} = "TRUE" ]; then
echo ${SLURM_SUBMIT_DIR}
cd $SLURM_SUBMIT_DIR
fi

Directory_to_work_from=`pwd`


cat ${path_to_scripts}PRS_arguments_script.sh


  if [[ ! -d "${training_set_name}_${validation_set_name}_output" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_output
  fi

  if [[ ! -d "${training_set_name}_${validation_set_name}_extrainfo" ]]; then
     mkdir ${training_set_name}_${validation_set_name}_extrainfo
  fi

