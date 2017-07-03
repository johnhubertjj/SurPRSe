#!/bin/bash

path_to_scripts = /home/$USER/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/PRS_set_whole_genome_pipeline/"
$1=location_of_dataset_names_to_change_in_arguments_script
$2=type_of_analysis


if [[ ${type_of_analysis} == "Whole_genome" ]] then;
# Need to figure out how to list the different names to loop through various job scripts
PREAMBLE = ${qsub ...}

for i in whatever : 

do 

FIRST=$(qsub ${path_to_scripts}POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh)
echo $FIRST



testing_jobs=$(qstat -u c1020109)
# test on a running job, match the name and obtain the jobID for each array job, or maybe just check if the variable is blank

a=0

while $a == 0: 
do
	echo "Waiting for job 1 to finish"
	sleep 30
		if [[ testing_jobs == #empty ]] 
			a=1     
		fi
done

# here put a kill switch, aka check if a filename still exists across a range of variables (eg 1:22 of a file you have removed at the end of the script)

SECOND=$(qsub -W depend=afterany:$FIRST job2.pbs)
echo $SECOND

THIRD=$(qsub -W depend=afterany:$SECOND job3.pbs)
echo $THIRD
done


if [[ ${type_of_analysis} == "Pathways" ]] then;

# different type of pipeline


