#!/bin/bash

path_to_scripts="/Users/johnhubert/Documents/PhD_scripts/Schizophrenia_pipeline_scripts/PRS_set_whole_genome_pipeline/"
Directory_to_work_from="/Users/johnhubert/CLOZUK_ALSPAC_DATA_PRS"
chromosomes=(`seq 1 22`)

sh Summary_stats_to_chromosome_converter ~/Documents/CLOZUK_ALSPAC_DATA_PRS

parallel ./POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh ::: ${chromosomes[@]}

# Use awk to extract the argument Extra_analysis from your output scripts in question #output the arguments to a temporary file in order to read it

if [[ Extra_analysis == TRUE ]]; then
