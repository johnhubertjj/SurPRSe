#! /bin/bash
# This script must be called with at least one parameter eg: sh Summary_stats_to_chromosome_converter.sh *.txt

./Taking_in_summary_stats_for_chromosome_conversion.R
./POLYGENIC_RISK_SCORE_ANALYSIS_CLOZUK_PRS_JOB_ARRAY.sh

 
