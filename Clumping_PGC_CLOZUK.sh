#!/bin/bash

#PBS -q batch-long
#PBS -P PR54
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime 10:00:00
#PBS -o PRS_CLOZUK_combine.txt
#PBS -e PRS_CLOZUK_combine_error.txt
#PBS -j oe
#PBS -J 1-22
#PBS -r y

#PGC_CLOZUK_SNP_table $PBS_ARRAY_INDEX .txt.gz


