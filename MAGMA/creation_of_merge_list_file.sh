#! /bin/bash

cd ~/Documents/testing_PRS_chromosome_22/test_chr5/output

chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

pathways=("FMRP_targets" "abnormal_behavior" "abnormal_nervous_system_electrophysiology" "abnormal_learning|memory|conditioning"  "abnormal_CNS_synaptic_transmission" "Cav2_channels" "abnormal_synaptic_transmission" "5HT_2C" "abnormal_long_term_potentiation" "abnormal_motor_capabilities|coordination|movement" "abnormal_behavioral_response_to_xenobiotic" "abnormal_associative_learning" "Lek2015_LoFintolerant_90" "BGS_top2_mean" "BGS_top2_max")

for a in `seq 0 14`;
do

if [ ! -d ${pathways[a]} ]; then
mkdir ${pathways[a]}
fi

if [ -f ./${pathways[a]}/make_full_plink_${pathways[a]}.txt ]; then
rm ./${pathways[a]}/make_full_plink_${pathways[a]}.txt
fi

shopt -s nullglob
number_of_files=(*${pathways[a]}SNPs_for_clumping.txt)
shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later
echo "${number_of_files[@]}"  # Note double-quotes to avoid extra parsing of funny characters in filenames
chromosome_number_total=${number_of_files}
length_of_array=`echo "${#number_of_files[@]}"`
num1=1
length_of_array=`echo "$((${length_of_array} - ${num1}))"`

        for i in `seq 0 ${length_of_array}` ;
        do
        chromosome_number_total=(`ls ${number_of_files[i]} | egrep -o [0-9]+`)
        number_of_files[i]=${chromosome_number_total[0]}
        echo ${number_of_files[i]}
	#printf "chromosome_${number_of_files[i]}_CLOZUK_PGC_${pathways[a]} %.0s" {1..2} >> ./${pathways[a]}/make_full_plink_${pathways[a]}.txt
	printf "chromosome_${number_of_files[i]}_CLOZUK_PGC_${pathways[a]}\n%.0s" {1} >> ./${pathways[a]}/make_full_plink_${pathways[a]}.txt  
        done
done

