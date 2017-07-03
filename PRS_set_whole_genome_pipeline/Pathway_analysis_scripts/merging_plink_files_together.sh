#! /bin/bash

cd ~/Documents/testing_PRS_chromosome_22/test_chr5/output

chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

pathways=("FMRP_targets" "abnormal_behavior" "abnormal_nervous_system_electrophysiology" "abnormal_learning|memory|conditioning"  "abnormal_CNS_synaptic_transmission" "Cav2_channels" "abnormal_synaptic_transmission" "5HT_2C" "abnormal_long_term_potentiation" "abnormal_motor_capabilities|coordination|movement" "abnormal_behavioral_response_to_xenobiotic" "abnormal_associative_learning" "Lek2015_LoFintolerant_90" "BGS_top2_mean" "BGS_top2_max")

for a in `seq 0 14`;
do
plink --merge-list ./${pathways[a]}/make_full_plink_${pathways[a]}.txt --make-bed --out ./${pathways[a]}/${pathways[a]}_Clumping_input
done

