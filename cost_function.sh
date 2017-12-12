#!/bin/bash


declare -a COVS=("0.000000" "0.010000" "0.020000" "0.030000" "0.040000" "0.050000" "0.060000" "0.070000" "0.080000" "0.090000" "0.100000" "0.200000" "0.300000" "0.400000" "0.500000" "0.600000" "0.700000" "0.800000" "0.900000" "1.000000" "")

declare -a COVS=("")
#declare -a COVS=("0.020000")


#declare -a INDIVIDUALS=("NA12890" "NA12717" "NA11892" "NA12283" "NA12748")
declare -a INDIVIDUALS=("NA12812")

for INDIVIDUAL in "${INDIVIDUALS[@]}"
do

 for COV in "${COVS[@]}"
 do
 #REFSET="CEU_ALL_"${INDIVIDUAL}
 #REFSET="YRI_ALL"
 #REFSET="CEU_10"
 REFSET="CEU_25"

# REFSET="CEU_50"

 ./phase --Ne $1 --individual $INDIVIDUAL --reference $REFSET --cov $COV

 done
done

echo ${INDIVIDUALS[*]}
echo ${#INDIVIDUALS[@]}

# Ne error refset numinds ind_1 .. ind_numinds numcovs cov1 .. cov_numcovs
../genoUtils/geno $1 $2 $REFSET ${#INDIVIDUALS[@]} ${INDIVIDUALS[*]} ${#COVS[@]} ${COVS[*]}

