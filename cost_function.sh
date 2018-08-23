#!/bin/bash

## redone to only take one ind and cov at a time
ne=$1
error=$2
ind=$3
cov=$4
niter=$5
split=$6

REFSET="CEU_"${ind}


declare -a COVS=("0.000000" "0.010000" "0.020000" "0.030000" "0.040000" "0.050000" "0.060000" "0.070000" "0.080000" "0.090000" "0.100000" "0.200000" "0.300000" "0.400000" "0.500000" "0.600000" "0.700000" "0.800000" "0.900000" "1.000000" "")
declare -a COVS=("0.500000")
#declare -a INDIVIDUALS=("NA12890" "NA12717" "NA11892" "NA12283" "NA12748")
declare -a INDIVIDUALS=("NA12890")



#for INDIVIDUAL in "${INDIVIDUALS[@]}"
#do

# for COV in "${COVS[@]}"
# do
# REFSET="CEU_"${INDIVIDUAL}
 #REFSET="YRI_ALL"
 #REFSET="CEU_10"
 #REFSET="CEU_"

# REFSET="CEU_50"

# ./phase --Ne $1 --individual $INDIVIDUAL --reference $REFSET --cov $COV

# done
#done

#echo ${INDIVIDUALS[*]}
#echo ${#INDIVIDUALS[@]}
#echo "now"
# Ne error refset numinds ind_1 .. ind_numinds numcovs cov1 .. cov_numcovs
#echo ../genoUtils/cost_function $ne $error $REFSET ${#INDIVIDUALS[@]} ${INDIVIDUALS[*]} ${#COVS[@]} ${COVS[*]}

#../genoUtils/cost_function $ne $error $REFSET ${#INDIVIDUALS[@]} ${INDIVIDUALS[*]} ${#COVS[@]} ${COVS[*]}


### redone to only take one ind and cov at a time

#./phase_lsmla --subset E --individual NA11892 --reference CEU_NA11892 --cov 0.000000 --Ne 25000 --seq_error 0.00100 --error 0.001 --split 5 --niter $niter
./phase_lsmla --subset E --individual $ind --reference $REFSET --cov $cov --Ne $ne --seq_error 0.00100 --error $2 --split $split --niter $niter

../genoUtils/cost_function $ne $error $REFSET 1 $ind 1 $cov $niter $split

