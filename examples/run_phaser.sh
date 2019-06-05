#!/bin/bash


# all paths are relative phaser/   (assumed this script is executed from phaser/examples)
cd ..

ne=11418
error=0.001
mapfile=examples/5_snps_interpolated_HapMap2_map_20
indir=examples/
sample_file=test_sample.vcf.gz
ref_file=examples/test_ref.vcf.gz
results_directory=examples/out/

# parameters for the linear_state  branch
algo=integrated
algo=separate
algo=fixed
niter=3



mkdir -p $results_directory

./create_template_vcf.sh $indir $sample_file

##### assuming the master branch has been compiled to this executable
./phase --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file --ne $ne --map $mapfile --error $error 

##### assuming the linear-state branch has been compiled to this executable
#./phase_linear_state --map $mapfile  --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file  --iterations  $niter --ne $ne --error $error --algorithm $algo





