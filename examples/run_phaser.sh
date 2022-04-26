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


mkdir -p $results_directory

./create_template_vcf.sh $indir $sample_file

##### assuming the master branch has been compiled to this executable
./phase_master --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file --ne $ne --map $mapfile --error $error 




