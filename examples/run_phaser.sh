#!/bin/bash

# all paths are relative prophaser/   (assumed this script is executed from phaser/examples)
cd ..

ne=11418
error=0.001
mapfile=examples/5_snps_interpolated_HapMap2_map_20
indir=examples/
sample_file=test_sample.vcf.gz
prefix=test_sample
ref_file=examples/test_ref.vcf.gz
results_directory=examples/out/


mkdir -p $results_directory

# Creates a template VCF that also contains a field for GP (posterir genotype probability).
./create_template_vcf_gtgp.sh $indir $sample_file

##### assuming the ompmod branch has been compiled to this executable
./phase_ompmod_cpu --outfile ${results_directory}${prefix}  --directory $indir --sample_file $sample_file --reference_file $ref_file --ne $ne --map $mapfile --error $error
