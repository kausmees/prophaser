#!/bin/bash


# all paths are relative phaser/   (assumed this script is executed from phaser/examples)

ne=11418
#error=0.001
error=0.00001
mapfile=examples/5_snps_interpolated_HapMap2_map_20
indir=Camille/
sample_file=IMP.chr20.pooled.snps.gl.chunk10000.vcf.gz
ref_file=$indir/REF.chr20.snps.gt.chunk10000.vcf.gz
results_directory=Camille_out_190822/
sample=$SLURM_ARRAY_TASK_ID

# parameters for the linear_state  branch
algo=integrated
algo=separate
algo=fixed
niter=3



mkdir -p $results_directory

#./create_template_vcf.sh $indir $sample_file

#export OMP_NUM_THREADS=8

##### assuming the master branch has been compiled to this executable
#perf stat -d -d -d -I 60000
./phase --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file --ne $ne --map $mapfile --error $error --sampleId $sample

##### assuming the linear-state branch has been compiled to this executable
#./phase_linear_state --map $mapfile  --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file  --iterations  $niter --ne $ne --error $error --algorithm $algo





