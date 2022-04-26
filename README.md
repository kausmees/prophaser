# phaser

Phasing and imputation using the Li & Stephens approximate coalescent algorithm for haplotype frequencies. 

Original reference for algorithm:
> [1] Modeling Linkage Disequilibrium and Identifying Recombination Hotspots Using Single-Nucleotide Polymorphism Data. N. Li and M. Stephens. 2013.

Another reference where the same underlying model is implemented. Some of the optimizations explained in this paper have also been implemented in this code:
> [2] MaCH: Using Sequence and Genotype Data to Estimate Haplotypes and Unobserved Genotypes. Y. Li, C. J. Willer, J. Ding, et al. 2010.


This software implementes a normalized forward-backward algorithm for hidden markov models explained in e.g.
>https://pdfs.semanticscholar.org/54dc/c2a758e7fa34b8c2ef19826f39f16c4d1731.pdf


## build ##

build using makefile

two external libraries are used:

* https://github.com/statgen/libStatGen
* http://eigen.tuxfamily.org/index.php?title=Main_Page

these need to be downloaded and the path to them updated in makefile


## branches ##
Two of the branches are current: *master* and *ompmod*:

### master ###
This branch implements the optimization based on collapsing of transition probabilities from [2], resulting in an effective computational complexity quadratic in the number of haplotypes.

This outputs two VCF files, one with phased haplotypes (.phased) and one with the most likely genotypes, calculated by summing together all the posterior state probabilities for each genotype (.genos). The genos tends to give much better imputation performance.

### ompmod ###
This branch implements an additional optimization for memory usage, using a two-level version of the checkpointing scheme from [2] to obtain a more cache-friendly algorithm that can be executed on GPUs. Uses OpenMP offloading pragmas and can be compiled for multicore CPUs and GPUs.

This outputs one VCF file, with the most likely genotypes and their posterior probabilities (.postgenos)

This branch corresponds to the implementation described in the paper "Achieving improved accuracy for imputation of ancient DNA".


# Running

Example scripts in [/examples](/examples)

```
./run_phaser.sh
```

# Notes

## genetic maps
Currently, the code assumes that all markers present in the reference file (= the markers that will be considered during imputation) are present in the map file. The directory [/examples](/examples) contains a script that performs linear interpolation of a genetic map to extend it to a specific set of genomic positions. If no mapfile is specified, the distance is set to 0.01 for all sites, which tends to give not so good performance.


## template VCF
I couldnt get the libStatGen utilities to work for creating new VCFs, so now a temporary solution is to use a template vcf. The imputation code assumes a template VCF file with a specific name exists for every input filename. Such a file is created by the script create_template_vcf.sh (or create_template_vcf_gtgp.sh for branch ompmod), which is called in the example script [/examples/run_phaser.sh](/examples/run_phaser.sh).

**The run will fail if:**
 * If the directory specified for output does not exist.
 * If a proper template does not exist.
 
**Don't do big runs before ensuring the above.**




