# phaser

phasing and imputation using Li & Stephens HMM algorithm 

Original reference for algorithm:
>Modeling Linkage Disequilibrium and Identifying Recombination Hotspots Using Single-Nucleotide Polymorphism Data
Na Li and Matthew Stephens

Another reference for the same algorithm. I've followed this one a lot, so e.g. transition probabilities and emission probabilities in the code correspond to the notation here:
>MaCH: Using Sequence and Genotype Data to Estimate Haplotypes and Unobserved Genotypes
Yun Li et.al.

Using normalized forward-backward algorithm for hidden markov models explained in e.g.
>https://pdfs.semanticscholar.org/54dc/c2a758e7fa34b8c2ef19826f39f16c4d1731.pdf


## build ##

build using makefile

two external libraries are used:

* https://github.com/statgen/libStatGen
* http://eigen.tuxfamily.org/index.php?title=Main_Page

these need to be downloaded and the path to them updated in makefile


## branches ##
two of the branches are "current": *master* and *linear_state*

### master ###
there are two versions, implemented in the classes HaplotypePhaser and HaplotypePhaserSym.

which one is used needs to be selected in main.cpp by commenting/uncommenting the pairs of lines:

```
//	HaplotypePhaser phaser;
//	string suffix = ".full";

	HaplotypePhaserSym phaser;
	string suffix = ".sym";
```

#### haplotypePhaser ####
The basic (=full) algorithm, with states being ordered pairs of reference haplotypes
(states represented by the class ChromosomePair)


#### haplotypePhaserSym ####
Removing symmetry in the state space by considering states as unordered pairs instead. 
(states represented by the class UnorderedChromosomePair, and the possible cases of transitioning from state to state (number of switches that occur) is just taken as the number of equal reference haplotypes between the states).

This rougly halves the size of the state space, and seems to give only slightly higher errors.

#### output ####

These both output two VCF files, one with phased haplotypes (.phased) and one with the most likely genotypes, calculated by summing together all the posterior state probabilities for each genotype (.genos). Genos gives much better imputation performance.


**If the directory specified for output does not exist, the run will fail. Thats why my example bash script creates the results directory.**


### linear_state ###
Branch where forward-backward is performed over one chromosome at a time. Much faster, but less correct. 

There are many different options for this code, mainly to do with how to handle the "other chromosome" while calculating probabilities for the current chromosome:

There are 3 ways implemented now:

These are specified using the command line argument --algorithm

#### fixed ####
 * Transition probabilities calculated considering both chromosomes (usual Li & Stephens algo) but with the other chromosome fixed at the state that was most ikely in the previous iteration.

* Emission probabilities also calculated considering both chromosomes, but with the other chromosomes value fixed at the _state_ with the highest posterior probability.

#### separate ####

* Transition probabilities calculated using pdf that only considers probability of switch happening along one chromosome.

* Emission probabilities calculated considering both chromosomes, but with the other chromosomes value fixed at the _allele_ with the highest posterior probability.

#### integrated ####
 * Transition probabilities calculated considering both chromosomes (usual Li & Stephens algo) and integrating over the possible states at the other chromosome, multiplying with their posterior probabilities from the previous iteration (slow)
  
 * Emission probabilities also calculated considering both chromosomes, but with the other chromosomes value fixed at the state with the highest posterior probability.

I'm not sure these are correct always in terms of the definition of pdfs and so on.
 
#### other command-line options ####

1. *--Geno_mls* which affects how the most likely genotypes are calulated. Each time after forward-backward has been done over a chromosome: posterior state and genotype probs are calculated. Thhis is done using the posterior state probabilites for the current chromosome and data for the "other" chromosome from the previous iteration. Turning this option on means the _state_ with highest posterior probability is used for the "other" chromosome, otherwise the most likely _allele_ is used.


2. *--iterations x* This is the number of times to alternate doing forward-backward over the chromosomes. Should be at least 2, but from what I've seen errors don't get better after 3 usually.


#### output ####

This branch outputs three VCF files:
1. .phased.MLS : phased using the state with highest posterior prob for each chromosome, for each marker
2. .phased.MLA : phased using the allele with highest posterior prob for each chromosome, for each marker
3. .genos : most likely genotypes, not phased.


## run ##

Example scripts in examples/

```
./run_phaser.sh
```


#### options ####
The parameters (also command-line options) *error*, *ne* (effective population size) and *map* affect performance a lot. The values in my script *run_phaser.sh* seem to work well on my test data. 

#### genetic map ####
Right now the code assumes all markers present in the reference file (= the markers that will be considered during imputation) are present in the map file. I created the example map file *examples/5_snps_interpolated_HapMap2_map_20* using linear interpolation from HapMap map files. If no mapfile is specified, the distance is set to 0.01 for all sites, which gives bad performance. Maybe lowering it will make it work ok with no map file.


#### template VCF ####
I couldnt get the libStatGen utilities to work for creating new VCFs, so now a temporary solution is to use a template vcf. The c++ code assumes a file with a specific name exists for every input filename. This is done by the script create_template_vcf.sh, also called by my example script. **Everything fails if a proper template does not exist, dont do big runs before ensuring the template thing works.**


