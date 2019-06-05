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

https://github.com/statgen/libStatGen
http://eigen.tuxfamily.org/index.php?title=Main_Page

these need to be downloaded and the path to them updated in makefile


## branches ##
two of the branches are "current": master and linear_state

### master ###
there are two versions, implemented in the classes HaplotypePhaser and HaplotypePhaserSym
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


These both output two VCF files, one with phased haplotypes (.phased) and one with the most likely genotypes, calculated by summing together all the posterior state probabilities for each genotype (.genos). Genos gives much better imputation performance.


**If the directory specified for output does not exist, the run will fail. Thats why my example bash script creates the results directory.**


### linear_state ###
Branch where forward-backward is performed over one chromosome at a time. Much faster, but less correct. 

There are many different options for this code, mainly to do with how to handle the "other chromosome" while calculating probabilities for the current chromosome. There are 3 ways implemented now:

These are specified using the command line argument --algorithm

#### fixed ####

#### separate ####

#### integrated ####

There is also the command line option *--Geno_mls* which affects how the most likely genotypes are calulated: turning this option on means  


I'm not sure these are correct always in terms of the definition of pdfs and so on.

Can also specify number of iterations using command line argument --iterations x.
This is the number of times to alternate doing forward-backward over the chromosomes. Should be at least 2, but from what I've seen errors don't get better after 3 usually.

There is also an option 


#### run ####
error
ne
data/map
examples

