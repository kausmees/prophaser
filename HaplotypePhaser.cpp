#include <math.h>
#include <string.h>
#include <omp.h>

#include "HaplotypePhaser.h"
#include "MemoryAllocators.h"


HaplotypePhaser::~HaplotypePhaser(){
	delete [] normalizers;

	FreeDoubleMatrix(s_forward, num_markers);
	FreeDoubleMatrix(s_backward, num_markers);

}


void HaplotypePhaser::AllocateMemory(){
	String::caseSensitive = false;

	num_markers = Pedigree::markerCount;

	pop_const = (4.0 * Ne) / 100.0;

	// Number of reference inds.
	num_ref_inds = ped.count;
	num_inds = num_ref_inds +1;

	printf("Num inds tot: %d \n", num_inds);


	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	num_haps = num_ref_inds*2;

	num_states = num_haps;

	L = MatrixXd(6,5);
	u = VectorXd(6);
	u << 0, 0, 0 , 0 , 0, 1;
	x = VectorXd(5);
	coefs_p0 = VectorXd(5);
	coefs_p1 = VectorXd(5);
	sample_gls.resize(num_markers*3);


	// default value for genetic map. using default does not give very good results.
	distances.resize(num_markers,0.01);

	all_posteriors = MatrixXd(num_markers,num_states);

	// haplotypes(h,m) = 0/1 representing allele of haplotype h at marker m
	haplotypes = MatrixXc(num_haps+2, num_markers);


	sample_gls.resize(num_markers*3);

	s_forward = AllocateDoubleMatrix(num_markers, num_states);
	s_backward = AllocateDoubleMatrix(num_markers, num_states);

	normalizers = new double[num_markers];

	ml_states_c1 = new int[num_markers];
	ml_states_c2 = new int[num_markers];

	ml_alleles_c1 = new int[num_markers];
	ml_alleles_c2 = new int[num_markers];

};

/**
 * Load reference data from specified file.
 *
 * Load genetic distances from map file if specified.
 *
 */
void HaplotypePhaser::LoadReferenceData(const String &ref_file, String &map_file){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadReferenceIndividuals(ped,ref_file);
	AllocateMemory();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);
	if(!map_file.IsEmpty()) {
		VcfUtils::LoadGeneticMap(map_file.c_str(), ped, distances);
	};
};

/**
 * Load sample data from vcf.
 *
 * Overwrites sample_gls, assumed handling one sample at a time.
 *
 */
void HaplotypePhaser::LoadSampleData(const String &sample_file, int sample_index){
	VcfUtils::LoadSampleIndividual(ped, sample_file, sample_index);
	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);

};

/**
 * Fill allele_counts.
 *
 * allele_counts(m) = number of haplotypes with allele 1 at marker m
 *
 */
void HaplotypePhaser::CalcAlleleCounts() {
	allele_counts = haplotypes.colwise().sum();
};




/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * for all hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 * for fixed reference haplotype at chromosome 1
 */

void HaplotypePhaser::CalcEmissionProbs(int marker, double * probs) {
	int h1;
	int h2;
	double sum;

	double case_1 = (pow(1 - error, 2) + pow(error, 2));
	double case_2 = 2 * (1 - error) * error;
	double case_3 = pow(1 - error, 2);
	double case_4 = (1 - error) * error;
	double case_5 = pow(error,2);

	int * ml_states_other;

	if (curr_chrom == 2) {
		ml_states_other = ml_states_c1;
	}
	if (curr_chrom == 1) {
		ml_states_other = ml_states_c2;
	}

	for (int state = 0; state < num_states; state++) {

		// Reference hapotype at chromosome 1 - fixed (0: REF 1: ALT)
		h1 = haplotypes(ml_states_other[marker], marker);

		// Reference hapotype at chromosome 2 (0: REF 1: ALT)
		h2 = haplotypes(state,marker);

		sum = 0.0;


		// case1: g = 0
		if(h1 == 0 and h2 == 0){
			sum += case_3 *  sample_gls[marker * 3];
		}
		else {
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){

				sum += case_4 * sample_gls[marker * 3];

			}
			else{
				sum +=  case_5 *  sample_gls[marker * 3];

			}
		}

		// case2: g = 1
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
			sum += case_1 *  sample_gls[marker * 3 + 1];

		}
		else{
			sum +=case_2 * sample_gls[marker * 3 + 1];

		}

		// case3: g = 2
		if(h1 == 1 and h2 == 1){
			sum += case_3 * sample_gls[marker * 3 + 2];

		} else{
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum += case_4 * sample_gls[marker * 3 + 2];

			}
			else{
				sum += case_5 * sample_gls[marker * 3 + 2];

			}
		}

		probs[state] = sum;
	}
}


void HaplotypePhaser::CalcEmissionProbsMla(int marker, double * probs) {
	int h1;
	int h2;
	double sum;

	double case_1 = (pow(1 - error, 2) + pow(error, 2));
	double case_2 = 2 * (1 - error) * error;
	double case_3 = pow(1 - error, 2);
	double case_4 = (1 - error) * error;
	double case_5 = pow(error,2);

	int * ml_alleles_other;

	if (curr_chrom == 2) {
		ml_alleles_other = ml_alleles_c1;

	}
	if (curr_chrom == 1) {
		ml_alleles_other = ml_alleles_c2;

	}


	// OPT1
	for (int state = 0; state < num_states; state++) {

		// Reference hapotype at chromosome 1 - fixed (0: REF 1: ALT)
		h1 = ml_alleles_other[marker];

		// Reference hapotype at chromosome 2 (0: REF 1: ALT)
		h2 = haplotypes(state,marker);

		sum = 0.0;
		// case1: g = 0

		if(h1 == 0 and h2 == 0){
			sum += case_3 *  sample_gls[marker * 3];
			//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
		}
		else {
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){

				sum += case_4 * sample_gls[marker * 3];
				//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

			}
			else{
				sum +=  case_5 *  sample_gls[marker * 3];
				//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

			}
		}

		// case2: g = 1
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
			sum += case_1 *  sample_gls[marker * 3 + 1];
			//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

		}
		else{
			sum +=case_2 * sample_gls[marker * 3 + 1];
			//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

		}

		// case3: g = 2
		if(h1 == 1 and h2 == 1){
			sum += case_3 * sample_gls[marker * 3 + 2];
			//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

		} else{
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum += case_4 * sample_gls[marker * 3 + 2];
				//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			}
			else{
				sum += case_5 * sample_gls[marker * 3 + 2];
				//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			}
		}

		probs[state] = sum;
	}
}



/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * forall hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 * for fixed reference haplotype at chromosome 1
 */

void HaplotypePhaser::CalcEmissionProbsIntegrated(int marker, double * probs) {
	int h1;
	int h2;
	double this_p;
	double case_1 = (pow(1 - error, 2) + pow(error, 2));
	double case_2 = 2 * (1 - error) * error;
	double case_3 = pow(1 - error, 2);
	double case_4 = (1 - error) * error;
	double case_5 = pow(error,2);

	double sum_geno0;
	double sum_geno1;
	double sum_geno2;



	int * ml_states_other;

	if (curr_chrom == 2) {
		ml_states_other = ml_states_c1;
	}
	if (curr_chrom == 1) {
		ml_states_other = ml_states_c2;
	}



	for (int state = 0; state < num_states; state++) {
		//		cout << "Doing emission probs state "<< state << endl;

		// First reference hapotype at state (0: REF 1: ALT)
		h1 = haplotypes(state, marker);

		sum_geno0 = 0.0;
		sum_geno1 = 0.0;
		sum_geno2 = 0.0;


		for (int state2 = 0; state2 < num_states; state2++) {

			// Second reference hapotype at state (0: REF 1: ALT)
			h2 = haplotypes(state2, marker);

			this_p = all_posteriors(marker,state2);

			// case1: g = 0

			if(h1 == 0 and h2 == 0){
				sum_geno0 += case_3 * this_p;
			}
			else {
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno0 += case_4 * this_p;

				}
				else{
					sum_geno0 += case_5 * this_p;

				}
			}

			// case2: g = 1
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum_geno1 += case_1 * this_p;

			}
			else{
				sum_geno1 += case_2 * this_p;

			}

			// case3: g = 2
			if(h1 == 1 and h2 == 1){
				sum_geno2 += case_3 * this_p;

			} else{
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno2 += case_4 * this_p;

				}
				else{
					sum_geno2 += case_5 * this_p;

				}
			}
		}

		if (abs(sum_geno0+sum_geno1+sum_geno2 -1 ) > 0.0001){
			cout << "Geno probs in emission: : "<< sum_geno0 << " " << sum_geno1 << " " << sum_geno2 << " ..... = "<< sum_geno0+sum_geno1+sum_geno2<<  endl;
		}

		probs[state] =
				sample_gls[marker * 3  ] * sum_geno0 +
				sample_gls[marker * 3+1] * sum_geno1 +
				sample_gls[marker * 3+2] * sum_geno2;
	}


}


/**
 *  Calculate marginal state probs: p(s) marginalized over emission probabilities for marker marker.
 *  p(s) = p0 if allele implied by s is 0
 *  p(s) = p1 if allele implied by s is 1
 */

void HaplotypePhaser::CalcMarginalStateProbs(int marker) {

	int q = allele_counts[marker];
	int p = num_haps - q;

	L << 0,pow((1-error),2),pow((1-error),2)-1,0,pow((1-error),2),
			0,2*(1-error)*error-1, 2*(1-error)*error, 0,2*(1-error)*error,
			0,pow(error,2),pow(error,2),0,(pow(error,2))-1,
			(1-error)*error,0,0,(1-error)*error*2-1,0,
			1-(pow(error,2))-pow(1-error,2),0,0,-2*(pow(1-error,2) + pow(error,2)),0,
			2*p*q,pow(p,2) +pow(q,2),pow(p,2) +pow(q,2),4*p*q,pow(p,2) +pow(q,2);


	x = L.bdcSvd(ComputeThinU | ComputeThinV).solve(u);
	//	cout << "cases primes = :\n" << x << endl;

	coefs_p0 << q, p, p , 2*q , p;
	coefs_p1 << p, q, q , 2*p , q;


	p0 = coefs_p0.dot(x);
	p1 = coefs_p1.dot(x);

	//	cout << "P0: " << p0 << endl;
	//	cout << "P1: " << p1 << endl;

}



/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * for haplotype 1 - marginalizing over haplotype 2
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 * where state space consists of only one haplotype (num_states = num_inds)
 */

void HaplotypePhaser::CalcEmissionProbsMarginalized(int marker, double * probs) {
	int h1;
	int h2;
	double sum_geno0;
	double sum_geno1;
	double sum_geno2;

	CalcMarginalStateProbs(marker);

	double this_p;

	double case_1 = (pow(1 - error, 2) + pow(error, 2));
	double case_2 = 2 * (1 - error) * error;
	double case_3 = pow(1 - error, 2);
	double case_4 = (1 - error) * error;
	double case_5 = pow(error,2);



	for (int state = 0; state < num_states; state++) {
		//		cout << "Doing emission probs state "<< state << endl;

		// First reference hapotype at state (0: REF 1: ALT)
		h1 = haplotypes(state, marker);

		sum_geno0 = 0.0;
		sum_geno1 = 0.0;
		sum_geno2 = 0.0;


		for (int state2 = 0; state2 < num_states; state2++) {

			// Second reference hapotype at state (0: REF 1: ALT)
			h2 = haplotypes(state2, marker);

			//			cout << "s1:"<< state  << "  s2:"<< state2 << endl;
			//			cout << "h1:"<< h1  << "  h2:"<< h2 << endl;

			// p(state2) depends on if the allele at this state is 0 or 1
			this_p = h2*p1 + (1-h2)*p0;

			// case1: g = 0

			if(h1 == 0 and h2 == 0){
				sum_geno0 += case_3 * this_p;
			}
			else {
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno0 += case_4 * this_p;

				}
				else{
					sum_geno0 += case_5 * this_p;

				}
			}

			// case2: g = 1
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum_geno1 += case_1 * this_p;

			}
			else{
				sum_geno1 += case_2 * this_p;

			}

			// case3: g = 2
			if(h1 == 1 and h2 == 1){
				sum_geno2 += case_3 * this_p;

			} else{
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno2 += case_4 * this_p;

				}
				else{
					sum_geno2 += case_5 * this_p;

				}
			}
		}
		//		cout << "Doing probs: "<< sum_geno0 << " " << sum_geno1 << " " << sum_geno2 << endl;
		probs[state] =
				sample_gls[marker * 3  ] * sum_geno0 +
				sample_gls[marker * 3+1] * sum_geno1 +
				sample_gls[marker * 3+2] * sum_geno2;
	}

}



void HaplotypePhaser::InitPriorScaledForward(){
	double * emission_probs = new double[num_states];
	double prior = 1.0 / num_states;
	double c1 = 0.0;

	CalcEmissionProbs(0, emission_probs);

	for(int s = 0; s < num_states; s++){
		c1 += emission_probs[s];
	};

	normalizers[0] = 1.0/(prior*c1);

	for(int s = 0; s < num_states; s++){
		s_forward[0][s] = (prior*emission_probs[s]) * normalizers[0];
	};

	delete [] emission_probs;
};


void HaplotypePhaser::InitPriorScaledForwardMarginalized(){
	double * emission_probs = new double[num_states];
	double prior = 1.0 / num_states;
	double c1 = 0.0;


	CalcEmissionProbsMarginalized(0, emission_probs);

	for(int s = 0; s < num_states; s++){
		c1 += emission_probs[s];
	};

	normalizers[0] = 1.0/(prior*c1);
	//	printf("First Normalizer = %e Prior = %e c1 = %e \n", normalizers[0], prior, c1);

	for(int s = 0; s < num_states; s++){
		s_forward[0][s] = (prior*emission_probs[s]) * normalizers[0];
	};

	delete [] emission_probs;
};




void HaplotypePhaser::InitPriorScaledBackward(){

	for(int s = 0; s < num_states; s++){
		s_backward[num_markers-1][s] = normalizers[num_markers-1];
	};
};


/*
 * Calculate scaled forward values for one haplotype (curr_chrom).
 *
 * Transition probabilities calculated considering both chromosomes (usual Li & Stephens algo) but with the other chromosome fixed at the state that was most
 * likely in the previous iteration.
 *
 * Emission probabilities also calculated considering both chromosomes, but with the other chromosomes value fixed at the state with the highest posterior probability.
 *
 *
 */
void HaplotypePhaser::CalcScaledForwardFixed(){
	double * emission_probs = new double[num_states];
	int * ml_states_other;
	double c;
	double probs[3];
	double scaled_dist;


	InitPriorScaledForward();


	if (curr_chrom == 2) {
		ml_states_other = ml_states_c1;
	}
	if (curr_chrom == 1) {
		ml_states_other = ml_states_c2;
	}


	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbs(m, emission_probs);
		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_haps);

		//both_switch
		probs[0] = pow(scaled_dist/num_haps, 2);
		//one switch
		probs[1] =  (((1 - scaled_dist) * scaled_dist) / num_haps) + pow(scaled_dist/num_haps, 2);
		// no switch
		probs[2] = pow(1 - scaled_dist, 2) + ((2 * (1 - scaled_dist) * scaled_dist) / num_haps) + (pow(scaled_dist,2) / pow(num_haps,2));

#pragma omp parallel for schedule(dynamic,32)
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				int chrom_case = ChromosomePair(ml_states_other[m],s).TransitionCase(ChromosomePair(ml_states_other[m-1],j));


				sum += s_forward[m-1][j] * probs[chrom_case];

			}
			s_forward[m][s] =  emission_probs[s] * sum;


		}

		for(int s = 0; s < num_states; s++){
			c += s_forward[m][s];
		}
		normalizers[m] = 1.0/c;

		for(int s = 0; s < num_states; s++){
			s_forward[m][s] = s_forward[m][s] * normalizers[m];
		}
	}

	delete [] emission_probs;

}


/*
 * Calculate scaled forward values for one haplotype (curr_chrom).
 *
 * Transition probabilities calculated using pdf that only considers probability of switch happening along one chromosome.
 * (Assumes independence of recombination along chromosomes)
 *
 * Emission probabilities calculated considering both chromosomes, but with the other chromosomes value fixed at the allele with the highest posterior probability.
 *
 *
 */
void HaplotypePhaser::CalcScaledForwardSeparate(){
	double * emission_probs = new double[num_states];
	double c;

	double scaled_dist;
	double c1,c2;


	InitPriorScaledForward();


	for(int m = 1; m < num_markers; m++){
		//		CalcEmissionProbsIntegrated(m, emission_probs);
		CalcEmissionProbsMla(m, emission_probs);

		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_haps);


		// different (switch happens)
		c2 = scaled_dist / num_haps;
		// same (switch does not happen)
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

		//#pragma omp parallel for schedule(dynamic,32)
#pragma omp parallel for
		for(int s = 0; s < num_states; s++){

			double sum = 0.0;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){


				if (j == s) {
					sum += s_forward[m-1][j] * c1;
				}
				else {
					sum += s_forward[m-1][j] * c2;

				}

			}
			s_forward[m][s] =  emission_probs[s] * sum;
		}

		for(int s = 0; s < num_states; s++){
			c+= s_forward[m][s];
		}
		normalizers[m] = 1.0/c;

		for(int s = 0; s < num_states; s++){
			s_forward[m][s] = s_forward[m][s] * normalizers[m];
		}
	}

	delete [] emission_probs;

}

/*
 * Calculate scaled forward values for one haplotype (curr_chrom).
 *
 * Transition probabilities calculated considering both chromosomes (usual Li & Stephens algo) and integrating over the possible states at the other chromosome,
 * multiplying with their posterior probabilities from the previous iteration.
 *
 * Emission probabilities also calculated considering both chromosomes, but with the other chromosomes value fixed at the state with the highest posterior probability.
 *
 * Not sure this is implemented correctly.
 *
 */
void HaplotypePhaser::CalcScaledForwardIntegrated(){
	double * emission_probs = new double[num_states];
	int * ml_states_other;
	double c;
	double case_probs[3];

	double scaled_dist;
	double c1,c2;

	if (curr_chrom == 2) {
		ml_states_other = ml_states_c1;
	}
	if (curr_chrom == 1) {
		ml_states_other = ml_states_c2;
	}


	InitPriorScaledForward();

	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbs(m, emission_probs);


		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_haps);

		c1 = scaled_dist/num_haps;
		c2 = 1 - scaled_dist;

		//both_switch
		case_probs[0] = pow(c1, 2);
		//one switch
		case_probs[1] =  (c2*c1) + pow(c1, 2);
		// no switch
		case_probs[2] = pow(c2, 2) + (2*c2*c1) + pow(c1, 2);


#pragma omp parallel for schedule(dynamic,32)
		for(int s = 0; s < num_states; s++){

			double sum = 0.0;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				double isum = 0.0;
				// Integrating over possible states on the other chromosome AT MARKER M

				for(int other_chrom_hap_at_m = 0; other_chrom_hap_at_m < num_states; other_chrom_hap_at_m++){

					// Still using ML state for other chromosome at marker m-1

					int chrom_case = ChromosomePair(other_chrom_hap_at_m,s).TransitionCase(ChromosomePair(ml_states_other[m-1],j));
					isum += case_probs[chrom_case] * all_posteriors(m, other_chrom_hap_at_m);

					// Also integrating over possible states at previous marker for other chromosome - slow

					//					for(int other_chrom_hap_at_prev = 0; other_chrom_hap_at_prev < num_states; other_chrom_hap_at_prev++){
					//						int chrom_case = ChromosomePair(other_chrom_hap_at_m,s).NumEquals2(ChromosomePair(other_chrom_hap_at_prev,j));
					//						isum += case_probs[chrom_case] * all_posteriors(m, other_chrom_hap_at_m) *  all_posteriors(m-1, other_chrom_hap_at_prev);
					//					}

				}
				sum += s_forward[m-1][j] * isum;
			}
			s_forward[m][s] =  emission_probs[s] * sum;
		}

		for(int s = 0; s < num_states; s++){
			c+= s_forward[m][s];
		}
		normalizers[m] = 1.0/c;

		for(int s = 0; s < num_states; s++){
			s_forward[m][s] = s_forward[m][s] * normalizers[m];
		}
	}

	delete [] emission_probs;

}

/*
 * Calculate scaled forward values for chromosome 1, marginalized over chromosome 2
 *
 * Transition probabilities calculated using pdf that only considers probability of switch happening along one chromosome.
 * (Assumes independence of recombination along chromosomes)
 *
 * Emission probabilities calculated considering both chromosomes, marginalizing over the possible states at chromosome 2.
 * P(observed genotype | state_at_chrom_1 ) = sum over state_at_chrom_2 : P(observed_genotype | state_at_chrom_1, state_at_chrom_2) * P(state_at_chrom_2)
 *
 *	where P(observed_genotype | state_at_chrom_1, state_at_chrom_2) is the ususal observation probs from Li & S
 *  and
 *  P(state_at_chrom_2) is calculated from above pdf and allele freqs
 *
 *  This is essentially a way to used a prior for chromosome 2 that depends on the allele freqs of the reference?
 *  Not sure its even better than uniform prior for chromosome 2.
 *
 *
 */
void HaplotypePhaser::CalcScaledForwardMarginalized(){
	double * emission_probs = new double[num_states];
	double c, c1, c2;
	double scaled_dist;



	InitPriorScaledForwardMarginalized();

	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbsMarginalized(m, emission_probs);

		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_haps);

		// different
		c2 = scaled_dist / num_haps;
		// same
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

#pragma omp parallel for schedule(dynamic,1)
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;

#pragma GCC ivdep
			for(int j = 0; j < num_states ; j++){

				if (j == s) {
					sum += s_forward[m-1][j] * c1;
				}
				else {
					sum += s_forward[m-1][j] * c2;

				}
			}
			s_forward[m][s] =  emission_probs[s] * sum;
		};


		for(int s = 0; s < num_states; s++){
			c+= s_forward[m][s];
		}

		// c can become 0 there are GLs that are value 0. GLs 0 should not be allowed.
		normalizers[m] = 1.0/c;

		for(int s = 0; s < num_states; s++){
			s_forward[m][s] = s_forward[m][s] * normalizers[m];
		};
	};

	delete [] emission_probs;
}


/*
 * Calculate scaled backward values for one chromosome (curr_chrom).
 *
 * Transition probabilities calculated considering both chromosomes (usual Li & Stephens algo) but with the other chromosome fixed at the state that was most
 * likely in the previous iteration.
 *
 * Emission probabilities also calculated considering both chromosomes, but with the other chromosomes value fixed at the state with the highest posterior probability.
 *
 *
 */
void HaplotypePhaser::CalcScaledBackwardFixed(){
	double * emission_probs = new double[num_states];
	int * ml_states_other;
	double probs[3];
	double scaled_dist;

	InitPriorScaledBackward();

	if (curr_chrom == 2) {
		ml_states_other = ml_states_c1;
	}
	if (curr_chrom == 1) {
		ml_states_other = ml_states_c2;
	}

	for(int m = num_markers-2; m >= 0; m--){
		scaled_dist = 1-exp(-(distances[m+1] * pop_const)/num_haps);

		//both_switch
		probs[0] = pow(scaled_dist/num_haps, 2);
		//one switch
		probs[1] =  (((1 - scaled_dist) * scaled_dist) / num_haps) + pow(scaled_dist/num_haps, 2);
		// no switch
		probs[2] = pow(1 - scaled_dist, 2) + ((2 * (1 - scaled_dist) * scaled_dist) / num_haps) + (pow(scaled_dist,2) / pow(num_haps,2));

		CalcEmissionProbs(m+1, emission_probs);

#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;

#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				int chrom_case = ChromosomePair(ml_states_other[m],s).TransitionCase(ChromosomePair(ml_states_other[m+1],j));

				sum +=  s_backward[m+1][j]* probs[chrom_case] * emission_probs[j];

			}
			s_backward[m][s] = sum * normalizers[m];
		}
	}
	delete [] emission_probs;
}



/*
 * Calculate scaled backward values for chromosome 1, marginalized over chromosome 2
 *
 * Transition probabilities calculated using pdf that only considers probability of switch happening along one chromosome.
 * (Assumes independence of recombination along chromosomes)
 *
 * Emission probabilities calculated considering both chromosomes, marginalizing over the possible states at chromosome 2.
 * P(observed genotype | state_at_chrom_1 ) = sum over state_at_chrom_2 : P(observed_genotype | state_at_chrom_1, state_at_chrom_2) * P(state_at_chrom_2)
 *
 *	where P(observed_genotype | state_at_chrom_1, state_at_chrom_2) is the ususal observation probs from Li & S
 *  and
 *  P(state_at_chrom_2) is calculated from above pdf and allele freqs
 *
 *  This is essentially a way to used a prior for chromosome 2 that depends on the allele freqs of the reference?
 *  Not sure its even better than uniform prior for chromosome 2.
 *
 *
 */
void HaplotypePhaser::CalcScaledBackwardMarginalized(){
	double * emission_probs = new double[num_states];

	double scaled_dist;
	double c1;
	double c2;



	InitPriorScaledBackward();


	for(int m = num_markers-2; m >= 0; m--){
		scaled_dist = 1.0 - exp(-(distances[m+1] * pop_const)/num_haps);

		// different
		c2 = scaled_dist / num_haps;
		// same
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

		CalcEmissionProbsMarginalized(m+1, emission_probs);


#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				if (s == j) {
					sum +=  s_backward[m+1][j]* c1 * emission_probs[j];
				}
				else {
					sum +=  s_backward[m+1][j]* c2 * emission_probs[j];

				}
			}
			s_backward[m][s] = sum * normalizers[m];

		}
	}
	delete [] emission_probs;
}



/*
 * Calculate scaled backward values for one haplotype (curr_chrom).
 *
 * Transition probabilities calculated considering both chromosomes (usual Li & Stephens algo) and integrating over the possible states at the other chromosome,
 * multiplying with their posterior probabilities from the previous iteration.
 *
 * Emission probabilities also calculated considering both chromosomes, but with the other chromosomes value fixed at the state with the highest posterior probability.
 *
 * Not sure this is implemented correctly.
 *
 */
void HaplotypePhaser::CalcScaledBackwardIntegrated(){
	double * emission_probs = new double[num_states];
	int num_h = 2*num_inds - 2;
	int * ml_states_other;

	double case_probs[3];


	if (curr_chrom == 2) {
		ml_states_other = ml_states_c1;
	}
	if (curr_chrom == 1) {
		ml_states_other = ml_states_c2;
	}



	double scaled_dist;
	double c1;
	double c2;



	InitPriorScaledBackward();

	for(int m = num_markers-2; m >= 0; m--){
		scaled_dist = 1.0 - exp(-(distances[m+1] * pop_const)/num_h);

		c1 = scaled_dist/num_h;
		c2 = 1 - scaled_dist;

		//both_switch
		case_probs[0] = pow(c1, 2);
		//one switch
		case_probs[1] =  (c2*c1) + pow(c1, 2);
		// no switch
		case_probs[2] = pow(c2, 2) + (2*c2*c1) + pow(c1, 2);

		CalcEmissionProbs(m+1, emission_probs);


#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;



#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				//				int chrom_case = (states[s]).NumEquals(states[j]);
				//				if(states[s].first == states[s].second) {
				//					chrom_case += 3;
				//				}

				double isum = 0.0;
				for(int other_chrom_hap_at_m = 0; other_chrom_hap_at_m < num_states; other_chrom_hap_at_m++){

					//stm1
					int chrom_case = ChromosomePair(other_chrom_hap_at_m,s).TransitionCase(ChromosomePair(ml_states_other[m+1],j));
					isum += case_probs[chrom_case] * all_posteriors(m,other_chrom_hap_at_m);

					//stm2
					//					for(int other_chrom_hap_at_next = 0; other_chrom_hap_at_next < num_states; other_chrom_hap_at_next++){
					//						int chrom_case = ChromosomePair(other_chrom_hap_at_m,s).NumEquals2(ChromosomePair(other_chrom_hap_at_next,j));
					//						isum += case_probs[chrom_case] * all_posteriors(m,other_chrom_hap_at_m) * all_posteriors(m+1,other_chrom_hap_at_next);
					//					}
				}
				sum +=  s_backward[m+1][j] * emission_probs[j] * isum;
			}
			s_backward[m][s] = sum * normalizers[m];
		}
	}
	delete [] emission_probs;
}


/*
 * Calculate scaled backward values for one haplotype (curr_chrom).
 *
 * Transition probabilities calculated using pdf that only considers probability of switch happening along one chromosome.
 * (Assumes independence of recombination along chromosomes)
 *
 * Emission probabilities calculated considering both chromosomes, but with the other chromosomes value fixed at the allele with the highest posterior probability.
 *
 *
 */
void HaplotypePhaser::CalcScaledBackwardSeparate(){
	double * emission_probs = new double[num_states];
	int num_h = 2*num_inds - 2;

	double scaled_dist;
	double c1;
	double c2;



	InitPriorScaledBackward();

	for(int m = num_markers-2; m >= 0; m--){
		scaled_dist = 1.0 - exp(-(distances[m+1] * pop_const)/num_h);

		// different (switch happens)
		c2 = scaled_dist / num_haps;
		// same (no switch)
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

		//		CalcEmissionProbsIntegrated(m+1, emission_probs);
		CalcEmissionProbsMla(m+1, emission_probs);


#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;

#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				if (s == j) {
					sum +=  s_backward[m+1][j]* c1 * emission_probs[j];
				}
				else {
					sum +=  s_backward[m+1][j]* c2 * emission_probs[j];

				}

			}

			s_backward[m][s] = sum * normalizers[m];

		}
	}
	delete [] emission_probs;
}



/**
 * For every marker m, get the state with highest posterior probability at that location, given the entire observation sequence
 * (not most likely sequence of states)
 *
 * for chromosome 1, marginalized over chromosome 2
 *
 *
 * Calculate array ml_states s.t.
 * ml_states[m] = s_i that maximises P(Q_m = s_i | O_1 ... O_num_markers)
 *
 * for m in {0 ... num_markers-1}
 *
 *
 *
 * print: print stats to filename, where stats contains one row for every marker
 *
 *
 * returns stats:
 * every row has 44 values:
 * 10 lowest posterior allele probabilities for chromosome 1
 * 10 states corresponding to lowest gallele probabilities for chromosome 1
 * to highest posterior allele probabilities for chromosome 1
 * 10 states corresponding to higest allele probabilities for chromosome 1
 * sum of posterior probs for all states at this marker / num states (expect this to be 1.0 / num_states)
 * posterior allele probability of allele 0
 * posterior allele probability of allele 1
 * -1
 *
 *
 * also fills in ml_states_h1 with the most likely state of chromosome 1 for each marker
 * and
 * also fills in ml_alles_h1 with the most likely allele for chromosome 1 for each marker
 * (allele probs calculated by summing over posterior state probs for all states implying a give allele)
 * and
 * also fills in all_posteriors with posterior probabilities of each state for each marker for chromosome 1
 *
 */
vector<vector<double>>  HaplotypePhaser::GetPosteriorStatsMarginalized(const char * filename, bool print){
	vector<vector<double>> stats;
	vector<vector<double>> allele_probs;
	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(44,-1.0);
	};

	for(int m = 0; m < num_markers; m++) {
		allele_probs.push_back({});
		allele_probs[m].resize(2,0.0);
	};

	for(int m = 0; m < num_markers; m++) {

		vector<double> posteriors;
		posteriors.resize(num_states, -1);

		double norm = 0.0;

		for(int i = 0; i < num_states; i++) {
			norm += s_forward[m][i] * s_backward[m][i];
		}

		double sum = 0.0;

		for(int s = 0; s < num_states; s++) {

			stringstream tmp;
			tmp << setprecision(prob_precision) << fixed << s_forward[m][s] * s_backward[m][s] / norm;
			posteriors[s] = stod(tmp.str());

			all_posteriors(m,s) = posteriors[s];


			sum += posteriors[s];

			int allele_code = haplotypes(s,m);
			allele_probs[m][allele_code] += posteriors[s];

		}

		vector<size_t> res = VcfUtils::sort_indexes(posteriors);


		//add lowest probabilities to stats[m] - in increasing order
		for(int i = 0; i < 10; i++) {
			stats[m][i] = posteriors[res[i]];
		}

		//add states corresponding to lowest probabilitoes to stats[m]
		for(int i = 0; i < 10; i++) {
			stats[m][10 + i] = res[i];
		}

		//add highest probabilities to stats[m]  - in increasing order (ml state gets index 39)
		for(int i = 0; i < 10; i++) {
			stats[m][20 + i] = posteriors[res[posteriors.size() - 10 + i]];
		}


		//add states corresponding to highest probabilitoes to stats[m]
		for(int i = 0; i < 10; i++) {
			stats[m][30 + i] = res[posteriors.size() - 10 + i];
		}

		stats[m][40] = sum / posteriors.size();
		stats[m][41] =  allele_probs[m][0];
		stats[m][42] =  allele_probs[m][1];
		stats[m][43] = -1;

		ml_states_c1[m] = stats[m][39];

		if(stats[m][41] >= stats[m][42]) {
			ml_alleles_c1[m] = 0;
		}
		else {
			ml_alleles_c1[m] = 1;
		}

	};




	if(print) {
		VcfUtils::writeVectorToCSV(filename, stats, "w");
	};

	return stats;
}



/**
 * For every marker m, get the state with highest posterior probability at that location, given the entire observation sequence
 * (not most likely sequence of states)
 *
 * Calculate array ml_states s.t.
 * ml_states[m] = s_i that maximises P(Q_m = s_i | O_1 ... O_num_markers)
 *
 * for m in {0 ... num_markers-1}
 *
 * print: print stats to filename, where stats contains one row for every marker
 *
 * ml_state: true:  use maximum likelihood state of other chromosome when calculating posterior genoytpe probabilities
 * 			 false: use maximum likelihood allele of other chromosome when calculating posterior genoytpe probabilities
 *
 * return stats:
 * every row has 44 values:
 * 10 lowest posterior genotype probabilities
 * 10 states corresponding to lowest genotype probabilities
 * to highest posterior genotype probabilities
 * 10 states corresponding to higest genotype probabilities
 * sum of posterior probs for all states at this marker / num states (expect this to be 1.0 / num_states)
 * posterior probability of genotype 0
 * posterior probability of genotype 1
 * posterior probability of genotype 2
 *
 *
 * current chromosome is the one implied by curr_chrom
 *
 * also fills in ml_states_h{current} with the most likely state of current chromosome for each marker
 * and
 * also fills in ml_alles_h{current} with the most likely allele for current chromosome for each marker
 * (allele probs calculated by summing over posterior state probs for all states implying a give allele)
 * and
 * also fills in all_posteriors with posterior probabilities of each state for each marker for current chromosome
 *
 *
 */
vector<vector<double>>  HaplotypePhaser::GetPosteriorStats(const char * filename, bool ml_state, bool print){
	vector<vector<double>> stats;
	vector<vector<double>> geno_probs;
	vector<vector<double>> allele_probs;

	int * ml_alleles_other;
	int * ml_alleles_this;
	int * ml_states_this;
	int * ml_states_other;

	if (curr_chrom == 2) {
		ml_alleles_other = ml_alleles_c1;
		ml_alleles_this = ml_alleles_c2;

		ml_states_other = ml_states_c1;
		ml_states_this = ml_states_c2;

	}
	if (curr_chrom == 1) {
		ml_alleles_other = ml_alleles_c2;
		ml_alleles_this = ml_alleles_c1;

		ml_states_other = ml_states_c2;
		ml_states_this = ml_states_c1;
	}



	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(44,-1.0);
	}

	for(int m = 0; m < num_markers; m++) {
		geno_probs.push_back({});
		geno_probs[m].resize(3,0.0);
	}

	for(int m = 0; m < num_markers; m++) {
		allele_probs.push_back({});
		allele_probs[m].resize(2,0.0);
	};




	for(int m = 0; m < num_markers; m++) {
		vector<double> posteriors;
		posteriors.resize(num_states, -1);

		double norm = 0.0;
		for(int i = 0; i < num_states; i++) {
			norm += s_forward[m][i] * s_backward[m][i];
		}

		double sum = 0.0;
		for(int s = 0; s < num_states; s++) {


			//			posteriors[s] = s_forward[m][s] * s_backward[m][s] / norm;
			stringstream tmp;
			tmp << setprecision(prob_precision) << fixed << s_forward[m][s] * s_backward[m][s] / norm;
			posteriors[s] = stod(tmp.str());

			all_posteriors(m,s) = posteriors[s];
			sum += posteriors[s];



			//////////genotype probability - retrieve which genotype this pair of ref haps implies /////////////////

			int ref_hap_this = s;

			// AGCT allele
			//			String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
			//			String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

			// 00, 01, 10, 11

			int hapcode_other;
			if (ml_state){
				// use ml state for genos
				hapcode_other = haplotypes(ml_states_other[m],m);

			} else{
				// use ml allele for genos
				hapcode_other = ml_alleles_other[m];
			}


			int hapcode_this = haplotypes(ref_hap_this,m);

			allele_probs[m][hapcode_this] += posteriors[s];




			int geno_code;

			if(hapcode_other != hapcode_this) {
				geno_code = 1;
			}
			else {
				geno_code = (hapcode_other == 0) ? 0 : 2;
			}


			geno_probs[m][geno_code] += posteriors[s];

		}

		// check that the sum of genotype probabilities adds up to 1
		float check_sum = 0.0;
		for(int i = 0; i < 3 ; i++) {
			check_sum += geno_probs[m][i];
		}
		if(abs(check_sum - 1.0) > 0.000001 ) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!Sum of all geno probs is %f at marker %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1\n ", check_sum, m);
		}


		vector<size_t> res = VcfUtils::sort_indexes(posteriors);

		//add lowest elements to stats[m]
		for(int i = 0; i < 10; i++) {
			stats[m][i] = posteriors[res[i]];
		}
		for(int i = 0; i < 10; i++) {
			stats[m][10 + i] = res[i];
		}
		for(int i = 0; i < 10; i++) {
			stats[m][20 + i] = posteriors[res[posteriors.size() - 10 + i]];
		}
		for(int i = 0; i < 10; i++) {
			stats[m][30 + i] = res[posteriors.size() - 10 + i];
		}

		stats[m][40] = sum / posteriors.size();
		stats[m][41] = geno_probs[m][0];
		stats[m][42] = geno_probs[m][1];
		stats[m][43] = geno_probs[m][2];

		ml_states_this[m] = stats[m][39];


		if(allele_probs[m][0] >= allele_probs[m][1]) {
			ml_alleles_this[m] = 0;
		}
		else {
			ml_alleles_this[m] = 1;
		}

	}


	if (print) {
		printf("Writing stats to %s \n", filename);
		VcfUtils::writeVectorToCSV(filename, stats, "w");

	}

	return stats;
}




/**
 * Read stats from filename
 */
vector<vector<double>>  HaplotypePhaser::ReadPosteriorStats(const char * filename){

	vector<vector<double>> stats;


	FILE * statstream = fopen(filename, "r");

	printf("Opened statsream for %s \n", filename);
	float low1;
	float low2;
	float low3;
	float low4;
	float low5;
	float low6;
	float low7;
	float low8;
	float low9;
	float low10;

	float low1_i;
	float low2_i;
	float low3_i;
	float low4_i;
	float low5_i;
	float low6_i;
	float low7_i;
	float low8_i;
	float low9_i;
	float low10_i;

	float high1;
	float high2;
	float high3;
	float high4;
	float high5;
	float high6;
	float high7;
	float high8;
	float high9;
	float high10;

	float high1_i;
	float high2_i;
	float high3_i;
	float high4_i;
	float high5_i;
	float high6_i;
	float high7_i;
	float high8_i;
	float high9_i;
	float high10_i;

	float average;
	int result;


	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(41,-1.0);
	}

	for(int m = 0; m < num_markers; m++) {
		printf("scanning for marker = %d \n", m);
		result = fscanf(statstream, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\n",
				&low1,&low2,&low3,&low4,&low5,&low6,&low7,&low8,&low9,&low10,&low1_i,&low2_i,&low3_i,&low4_i,&low5_i,&low6_i,&low7_i,&low8_i,&low9_i,&low10_i,
				&high1,&high2,&high3,&high4,&high5,&high6,&high7,&high8,&high9,&high10,&high1_i,&high2_i,&high3_i,&high4_i,&high5_i,&high6_i,&high7_i,&high8_i,&high9_i,&high10_i,&average);


		printf("result = %d \n", result);

		printf("scanned for marker = %d \n", m);
		stats[m][0] = low1;
		stats[m][1] = low2;
		stats[m][2] = low3;
		stats[m][3] = low4;
		stats[m][4] = low5;
		stats[m][5] = low6;
		stats[m][6] = low7;
		stats[m][7] = low8;
		stats[m][8] = low9;
		stats[m][9] = low10;

		stats[m][10] = low1_i;
		stats[m][11] = low2_i;
		stats[m][12] = low3_i;
		stats[m][13] = low4_i;
		stats[m][14] = low5_i;
		stats[m][15] = low6_i;
		stats[m][16] = low7_i;
		stats[m][17] = low8_i;
		stats[m][18] = low9_i;
		stats[m][19] = low10_i;

		stats[m][20] = high1;
		stats[m][21] = high2;
		stats[m][22] = high3;
		stats[m][23] = high4;
		stats[m][24] = high5;
		stats[m][25] = high6;
		stats[m][26] = high7;
		stats[m][27] = high8;
		stats[m][28] = high9;
		stats[m][29] = high10;

		stats[m][30] = high1_i;
		stats[m][31] = high2_i;
		stats[m][32] = high3_i;
		stats[m][33] = high4_i;
		stats[m][34] = high5_i;
		stats[m][35] = high6_i;
		stats[m][36] = high7_i;
		stats[m][37] = high8_i;
		stats[m][38] = high9_i;
		stats[m][39] = high10_i;

		stats[m][40] = average;


	}
	return stats;
}


/**
 * Print the given genotypes to a VCF.
 *
 */
void HaplotypePhaser::PrintGenotypesToVCF(vector<vector<int>> & genotypes, const char * out_file, const char * sample_file, const char * vcf_template ){
	VcfRecord record_template;
	VcfFileReader reader;
	VcfHeader header_read;

	reader.open(vcf_template, header_read);
	reader.readRecord(record_template);

	int num_samples = header_read.getNumSamples();

	if (num_samples != record_template.getNumSamples()) {
		printf("!!!!! INCONSISTENT NUMBER OF SAMPLES when writing to VCF !!!!!\nSomething is probably wrong with the vcf template \n");
		printf("num samples from header:  %d \n", num_samples);
		printf("num samples in record:  %d \n", record_template.getNumSamples());
	}
	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_read, InputFile::BGZF);

	for(int m = 0; m < num_markers; m++) {

		MarkerInfo* markerinfo = Pedigree::GetMarkerInfo(m);
		std::string marker_name = markerinfo->name.c_str();
		std::size_t delim = marker_name.find(":");
		string chrom =  marker_name.substr(0,delim);
		int pos =  std::stoi(marker_name.substr(delim+1));
		record_template.setChrom(chrom.c_str());
		record_template.set1BasedPosition(pos);
		record_template.setID(marker_name.c_str());
		record_template.setRef((markerinfo->GetAlleleLabel(1)).c_str());
		record_template.setAlt((markerinfo->GetAlleleLabel(2)).c_str());
		record_template.setQual(".");

		for(int sample = 0; sample < num_samples; sample++) {

			int succ;

			//			if (genotypes[sample][m] == 0) {
			//				succ = record_template.getGenotypeInfo().setString("GL",sample, "1,0,0");
			//			}
			//			if (genotypes[sample][m] == 1) {
			//				succ = record_template.getGenotypeInfo().setString("GL",sample, "0,1,0");
			//			}
			//			if (genotypes[sample][m] == 2) {
			//				succ = record_template.getGenotypeInfo().setString("GL",sample, "0,0,1");
			//			}

			if (genotypes[sample][m] == 0) {
				succ = record_template.getGenotypeInfo().setString("GT",sample, "0/0");
			}
			if (genotypes[sample][m] == 1) {
				succ = record_template.getGenotypeInfo().setString("GT",sample, "0/1");
			}
			if (genotypes[sample][m] == 2) {
				succ = record_template.getGenotypeInfo().setString("GT",sample, "1/1");
			}


			if (!succ) {
				printf("ERROR IN WRITING TO VCF for marker %d for ind %d \n", m, sample);
			}
		}

		writer.writeRecord(record_template);

	}

	printf("WROTE GENOS TO %s \n", out_file);

	reader.close();
	writer.close();

}




/**
 *
 * Print the phased haplotypes indicated by the maximum likelihood states ml_states for each chromosome to vcf.
 *
 * ml_states_h{z}: n_samples x n_markers vector of most likely states for each chromosome z
 *
 */

void HaplotypePhaser::PrintMLSHaplotypesToVCF(vector<vector<int>> & ml_states_h1,vector<vector<int>> & ml_states_h2,  const char * out_file, const char * sample_file, const char * vcf_template){

	std::vector<String> h1;
	std::vector<String> h2;

	int ref_hap1;
	int ref_hap2;


	VcfRecord record_template;
	VcfFileReader reader;
	VcfHeader header_read;

	reader.open(vcf_template, header_read);
	reader.readRecord(record_template);


	int num_samples = header_read.getNumSamples();

	if (num_samples != record_template.getNumSamples()) {
		printf("!!!!! INCONSISTENT NUMBER OF SAMPLES when writing to VCF !!!!!\nSomething is probably wrong with the vcf template \n");
		printf("num samples from header:  %d \n", num_samples);
		printf("num samples in record:  %d \n", record_template.getNumSamples());
	}

	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_read, InputFile::BGZF);

	for(int m = 0; m < num_markers; m++) {
		MarkerInfo* markerinfo = Pedigree::GetMarkerInfo(m);
		std::string marker_name = markerinfo->name.c_str();
		std::size_t delim = marker_name.find(":");
		string chrom =  marker_name.substr(0,delim);
		int pos =  std::stoi(marker_name.substr(delim+1));
		record_template.setChrom(chrom.c_str());
		record_template.set1BasedPosition(pos);
		record_template.setID(marker_name.c_str());
		record_template.setRef((markerinfo->GetAlleleLabel(1)).c_str());
		record_template.setAlt((markerinfo->GetAlleleLabel(2)).c_str());
		record_template.setQual(".");

		for(int sample = 0; sample < num_samples; sample++) {

			int succ;

			ref_hap1 = ml_states_h1[sample][m];
			ref_hap2 = ml_states_h2[sample][m];;


			std::stringstream ss;
			ss << to_string(int(haplotypes(ref_hap1,m))) << "|" << to_string(int(haplotypes(ref_hap2,m)));
			string GTstring = ss.str();

			succ = record_template.getGenotypeInfo().setString("GT",sample, GTstring);

			if (!succ) {
				printf("ERROR IN WRITING TO VCF for marker %d for ind %d \n", m, sample);
			}
		}

		writer.writeRecord(record_template);

	}

	printf("WROTE PHASED (most likely state) TO %s \n", out_file);

	reader.close();
	writer.close();

}


/**
 *
 * Print the phased haplotypes indicated by the maximum likelihood alleles ml_alleles for each chromosome to vcf.
 *
 * ml_alleles_h{z}: n_samples x n_markers vector of most likely alleles for each chromosome z
 *
 */
void HaplotypePhaser::PrintMLAHaplotypesToVCF(vector<vector<int>> & ml_alleles_h1, vector<vector<int>> & ml_alleles_h2,  const char * out_file, const char * sample_file, const char * vcf_template){

	std::vector<String> h1;
	std::vector<String> h2;


	VcfRecord record_template;
	VcfFileReader reader;
	VcfHeader header_read;

	reader.open(vcf_template, header_read);
	reader.readRecord(record_template);


	int num_samples = header_read.getNumSamples();

	if (num_samples != record_template.getNumSamples()) {
		printf("!!!!! INCONSISTENT NUMBER OF SAMPLES when writing to VCF !!!!!\nSomething is probably wrong with the vcf template \n");
		printf("num samples from header:  %d \n", num_samples);
		printf("num samples in record:  %d \n", record_template.getNumSamples());
	}

	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_read, InputFile::BGZF);

	for(int m = 0; m < num_markers; m++) {
		MarkerInfo* markerinfo = Pedigree::GetMarkerInfo(m);
		std::string marker_name = markerinfo->name.c_str();
		std::size_t delim = marker_name.find(":");
		string chrom =  marker_name.substr(0,delim);
		int pos =  std::stoi(marker_name.substr(delim+1));
		record_template.setChrom(chrom.c_str());
		record_template.set1BasedPosition(pos);
		record_template.setID(marker_name.c_str());
		record_template.setRef((markerinfo->GetAlleleLabel(1)).c_str());
		record_template.setAlt((markerinfo->GetAlleleLabel(2)).c_str());
		record_template.setQual(".");

		for(int sample = 0; sample < num_samples; sample++) {

			int succ;

			std::stringstream ss;
			ss << to_string(ml_alleles_h1[sample][m]) << "|" << to_string(ml_alleles_h2[sample][m]);
			string GTstring = ss.str();

			succ = record_template.getGenotypeInfo().setString("GT",sample, GTstring);

			if (!succ) {
				printf("ERROR IN WRITING TO VCF for marker %d for ind %d \n", m, sample);
			}
		}

		writer.writeRecord(record_template);

	}

	printf("WROTE PHASED (most likely allele) TO %s \n", out_file);

	reader.close();
	writer.close();

}
