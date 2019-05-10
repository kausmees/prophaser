
#include <string.h>
#include <omp.h>

#include "HaplotypePhaserSym.h"
#include "MemoryAllocators.h"




HaplotypePhaserSym::~HaplotypePhaserSym(){
	//	FreeCharMatrix(haplotypes, ped.count*2);
	//	FreeCharMatrix(genotypes, ped.count);
	delete [] phred_probs;
	FreeDoubleMatrix(s_forward, num_markers);
	FreeDoubleMatrix(s_backward, num_markers);

	//TODO need to delete ped?

}


void HaplotypePhaserSym::AllocateMemory(){
	String::caseSensitive = false;

	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	num_haps = (ped.count-1)*2;

	num_states = ((num_haps)*(num_haps+1)) / 2;

	for(int i = 0; i < num_haps; i++) {
		for(int j = i; j < num_haps; j++) {
			states.push_back(ChromosomePair(i,j));
		}
	}


	num_markers = Pedigree::markerCount;

	// Total number of individuals. Samples + references.
	num_inds = ped.count;

	//	error = 0.01;
	//	theta = 0.01;

	distances.resize(num_markers,0.01);
	//	errors.resize(num_markers,0.01);

	phred_probs = new float[256];

	for(int i = 0; i < 256; i++){
		phred_probs[i] = pow(10.0, -i * 0.1);
	};

	//	haplotypes = AllocateCharMatrix(num_inds*2, num_markers);
	//	haplotypes = MatrixXi(num_inds*2, num_markers);

	// haplotypes(h,m) = 0/1 representing allele of haplotype h at marker m
	// number of rows is num_haps + 2 (if there is one sample)
	haplotypes = MatrixXc(num_inds*2, num_markers);

	// allele_counts(m) = number of haplotypes with allele 1 at marker m
	allele_counts = VectorXi(num_markers);

	L = MatrixXd(6,5);
	u = VectorXd(6);
	u << 0, 0, 0 , 0 , 0, 1;
	x = VectorXd(5);
	coefs_p0 = VectorXd(5);
	coefs_p1 = VectorXd(5);
	sample_gls.resize(num_markers*3);

	all_posteriors = MatrixXd(num_markers,num_states);


	s_forward = AllocateDoubleMatrix(num_markers, num_states);
	s_backward = AllocateDoubleMatrix(num_markers, num_states);
	normalizers = new double[num_markers];
	ml_states_h1 = new int[num_markers];
	ml_states_h2 = new int[num_markers];

	ml_alleles_h1 = new int[num_markers];
	ml_alleles_h2 = new int[num_markers];


};

void HaplotypePhaserSym::DeAllocateMemory(){
	//	FreeCharMatrix(haplotypes, ped.count*2);
};
void HaplotypePhaserSym::setDistanceCode(int c) {
	distance_code = c;
};

/**
 * Load vcf data from files.
 * Reference haplotypes are loaded from ref_file.
 * Unphased sample is loaded from sample_file, the individual with sample_index.
 *
 */
void HaplotypePhaserSym::LoadData(const String &ref_file, const String &sample_file, int sample_index, const String &map_file){
	VcfUtils::LoadReferenceMarkers(ref_file);
	printf("1\n");
	VcfUtils::LoadReferenceIndividuals(ped, ref_file);
	VcfUtils::LoadSampleIndividual(ped, sample_file, sample_index);

	printf("2\n");
	AllocateMemory();
	printf("3\n");
	//	CalcCases();
	printf("4\n");
	// Loads the reference haplotypes
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);

	// (Fill the initial sample haplotypes - not used for single-sample: but prep for multi-sample)
	VcfUtils::FillSampleHaplotypes(ped, haplotypes, num_haps);



	printf("5\n");
	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);
	printf("6\n");
	//	VcfUtils::LoadGeneticMap(map_file, ped, distances);
	printf("7\n");

};

/**
 * Load vcf data from files. Only reference data.
 *
 */
void HaplotypePhaserSym::LoadReferenceData(const String &ref_file){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadReferenceIndividuals(ped,ref_file);
	AllocateMemory();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);
};

/**
 * Load vcf data from files. Only reference data.
 *
 */
void HaplotypePhaserSym::LoadSampleData(const String &sample_file, int sample_index){
	VcfUtils::LoadSampleIndividual(ped, sample_file, sample_index);
	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);
	//	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
	//		VcfUtils::LoadGeneticMap("data/chr20.OMNI.interpolated_genetic_map", ped, distances);

};


//void HaplotypePhaserSym::CalcCases(){
//	for(int j = 0; j < num_states; j ++) {
//		for(int s = 0; s < num_states; s ++) {
//			int chrom_case = (states[s]).NumEquals2(states[j]);
//			if(states[j].first == states[j].second) {
//				chrom_case += 3;
//			}
//			pdf_cases(j,s) = chrom_case;
//		}
//	}
//
////	pdf_cases_rowmaj = pdf_cases;
//
////	std::cout << "Pdf cases: m:\n" << pdf_cases << std::endl;
//
//};


/**
 * Fill allele_counts.
 *
 * allele_counts(m) = number of haplotypes with allele 1 at marker m
 *
 */
void HaplotypePhaserSym::CalcAlleleCounts() {
	allele_counts = haplotypes.colwise().sum();
};



/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * forall hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 * for fixed reference haplotype at chromosome 1
 */

void HaplotypePhaserSym::CalcEmissionProbs(int marker, double * probs) {
	int h1;
	int h2;
	double sum;

	double case_1 = (pow(1 - error, 2) + pow(error, 2));
	double case_2 = 2 * (1 - error) * error;
	double case_3 = pow(1 - error, 2);
	double case_4 = (1 - error) * error;
	double case_5 = pow(error,2);

	int * ml_states_other;

	if (curr_hap == 2) {
		ml_states_other = ml_states_h1;
	}
	if (curr_hap == 1) {
		ml_states_other = ml_states_h2;
	}


	// OPT1
	for (int state = 0; state < num_states; state++) {

		// Reference hapotype at chromosome 1 - fixed (0: REF 1: ALT)
		h1 = haplotypes(ml_states_other[marker], marker);

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

void HaplotypePhaserSym::CalcEmissionProbsMla(int marker, double * probs) {
	int h1;
	int h2;
	double sum;

	double case_1 = (pow(1 - error, 2) + pow(error, 2));
	double case_2 = 2 * (1 - error) * error;
	double case_3 = pow(1 - error, 2);
	double case_4 = (1 - error) * error;
	double case_5 = pow(error,2);

	int * ml_alleles_other;

	if (curr_hap == 2) {
		ml_alleles_other = ml_alleles_h1;

	}
	if (curr_hap == 1) {
		ml_alleles_other = ml_alleles_h2;

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

void HaplotypePhaserSym::CalcEmissionProbsIntegrated(int marker, double * probs) {
	int h1;
	int h2;
	double sum;
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

	if (curr_hap == 2) {
		ml_states_other = ml_states_h1;
	}
	if (curr_hap == 1) {
		ml_states_other = ml_states_h2;
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
				//				sum_states += case_3 *  sample_gls[marker * 3];
				//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
			}
			else {
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno0 += case_4 * this_p;
					//					sum_states += case_4 * sample_gls[marker * 3];
					//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

				}
				else{
					sum_geno0 += case_5 * this_p;
					//					sum_states +=  case_5 *  sample_gls[marker * 3];
					//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

				}
			}

			// case2: g = 1
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum_geno1 += case_1 * this_p;
				//				sum_states += case_1 *  sample_gls[marker * 3 + 1];
				//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

			}
			else{
				sum_geno1 += case_2 * this_p;
				//				sum_states +=case_2 * sample_gls[marker * 3 + 1];
				//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

			}

			// case3: g = 2
			if(h1 == 1 and h2 == 1){
				sum_geno2 += case_3 * this_p;
				//				sum_states += case_3 * sample_gls[marker * 3 + 2];
				//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			} else{
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno2 += case_4 * this_p;
					//					sum_states += case_4 * sample_gls[marker * 3 + 2];
					//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

				}
				else{
					sum_geno2 += case_5 * this_p;
					//					sum_states += case_5 * sample_gls[marker * 3 + 2];
					//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

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
 *  Fill marginal_state_probs with p(s) marginalized over emission probabilities for marker marker.
 *
 */

void HaplotypePhaserSym::CalcMarginalStateProbs(int marker) {

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

void HaplotypePhaserSym::CalcEmissionProbsMarginalized(int marker, double * probs) {
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

			this_p = h2*p1 + (1-h2)*p0;

			// case1: g = 0

			if(h1 == 0 and h2 == 0){
				sum_geno0 += case_3 * this_p;
				//				sum_states += case_3 *  sample_gls[marker * 3];
				//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
			}
			else {
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno0 += case_4 * this_p;
					//					sum_states += case_4 * sample_gls[marker * 3];
					//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

				}
				else{
					sum_geno0 += case_5 * this_p;
					//					sum_states +=  case_5 *  sample_gls[marker * 3];
					//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

				}
			}

			// case2: g = 1
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum_geno1 += case_1 * this_p;
				//				sum_states += case_1 *  sample_gls[marker * 3 + 1];
				//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

			}
			else{
				sum_geno1 += case_2 * this_p;
				//				sum_states +=case_2 * sample_gls[marker * 3 + 1];
				//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

			}

			// case3: g = 2
			if(h1 == 1 and h2 == 1){
				sum_geno2 += case_3 * this_p;
				//				sum_states += case_3 * sample_gls[marker * 3 + 2];
				//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			} else{
				if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
					sum_geno2 += case_4 * this_p;
					//					sum_states += case_4 * sample_gls[marker * 3 + 2];
					//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

				}
				else{
					sum_geno2 += case_5 * this_p;
					//					sum_states += case_5 * sample_gls[marker * 3 + 2];
					//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

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




void HaplotypePhaserSym::InitPriorScaledForward(){
	double * emission_probs = new double[num_states];
	double prior = 1.0 / num_states;
	double c1 = 0.0;


	//	CalcEmissionProbsIntegrated(0, emission_probs);
	CalcEmissionProbsMla(0, emission_probs);


	//		printf("Emission probs[0] h2: \n");
	//		for(int s = 0; s < num_states; s++){
	//			printf("%e \n", emission_probs[s]);
	//		};


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


void HaplotypePhaserSym::InitPriorScaledForwardMarginalized(){
	double * emission_probs = new double[num_states];
	double prior = 1.0 / num_states;
	double c1 = 0.0;


	CalcEmissionProbsMarginalized(0, emission_probs);

	//
	//	printf("Emission probs[0] h1: \n");
	//	for(int s = 0; s < num_states; s++){
	//		printf("%e \n", emission_probs[s]);
	//	};


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



void HaplotypePhaserSym::InitPriorScaledBackward(){

	for(int s = 0; s < num_states; s++){
		s_backward[num_markers-1][s] = normalizers[num_markers-1];
	};
};


void HaplotypePhaserSym::CalcScaledForwardSeparate(){
	double * emission_probs = new double[num_states];
	double c;

	double scaled_dist;
	double c1,c2;


	printf("Num states = %d \n", num_states);
	printf("Num h = %d \n", num_haps);
	printf("Num markers = %d \n", num_markers);


	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledForward();


	for(int m = 1; m < num_markers; m++){
		//		CalcEmissionProbsIntegrated(m, emission_probs);
		CalcEmissionProbsMla(m, emission_probs);

		//		if(m==3) {
		//			printf("Emission probs[3] h2: \n");
		//			for(int s = 0; s < num_states; s++){
		//				printf("%e \n", emission_probs[s]);
		//			};
		//		};

		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_haps);


		// different
		c2 = scaled_dist / num_haps;
		// same
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

		//#pragma omp parallel for schedule(dynamic,32)
#pragma omp parallel for
		for(int s = 0; s < num_states; s++){

			double sum = 0.0;
			//			int marker_c1 = s / num_h;
			//			int marker_c2 = s % num_h;


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


void HaplotypePhaserSym::CalcScaledForward(){
	double * emission_probs = new double[num_states];
	int * ml_states_other;
	double c;
	double case_probs[3];

	double scaled_dist;
	double c1,c2;

	if (curr_hap == 2) {
		ml_states_other = ml_states_h1;
	}
	if (curr_hap == 1) {
		ml_states_other = ml_states_h2;
	}




	printf("Num states = %d \n", num_states);
	printf("Num h = %d \n", num_haps);
	printf("Num markers = %d \n", num_markers);


	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledForward();


	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbs(m, emission_probs);

		//		if(m==3) {
		//			printf("Emission probs[3] h2: \n");
		//			for(int s = 0; s < num_states; s++){
		//				printf("%e \n", emission_probs[s]);
		//			};
		//		};

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
			//			int marker_c1 = s / num_h;
			//			int marker_c2 = s % num_h;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				int chrom_case = ChromosomePair(ml_states_other[m],s).NumEquals2(ChromosomePair(ml_states_other[m-1],j));

				sum += s_forward[m-1][j] * case_probs[chrom_case];

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




void HaplotypePhaserSym::CalcScaledForwardIntegrated(){
	double * emission_probs = new double[num_states];
	int * ml_states_other;
	double c;
	double case_probs[3];

	double scaled_dist;
	double c1,c2;

	if (curr_hap == 2) {
		ml_states_other = ml_states_h1;
	}
	if (curr_hap == 1) {
		ml_states_other = ml_states_h2;
	}




	printf("Num states = %d \n", num_states);
	printf("Num h = %d \n", num_haps);
	printf("Num markers = %d \n", num_markers);


	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledForward();


	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbs(m, emission_probs);

		//		if(m==3) {
		//			printf("Emission probs[3] h2: \n");
		//			for(int s = 0; s < num_states; s++){
		//				printf("%e \n", emission_probs[s]);
		//			};
		//		};

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
			//			int marker_c1 = s / num_h;
			//			int marker_c2 = s % num_h;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				double isum = 0.0;
				// Integrating over possible states on the other chromosome AT MARKER M
				// Still using ML state for other chromosome at marker m-1

				for(int other_chrom_hap_at_m = 0; other_chrom_hap_at_m < num_states; other_chrom_hap_at_m++){

					//stm1
					int chrom_case = ChromosomePair(other_chrom_hap_at_m,s).NumEquals2(ChromosomePair(ml_states_other[m-1],j));
					isum += case_probs[chrom_case] * all_posteriors(m, other_chrom_hap_at_m);

					//stm2
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


void HaplotypePhaserSym::CalcScaledForwardMarginalized(){
	double * emission_probs = new double[num_states];
	double c, c1, c2;
	double scaled_dist;

	printf("Num states = %d \n", num_states);
	printf("Num h = %d \n", num_haps);
	printf("Num markers = %d \n", num_markers);


	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledForwardMarginalized();
//	printf("Done initi prior \n");

	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbsMarginalized(m, emission_probs);
//		printf("Done Calc Emission probs \n");

		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_haps);

		// different
		c2 = scaled_dist / num_haps;
		// same
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

//		omp_set_num_threads(4);
		// not 4 threads on local computer
		//#pragma omp parallel for schedule(dynamic,32)

		//#pragma omp parallel for
		//#pragma omp parallel for schedule(dynamic)
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
//			cout << "\n Loop Done : Thread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<"\n";
		};


		for(int s = 0; s < num_states; s++){
			c+= s_forward[m][s];
		}

		// c can become 0 there are GLs that are value 0. GLs 0 should not be allowed.
		normalizers[m] = 1.0/c;

		for(int s = 0; s < num_states; s++){
			s_forward[m][s] = s_forward[m][s] * normalizers[m];
		}

		//		cout << "Normalizer: " << normalizers[m] << endl;
		//		cout << "\n Loop Done : Thread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<"\n";

	};

	delete [] emission_probs;
}


void HaplotypePhaserSym::CalcScaledBackwardMarginalized(){
	double * emission_probs = new double[num_states];

	double scaled_dist;
	double c1;
	double c2;


	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledBackward();

	//	for(int s = 0; s < num_states; s++){
	//		cout << "sblast: " << s_backward[num_markers-1][s] << endl;
	//	};


	for(int m = num_markers-2; m >= 0; m--){
		scaled_dist = 1.0 - exp(-(distances[m+1] * pop_const)/num_haps);

		// different
		c2 = scaled_dist / num_haps;
		// same
		c1 = ((1 - scaled_dist)*num_haps + scaled_dist) / num_haps;

		CalcEmissionProbsMarginalized(m+1, emission_probs);

		//		if (m > num_markers -100) {
		//			cout << "eps: ";
		//			for(int s = 0; s < num_states; s++){
		//				cout << " " << emission_probs[s] << " ";
		//			};
		//			cout << endl;
		//
		//		}


#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;

//			printf("emission p: %f \n ", emission_probs[s]);

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


void HaplotypePhaserSym::CalcScaledBackwardSeparate(){
	double * emission_probs = new double[num_states];
	int num_h = 2*num_inds - 2;

	double scaled_dist;
	double c1;
	double c2;


	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledBackward();

	for(int m = num_markers-2; m >= 0; m--){
		scaled_dist = 1.0 - exp(-(distances[m+1] * pop_const)/num_h);

		// different
		c2 = scaled_dist / num_haps;
		// same
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



void HaplotypePhaserSym::CalcScaledBackward(){
	double * emission_probs = new double[num_states];
	int num_h = 2*num_inds - 2;
	int * ml_states_other;

	double case_probs[3];


	if (curr_hap == 2) {
		ml_states_other = ml_states_h1;
	}
	if (curr_hap == 1) {
		ml_states_other = ml_states_h2;
	}



	double scaled_dist;
	double c1;
	double c2;


	double pop_const = (4.0 * Ne) / 100.0;

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

				int chrom_case = ChromosomePair(ml_states_other[m],s).NumEquals2(ChromosomePair(ml_states_other[m+1],j));

				sum +=  s_backward[m+1][j]* case_probs[chrom_case] * emission_probs[j];

			}

			s_backward[m][s] = sum * normalizers[m];

		}
	}
	delete [] emission_probs;
}


void HaplotypePhaserSym::CalcScaledBackwardIntegrated(){
	double * emission_probs = new double[num_states];
	int num_h = 2*num_inds - 2;
	int * ml_states_other;

	double case_probs[3];


	if (curr_hap == 2) {
		ml_states_other = ml_states_h1;
	}
	if (curr_hap == 1) {
		ml_states_other = ml_states_h2;
	}



	double scaled_dist;
	double c1;
	double c2;


	double pop_const = (4.0 * Ne) / 100.0;

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
					int chrom_case = ChromosomePair(other_chrom_hap_at_m,s).NumEquals2(ChromosomePair(ml_states_other[m+1],j));
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





/**
 * For every marker m, get the state with highest posterior probability at that location, given the entire observation sequence
 * (not most likely sequence of states)
 *
 * Calculate array ml_states s.t.
 * ml_states[m] = s_i that maximises P(Q_m = s_i | O_1 ... O_num_markers)
 *
 * for m in {0 ... num_markers-1}
 *
 */
void HaplotypePhaserSym::GetMLHaplotypes(int * ml_states){


#pragma omp parallel for
	for(int m = 0; m < num_markers; m++) {
		double posterior;
		double max_posterior = 0.0;
		int ml_state = -1;
		double norm = 0.0;
#pragma GCC ivdep
		for(int i = 0; i < num_states; i++) {
			norm += s_forward[m][i] * s_backward[m][i];
		}
#pragma GCC ivdep
		for(int s = 0; s < num_states; s++) {
			posterior = s_forward[m][s] * s_backward[m][s] / norm;

			if(posterior > max_posterior) {
				ml_state = s;
				max_posterior = posterior;
			}
		}
		ml_states[m] = ml_state;

	}
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
 * for haplotype 1, marginalized over haplotype 2
 */
vector<vector<double>>  HaplotypePhaserSym::GetPosteriorStatsMarginalized(const char * filename){
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
		//		cout << "m = " << m <<  endl;
		//		for(int s = 0; s < num_states; s++) {
		//			cout << "f: " << s_forward[m][s] << " b: " << s_backward[m][s] << endl;
		//		}

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

			int allele_code = haplotypes(s,m);
			allele_probs[m][allele_code] += posteriors[s];

		}

		vector<size_t> res = sort_indexes(posteriors);

		//		if(m==1 || m==28007 || m==28008 ) {
		//			printf("posteriors at %d: \n", m);
		//			for (int i = 0; i < num_states; i++) {
		//				printf("%d : %f \n", res[i], posteriors[res[i]]);
		//			}
		//		}



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

	}

	//	writeVectorToCSV(filename, stats, "w");
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
 */
vector<vector<double>>  HaplotypePhaserSym::GetPosteriorStats(const char * filename){
	vector<vector<double>> stats;
	vector<vector<double>> geno_probs;
	vector<vector<double>> allele_probs;
	int * ml_alleles_other;
	int * ml_alleles_this;
	int * ml_states_this;
	int * ml_states_other;


	if (curr_hap == 2) {
		ml_alleles_other = ml_alleles_h1;
		ml_alleles_this = ml_alleles_h2;

		ml_states_other = ml_states_h1;
		ml_states_this = ml_states_h2;

	}
	if (curr_hap == 1) {
		ml_alleles_other = ml_alleles_h2;
		ml_alleles_this = ml_alleles_h1;

		ml_states_other = ml_states_h2;
		ml_states_this = ml_states_h1;
	}


	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(44,-1.0);
	}

	for(int m = 0; m < num_markers; m++) {
		allele_probs.push_back({});
		allele_probs[m].resize(2,0.0);
	};

	for(int m = 0; m < num_markers; m++) {
		geno_probs.push_back({});
		geno_probs[m].resize(3,0.0);
	}


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


			sum += posteriors[s];



			//			cout << "Wtf3 " <<  posteriors[s] << " " <<  s_forward[m][s] << " " << s_backward[m][s] << "\n";

			//////////genotype probability - retrieve which genotype this pair of ref haps implies /////////////////

			int ref_hap_this = s;

			// AGCT allele
			//			String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
			//			String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

			// 00, 01, 10, 11


			// use ml allele for genos
			int hapcode_other = ml_alleles_other[m];

			// use ml state for genos
			//			int hapcode_other = haplotypes(ml_states_other[m],m);
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

		float check_sum = 0.0;
		for(int i = 0; i < 3 ; i++) {
			check_sum += geno_probs[m][i];
		}

		//		if(abs(check_sum - 1.0) > 0.000000001 ) {
		//			printf("!!!!!!!!!!!!!!!!!!!!!!!!!Sum of all geno probs is %f at marker %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1\n ", check_sum, m);
		//		}


		vector<size_t> res = sort_indexes(posteriors);

		//		if(m==1 || m==28007 || m==28008 ) {
		//			printf("posteriors at %d: \n", m);
		//			for (int i = 0; i < num_states; i++) {
		//				printf("%d : %f \n", res[i], posteriors[res[i]]);
		//			}
		//		}

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

		if(allele_probs[m][0] >= allele_probs[m][1]) {
			ml_alleles_this[m] = 0;
		}
		else {
			ml_alleles_this[m] = 1;
		}
	}



	for (int i = 0; i < num_markers; i++) {
		ml_states_this[i] = stats[i][39];
	}

	//writeVectorToCSV(filename, stats, "w");


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
 */
vector<vector<double>>  HaplotypePhaserSym::GetPosteriorStatsSpecial(const char * filename){
	vector<vector<double>> stats;
	vector<vector<double>> geno_probs;
	vector<vector<double>> allele_probs;
	int * ml_alleles_other;
	int * ml_alleles_this;
	int * ml_states_this;
	int * ml_states_other;


	if (curr_hap == 2) {
		ml_alleles_other = ml_alleles_h1;
		ml_alleles_this = ml_alleles_h2;

		ml_states_other = ml_states_h1;
		ml_states_this = ml_states_h2;

	}
	if (curr_hap == 1) {
		ml_alleles_other = ml_alleles_h2;
		ml_alleles_this = ml_alleles_h1;

		ml_states_other = ml_states_h2;
		ml_states_this = ml_states_h1;
	}


	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(44,-1.0);
	}

	for(int m = 0; m < num_markers; m++) {
		allele_probs.push_back({});
		allele_probs[m].resize(2,0.0);
	};

	for(int m = 0; m < num_markers; m++) {
		geno_probs.push_back({});
		geno_probs[m].resize(3,0.0);
	}

	//	cout << "1 " <<  "\n";


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

			//			cout << "2 " <<  "\n";

			all_posteriors(m,s) = posteriors[s];

			sum += posteriors[s];

			//			if(m==0 ) {
			//					printf("%d : %f  %f \n", s, s_forward[m][s], s_backward[m][s]);
			//			}



			//			cout << "Wtf3 " <<  posteriors[s] << " " <<  s_forward[m][s] << " " << s_backward[m][s] << "\n";

			//////////genotype probability - retrieve which genotype this pair of ref haps implies /////////////////

			int ref_hap_this = s;

			// AGCT allele
			//			String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
			//			String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

			// 00, 01, 10, 11


			// use ml allele for genos
			//			int hapcode_other = ml_alleles_other[m];

			// use ml state for genos
			int hapcode_other = haplotypes(ml_states_other[m],m);
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

		//		if(m==0 || m==28007 || m==28008 ) {
		//			printf("posteriors at %d: \n", m);
		//			for (int i = 0; i < num_states; i++) {
		//				printf("%d : %f \n", i, posteriors[i]);
		//			}
		//		}


		//		cout << "3 " <<  "\n";

		float check_sum = 0.0;
		for(int i = 0; i < 3 ; i++) {
			check_sum += geno_probs[m][i];
		}

		if(abs(check_sum - 1.0) > 0.0001 ) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!Sum of all geno probs is %f at marker %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1\n ", check_sum, m);
		}


		vector<size_t> res = sort_indexes(posteriors);
		//		cout << "4 " <<  "\n";


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

		//		cout << "5 " <<  "\n";

		if(allele_probs[m][0] >= allele_probs[m][1]) {
			ml_alleles_this[m] = 0;
		}
		else {
			ml_alleles_this[m] = 1;
		}
		if (abs(sum -1 ) > 0.0001){

			cout << "SUM OF POSTERIORS =  " << sum <<  "\n";
		}
	}



	for (int i = 0; i < num_markers; i++) {
		ml_states_this[i] = stats[i][39];
	}

	//	writeVectorToCSV(filename, stats, "w");
	//	cout << "7 " <<  "\n";


	return stats;
}




/**
 * Read stats from filename
 */
vector<vector<double>>  HaplotypePhaserSym::ReadPosteriorStats(const char * filename){

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
 * Translate genotype probabilities to most likely genotypes (represented as haplotype pair).
 *
 * Haplotypes printed to out_file and returned.
 *
 */
HaplotypePair HaplotypePhaserSym::PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file, const char * sample_file){
	std::vector<String> h1;
	std::vector<String> h2;

	//hapcodes 00, 10, 11  represent geno_codes 0,1,2
	int hapcode1;
	int hapcode2;
	int max_geno_code;
	float max_geno_prob;

	VcfFileReader reader;
	VcfHeader header_read;
	reader.open(sample_file, header_read);
	string sample_name = header_read.getSampleName(0);
	reader.close();


	reader.open("template.vcf", header_read);
	VcfRecord record_template;
	record_template.getGenotypeInfo().addStoreField("GT");
	reader.readRecord(record_template);
	reader.close();




	VcfHeader header_new;
	header_new.appendMetaLine("##fileformat=VCFv4.2");
	header_new.appendMetaLine("##FILTER=<ID=PASS,Description=\"All filters passed\">");
	header_new.appendMetaLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	header_new.addHeaderLine(("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	" + sample_name).c_str());

	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_new, InputFile::BGZF);





	for(int m = 0; m < num_markers; m++) {

		max_geno_prob = 0.0;
		for(int i=0; i<3;i++) {
			if(stats[m][41+i] > max_geno_prob) {
				max_geno_prob = stats[m][41+i];
				max_geno_code = i;
			}
		}

		if(max_geno_code == 0) {
			hapcode1 = 0;
			hapcode2 = 0;
		}
		if(max_geno_code == 1) {
			hapcode1 = 0;
			hapcode2 = 1;
		}
		if(max_geno_code == 2) {
			hapcode1 = 1;
			hapcode2 = 1;
		}

		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(hapcode1+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(hapcode2+1));


		/////////////////////

		MarkerInfo* markerinfo = Pedigree::GetMarkerInfo(m);
		std::string marker_name = markerinfo->name.c_str();
		std::size_t delim = marker_name.find(":");
		string chrom =  marker_name.substr(0,delim);
		int pos =  std::stoi(marker_name.substr(delim+1));
		//	printf("%s %d %s \n", chrom.c_str(), pos, (markerinfo->name).c_str());
		record_template.setChrom(chrom.c_str());
		record_template.set1BasedPosition(pos);
		record_template.setID(marker_name.c_str());
		record_template.setRef((markerinfo->GetAlleleLabel(1)).c_str());
		record_template.setAlt((markerinfo->GetAlleleLabel(2)).c_str());
		record_template.setQual(".");

		std::stringstream ss;
		ss << to_string(hapcode1) << "|" << to_string(hapcode2);
		string GTstring = ss.str();

		int succ = record_template.getGenotypeInfo().setString("GT",0, GTstring.c_str());

		if (succ == 0) {
			printf("ERROR IN WRITING TO VCF %s \n", GTstring.c_str());
		}
		else {
			//			printf("Wrote %s \n", GTstring.c_str());
			writer.writeRecord(record_template);
		}


		///////////////////////



	}
	printf("WROTE GENOS TO VCF \n");

	writer.close();
	HaplotypePair hp(h1,h2);
	hp.printToFile(out_file);
	return hp;
	//	Print to file
}


/**
 * Print the given genotypes to a VCF.
 *
 */
void HaplotypePhaserSym::PrintGenotypesToVCF(vector<vector<int>> & genotypes, const char * out_file, const char * sample_file){
//	std::vector<String> h1;
//	std::vector<String> h2;

	//hapcodes 00, 10, 11  represent geno_codes 0,1,2
//	int hapcode1;
//	int hapcode2;
//	int max_geno_code;
//	float max_geno_prob;

	VcfFileReader reader;
	VcfHeader header_read;

	reader.open(sample_file, header_read);
	VcfRecord record_template;

	int num_samples = header_read.getNumSamples();
	reader.readRecord(record_template);

//	record_template.getGenotypeInfo().addStoreField("GT");


//	printf("BEFORE pos %d  num samples: %d\n", record_template.get1BasedPosition(), record_template.getGenotypeInfo().getNumSamples());
	//	reader.close();
	//	record_template.reset();
	//	record_template.getGenotypeInfo().setGT(0,0,0);
	//	printf("AFTER pos %d  num samples: %d\n", record_template.get1BasedPosition(), record_template.getGenotypeInfo().getNumSamples());


//	header_read.appendMeta	Line("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotypez\">");

	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_read, InputFile::BGZF);

	for(int m = 0; m < num_markers; m++) {

		/////////////////////

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
			if (genotypes[sample][m] == 0) {
				succ = record_template.getGenotypeInfo().setString("GL",sample, "1,0,0");
			}
			if (genotypes[sample][m] == 1) {
				succ = record_template.getGenotypeInfo().setString("GL",sample, "0,1,0");
			}
			if (genotypes[sample][m] == 2) {
				succ = record_template.getGenotypeInfo().setString("GL",sample, "0,0,1");
			}


			if (!succ) {
				printf("ERROR IN WRITING TO VCF for marker %d \n", m);
			}
		}

		writer.writeRecord(record_template);

	}

	printf("WROTE GENOS TO %s \n", out_file);

	reader.close();
	writer.close();

}


/**
 * Translate maximum likelihood states to haplotypes.
 *
 * Haplotypes printed to out_file and returned.
 * Also printed to out_file.vcf.gz
 *
 */
HaplotypePair HaplotypePhaserSym::PrintHaplotypesToFile(const char * out_file, const char * sample_file){
	std::vector<String> h1;
	std::vector<String> h2;

	int ref_hap1;
	int ref_hap2;
	//	int prev_ref_hap1;
	//	int prev_ref_hap2;

	VcfFileReader reader;
	VcfHeader header_read;
	reader.open(sample_file, header_read);
	string sample_name = header_read.getSampleName(0);
	reader.close();


	reader.open("template.vcf", header_read);
	VcfRecord record_template;
	record_template.getGenotypeInfo().addStoreField("GT");
	reader.readRecord(record_template);
	reader.close();




	VcfHeader header_new;
	header_new.appendMetaLine("##fileformat=VCFv4.2");
	header_new.appendMetaLine("##FILTER=<ID=PASS,Description=\"All filters passed\">");
	header_new.appendMetaLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	header_new.addHeaderLine(("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	" + sample_name).c_str());

	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_new, InputFile::BGZF);




	// First marker. order does not matter
	ref_hap1 = ml_states_h1[0];
	ref_hap2 = ml_states_h2[0];

	h1.push_back(Pedigree::GetMarkerInfo(0)->GetAlleleLabel(int(haplotypes(ref_hap1,0))+1));
	h2.push_back(Pedigree::GetMarkerInfo(0)->GetAlleleLabel(int(haplotypes(ref_hap2,0))+1));


	MarkerInfo* markerinfo = Pedigree::GetMarkerInfo(0);
	std::string marker_name = markerinfo->name.c_str();
	std::size_t delim = marker_name.find(":");
	string chrom =  marker_name.substr(0,delim);
	int pos =  std::stoi(marker_name.substr(delim+1));
	//	printf("%s %d %s \n", chrom.c_str(), pos, (markerinfo->name).c_str());
	record_template.setChrom(chrom.c_str());
	record_template.set1BasedPosition(pos);
	record_template.setID(marker_name.c_str());
	record_template.setRef((markerinfo->GetAlleleLabel(1)).c_str());
	record_template.setAlt((markerinfo->GetAlleleLabel(2)).c_str());
	record_template.setQual(".");

	std::stringstream ss;
	ss << to_string(int(haplotypes(ref_hap1,0))) << "|" << to_string(int(haplotypes(ref_hap2,0)));
	string GTstring = ss.str();

	int succ = record_template.getGenotypeInfo().setString("GT",0, GTstring.c_str());
	//	int succ = record_template.getGenotypeInfo().setString("GT",0, "1|1");

	//	record_template.getGenotypeInfo().addStoreField("GT");

	if (succ == 0) {
		printf("ERROR IN WRITING TO VCF %s \n", GTstring.c_str());
	}
	else {
		//		printf("Wrote %s \n", GTstring.c_str());
		writer.writeRecord(record_template);
	}


	for(int m = 1; m < num_markers; m++) {


		ref_hap1 = ml_states_h1[m];
		ref_hap2 = ml_states_h2[m];

		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(int(haplotypes(ref_hap1,m))+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(int(haplotypes(ref_hap2,m))+1));

		////////////////////////

		markerinfo = Pedigree::GetMarkerInfo(m);
		marker_name = markerinfo->name.c_str();
		delim = marker_name.find(":");
		chrom =  marker_name.substr(0,delim);
		pos =  std::stoi(marker_name.substr(delim+1));
		//	printf("%s %d %s \n", chrom.c_str(), pos, (markerinfo->name).c_str());
		record_template.setChrom(chrom.c_str());
		record_template.set1BasedPosition(pos);
		record_template.setID(marker_name.c_str());
		record_template.setRef((markerinfo->GetAlleleLabel(1)).c_str());
		record_template.setAlt((markerinfo->GetAlleleLabel(2)).c_str());

		ss.str("");
		ss << to_string(int(haplotypes(ref_hap1,m))) << "|" << to_string(int(haplotypes(ref_hap2,m)));
		GTstring = ss.str();

		int succ = record_template.getGenotypeInfo().setString("GT",0, GTstring.c_str());
		//		int succ = record_template.getGenotypeInfo().setString("GT",0, "0|1");

		if (succ == 0 && m==1) {
			printf("ERROR IN WRITING TO VCF %s %d \n", GTstring.c_str(), haplotypes(ref_hap2,m));

		}
		else {
			writer.writeRecord(record_template);
		}


		///////////////////////


	}
	writer.close();

	HaplotypePair hp(h1,h2);
	hp.printToFile(out_file);
	return hp;
	//	Print to file
}



/**
 * Translate maximum likelihood states to corresponding reference haplotypes.
 * Reference haplotypes at each position printed to stdout and
 * to out_file_ref_haps
 *
 *
 */
void HaplotypePhaserSym::PrintReferenceHaplotypes(int * ml_states, const char * out_file){


	//	std::vector<int> ref1;
	//	std::vector<int> ref2;
	//
	//
	//	int num_h = 2*num_inds - 2;
	//	int ref_hap1;
	//	int ref_hap2;
	//	int prev_ref_hap1;
	//	int prev_ref_hap2;
	//
	//
	//	std::vector<char> codes;
	//
	//	for(int i = 65; i < 65+num_h; i++) {
	//		codes.push_back(i);
	//	}
	//
	//	ref_hap1 = states[ml_states[0]].first;
	//	ref_hap2 = states[ml_states[0]].second;
	//
	//	prev_ref_hap1 = ref_hap1;
	//	prev_ref_hap2 = ref_hap2;
	//	ref1.push_back(ref_hap1);
	//	ref2.push_back(ref_hap2);
	//	for(int m = 1; m < num_markers; m++) {
	//
	//		int first = states[ml_states[m]].first;
	//		int second = states[ml_states[m]].second;
	//
	//		// No switch
	//		if(states[ml_states[m]].NumEquals(states[ml_states[m-1]]) == 2) {
	//			ref_hap1 = prev_ref_hap1;
	//			ref_hap2 = prev_ref_hap2;
	//		}
	//		else {
	//			// One switch
	//			if(states[ml_states[m]].NumEquals(states[ml_states[m-1]]) == 1) {
	//
	//				if(first == prev_ref_hap1) {
	//					ref_hap1 = first;
	//					ref_hap2 = second;
	//				}
	//				else {
	//					if(first == prev_ref_hap2) {
	//						ref_hap1 = second;
	//						ref_hap2 = first;
	//					}
	//					else {
	//						if(second == prev_ref_hap1) {
	//							ref_hap1 = second;
	//							ref_hap2 = first;
	//						}
	//						// here we know s == prev_ref_hap2
	//						else {
	//							ref_hap1 = first;
	//							ref_hap2 = second;
	//						}
	//					}
	//				}
	//			}
	//			// if we have two switches then we can use the original order - should not matter?
	//			else {
	//				ref_hap1 = first;
	//				ref_hap2 = second;
	//			}
	//		}
	//		prev_ref_hap1 = ref_hap1;
	//		prev_ref_hap2 = ref_hap2;
	//
	//		ref1.push_back(ref_hap1);
	//		ref2.push_back(ref_hap2);
	//	}
	//
	//	FILE * hapout = fopen(out_file, "w");
	//
	//	for(auto h : ref1) {
	//		fprintf(hapout,"%d\n",h);
	//	}
	//	putc('\n', hapout);
	//	for(auto h : ref2) {
	//		fprintf(hapout,"%d\n",h);
	//	}
	//
	//	fclose(hapout);

	//	return refs;
	//	Print to file
}



/**
 * Translate maximum likelihood states to corresponding reference haplotypes.
 * Reference haplotypes at each position printed to stdout and
 * to out_file_ref_haps
 *
 *
 */
void HaplotypePhaserSym::PrintReferenceHaplotypesMarginalized(const char * out_file){


	std::vector<int> ref1;

	int ref_hap1;


	std::vector<char> codes;

	for(int i = 65; i < 65+num_haps; i++) {
		codes.push_back(i);
	}

	ref_hap1 = ml_states_h1[0];

	ref1.push_back(ref_hap1);

	for(int m = 1; m < num_markers; m++) {

		ref_hap1 = ml_states_h1[m];

		ref1.push_back(ref_hap1);
	}

	FILE * hapout = fopen(out_file, "w");

	for(auto h : ref1) {
		fprintf(hapout,"%d\n",h);
	}
	putc('\n', hapout);

	fclose(hapout);

}

