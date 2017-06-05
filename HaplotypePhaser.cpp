
#include <math.h>
#include <string.h>
#include <omp.h>

#include "HaplotypePhaser.h"
#include "MemoryAllocators.h"


HaplotypePhaser::~HaplotypePhaser(){
	FreeCharMatrix(haplotypes, ped.count*2);
	//	FreeCharMatrix(genotypes, ped.count);
	delete [] phred_probs;
	FreeFloatMatrix(forward, num_markers);
	FreeFloatMatrix(backward, num_markers);
	FreeDoubleMatrix(s_forward, num_markers);
	FreeDoubleMatrix(s_backward, num_markers);

	//TODO need to delete ped?

}


void HaplotypePhaser::AllocateMemory(){
	String::caseSensitive = false;

	num_states = pow((ped.count-1)*2,2);
	num_markers = Pedigree::markerCount;
	num_inds = ped.count;

	error = 0.01;
	theta = 0.01;
	Ne = 11418;


	distances.resize(num_markers,0.01);
	errors.resize(num_markers,0.01);

	phred_probs = new float[256];

	for(int i = 0; i < 256; i++){
		phred_probs[i] = pow(10.0, -i * 0.1);
	};

	haplotypes = AllocateCharMatrix(num_inds*2, num_markers);
	sample_gls.resize(num_markers*3);
	//	genotypes = AllocateCharMatrix(num_inds, num_markers*3);

	forward = AllocateFloatMatrix(num_markers, num_states);
	backward = AllocateFloatMatrix(num_markers, num_states);

	s_forward = AllocateDoubleMatrix(num_markers, num_states);
	s_backward = AllocateDoubleMatrix(num_markers, num_states);
	normalizers = new double[num_markers];

};

void HaplotypePhaser::DeAllocateMemory(){
	FreeCharMatrix(haplotypes, ped.count*2);
	//	FreeCharMatrix(genotypes, ped.count);
};


/**
 * Load vcf data from files.
 * Reference haplotypes are loaded from ref_file.
 * Unphased sample is loaded from sample_file, the individual with sample_index.
 *
 */
void HaplotypePhaser::LoadData(const String &ref_file, const String &sample_file, int sample_index){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadIndividuals(ped,ref_file, sample_file, sample_index);
	AllocateMemory();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);

	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);
	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
};

/**
 * Load vcf data from files. Only reference data.
 *
 */
void HaplotypePhaser::LoadReferenceData(const String &ref_file, const String &sample_file, int sample_index){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadIndividuals(ped,ref_file, sample_file, sample_index);
	AllocateMemory();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);

	//	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, genotypes, sample_index);
	//	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
};

/**
 * Load vcf data from files. Only reference data.
 *
 */
void HaplotypePhaser::LoadSampleData(const String &ref_file, const String &sample_file, int sample_index){
	//	VcfUtils::LoadReferenceMarkers(ref_file);
	//	VcfUtils::LoadIndividuals(ped,ref_file, sample_file, sample_index);
	//	AllocateMemory();
	//	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);

	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);
	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
};



/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * forall hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 *
 */

void HaplotypePhaser::CalcEmissionProbs(int marker, double * probs) {
	int num_h = 2*num_inds - 2;
	int h1;
	int h2;
	double sum;
	double max = 0.0;
	double min = 10000;


	for (int state = 0; state < num_states; state++) {

		h1 = haplotypes[state / num_h][marker];
		h2 = haplotypes[state % num_h][marker];
		//	printf("h1 = %d  h2 = %d \n", h1, h2);


		sum = 0.0;

		//TODO maybe use log here
		//	printf(" GetEmissionProb: state = %d marker = %d \n", state, marker);
		//	int ph1 = (unsigned char) genotypes[num_inds-1][marker * 3];
		//	int ph2 = (unsigned char) genotypes[num_inds-1][marker * 3 + 1];
		//	int ph3 = (unsigned char) genotypes[num_inds-1][marker * 3 + 2];
		//
		//	printf("indices : %d %d %d \n", ph1, ph2, ph3);


		// case1: g = 0
		if(h1 == 0 and h2 == 0){
			sum += pow(1 - errors[marker], 2) *  sample_gls[marker * 3];
			//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
		}
		else {
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){

				sum += (1 - errors[marker]) * error * sample_gls[marker * 3];
				//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

			}
			else{
				sum +=  pow(errors[marker],2) *  sample_gls[marker * 3];
				//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

			}
		}

		// case2: g = 1
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
			sum += (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) *  sample_gls[marker * 3 + 1];
			//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

		}
		else{
			sum += 2 * (1 - errors[marker]) * errors[marker] * sample_gls[marker * 3 + 1];
			//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

		}

		// case3: g = 2
		if(h1 == 1 and h2 == 1){
			sum += pow(1 - errors[marker], 2) * sample_gls[marker * 3 + 2];
			//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

		} else{
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum += (1 - errors[marker]) * error * sample_gls[marker * 3 + 2];
				//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			}
			else{
				sum += pow(errors[marker],2) * sample_gls[marker * 3 + 2];
				//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			}
		}
		probs[state] = sum;
		//		if(sum > max) {
		//			max = sum;
		//		}
		//		if(sum < min) {
		//			min = sum;
		//		}
	}
	//	printf("Marker: %d min = %f max  = %f \n", marker, min,max);

	//	sum = 0;
	//	for (int state = 0; state < num_states; state++) {
	//		sum += probs[state];
	//
	//	}
	//	printf("Marker: %d Sum over possibilities = %f \n", marker, sum);

}


/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * for the given hidden state
 *
 * P(GLs(marker)|state)
 */
//double HaplotypePhaser::GetEmissionProb(int state, int marker){
//
//	int num_h = 2*num_inds - 2;
//	int h1 = haplotypes[state / num_h][marker];
//	int h2 = haplotypes[state % num_h][marker];
//	//	printf("h1 = %d  h2 = %d \n", h1, h2);
//
//
//	float sum = 0.0;

//TODO maybe use log here
//	printf(" GetEmissionProb: state = %d marker = %d \n", state, marker);
//	int ph1 = (unsigned char) genotypes[num_inds-1][marker * 3];
//	int ph2 = (unsigned char) genotypes[num_inds-1][marker * 3 + 1];
//	int ph3 = (unsigned char) genotypes[num_inds-1][marker * 3 + 2];
//
//	printf("indices : %d %d %d \n", ph1, ph2, ph3);

//
//	// g = 00
//	if(h1 == 0 and h2 == 0){
//		sum += pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
//		//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
//	}
//	else {
//		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
//
//			sum += (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
//			//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
//
//		}
//		else{
//			sum +=  pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
//			//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
//
//		}
//	}
//
//	// g = 01
//	if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
//		sum += (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]];
//		//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);
//
//	}
//	else{
//		sum += 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]];
//		//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);
//
//	}
//
//	// g = 11
//	if(h1 == 1 and h2 == 1){
//		sum += pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
//		//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);
//
//	} else{
//		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
//			sum += (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
//			//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);
//
//		}
//		else{
//			sum += pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
//			//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);
//
//		}
//	}


/*
	for(int i = 0; i < 3; i++){
		sum += geno_emission_probs[state*3 + i]*genotypes[num_inds-1][marker*3 + i];
	}
 */
//	return sum;
//}


/**
 * Calculate transition probabilities.
 * Array probs is given values:
 *
 * P(S_marker = marker_state | S_marker-1 = x) for x in {0...num_states-1} for marker_state in {0...num_states-1}
  *
 * <=>
 *
 * P(S_marker = x | S_marker-1 = marker_state) for x in {0...num_states-1} for marker_state in {0...num_states-1}
 *
 *
 *
 */
void HaplotypePhaser::CalcTransitionProbs(int marker, double ** probs){
	int num_h = 2*num_inds - 2;

#pragma omp parallel for
	for(int marker_state = 0; marker_state < num_states; marker_state++) {

		int marker_c1 = marker_state / num_h;
		int marker_c2 = marker_state % num_h;

		double scaled_dist = 1-exp(-distances[marker]/num_h);
//		double scaled_dist = 1-exp(-distances[marker]);
//		double scaled_dist = 0.01;
#pragma GCC ivdep
		for(int i = 0; i < num_states; i++){
			int other_c1 = i / num_h;
			int other_c2 = i % num_h;

			//no switch
			if(other_c1 - marker_c1 == 0 and other_c2 - marker_c2 == 0) {
				probs[marker_state][i] = pow(1 - scaled_dist, 2) + ((2 * (1 - scaled_dist) * scaled_dist) / num_h) + (pow(scaled_dist,2) / pow(num_h,2));
			}
			else {
				// both switch
				if(other_c1-marker_c1 != 0 and other_c2 - marker_c2 != 0) {
					probs[marker_state][i] = pow(scaled_dist/num_h, 2);

				}
				// one switch
				else{
					probs[marker_state][i] = (((1 - scaled_dist) * scaled_dist) / num_h) + pow(scaled_dist/num_h, 2);

				}
			}
		}
	}
}




/**
 * Calculate transition probabilities.
 * Array probs is given values:
 *
 * P(S_marker = marker_state | S_marker-1 = x) for x in {0...num_states-1}
 *
 * <=>
 *
 * P(S_marker = x | S_marker-1 = marker_state) for x in {0...num_states-1}
 *
 *
 *
 */
void HaplotypePhaser::CalcTransitionProbs(int marker, int marker_state, double * probs){
	int num_h = 2*num_inds - 2;

	int marker_c1 = marker_state / num_h;
	int marker_c2 = marker_state % num_h;
	int other_c1;
	int other_c2;

	double scaled_dist = 1-exp(-distances[marker]/num_h);
	//	double scaled_dist = 1-exp(-distances[marker]);
	//			double scaled_dist = 0.01;


	int occurrences_counter_p1 = 0;
	int occurrences_counter_p2 = 0;
	int occurrences_counter_p3 = 0;
	double prob_sum = 0;

	double max = 0;
	double min = 100;

#pragma omp parallel for
	for(int i = 0; i < num_states; i++){
		other_c1 = i / num_h;
		other_c2 = i % num_h;

		//no switch
		if(other_c1 - marker_c1 == 0 and other_c2 - marker_c2 == 0) {
			probs[i] = pow(1 - scaled_dist, 2) + ((2 * (1 - scaled_dist) * scaled_dist) / num_h) + (pow(scaled_dist,2) / pow(num_h,2));
			prob_sum += probs[i];
			//						printf("case 3 : %e \n", probs[i]);
			//
			//						occurrences_counter_p3 ++;
		}
		else {
			// both switch
			if(other_c1-marker_c1 != 0 and other_c2 - marker_c2 != 0) {
				probs[i] = pow(scaled_dist/num_h, 2);
				//								occurrences_counter_p1 ++;
				prob_sum += probs[i];
				//								printf("case 1 : %e \n", probs[i]);

			}
			// one switch
			else{
				probs[i] = (((1 - scaled_dist) * scaled_dist) / num_h) + pow(scaled_dist/num_h, 2);
				//								occurrences_counter_p2 ++;
				prob_sum += probs[i];
				//								printf("case 2 : %e \n", probs[i]);

			}
		}

		if(probs[i] > max) {
			max = probs[i];
		}
		if(probs[i] < min) {
			min = probs[i];
		}



	}


	//	printf("num states = %d Marker: %d marker state = %d min = %e max = %e probs sum = %f \n", num_states, marker, marker_state, min, max, prob_sum);



	//	for(int i = 0; i < num_states; i++){
	//		other_c1 = i / num_h;
	//		other_c2 = i % num_h;
	//
	//		if(other_c1 - marker_c1 == 0 and other_c2 - marker_c2 == 0) {
	//			probs[i] = pow(1 - distances[marker], 2) + ((2 * (1 - distances[marker]) * distances[marker]) / num_h) + (pow(distances[marker],2) / pow(num_h,2));
	//		}
	//		else {
	//			if(other_c1-marker_c1 != 0 and other_c2 - marker_c2 != 0) {
	//				probs[i] = pow(distances[marker]/num_h, 2);
	//			}
	//			else{
	//				probs[i] = (((1 - distances[marker]) * distances[marker]) / num_h) + pow(distances[marker]/num_h, 2);
	//			}
	//		}
	//	}

}

//void HaplotypePhaser::InitPriorForward(){
//	double emission_prob;
//	float prior = 1.0 / num_states;
//	for(int s = 0; s < num_states; s++){
//		emission_prob = GetEmissionProb(s,0);
//		printf("s = %d emission prob = %f \n", s,emission_prob);
//		forward[0][s] = prior*emission_prob;
//	};
//};
//
//void HaplotypePhaser::InitPriorBackward(){
//
//	for(int s = 0; s < num_states; s++){
//		backward[num_markers-1][s] = 1.0;
//	};
//};
//
//
//void HaplotypePhaser::CalcForward(){
//	double emission_prob;
//	// TODO is this an efficient way to handle transition probs??
//	double * transition_probs = new double[num_states];
//	float sum;
//	InitPriorForward();
//
//	for(int m = 1; m < num_markers; m++){
//		for(int s = 0; s < num_states; s++){
//			emission_prob = GetEmissionProb(s,m);
//			printf("marker = %d state = %d emissionprob = %e\n", m,s,emission_prob);
//			CalcTransitionProbs(m, s, transition_probs);
//			sum = 0.0;
//			for(int i = 0; i < num_states; i++){
//				sum += forward[m-1][i] * transition_probs[i] * emission_prob;
//			}
//			forward[m][s] = sum;
//		}
//	}
//	delete [] transition_probs;
//}
//
//
//// TODO do the smoothing at the same time as backward
//void HaplotypePhaser::CalcBackward(){
//	double emission_prob;
//
//	double * transition_probs = new double[num_states];
//	float sum;
//	InitPriorBackward();
//
//	for(int m = num_markers-2; m >= 0; m--){
//		for(int s = 0; s < num_states; s++){
//			CalcTransitionProbs(m+1, s, transition_probs);
//			sum = 0.0;
//			for(int i = 0; i < num_states; i++){
//				emission_prob = GetEmissionProb(i,m+1);
//				sum += backward[m+1][i] * transition_probs[i] * emission_prob;
//			}
//			backward[m][s] = sum;
//			printf("backward %d %d = %f", m,s,sum);
//		}
//	}
//	delete [] transition_probs;
//}
//
//void HaplotypePhaser::CalcPosteriorOld(){
//	float ** posteriors = AllocateFloatMatrix(num_markers, num_states);
//
//	float norm1 = 0;
//	for(int i = 0; i < num_states; i++) {
//		norm1 += forward[num_markers-1][i];
//	}
//
//	float norm2 = 0;
//	for(int i = 0; i < num_states; i++) {
//		norm2 += backward[0][i];
//	}
//
//	float norm3 = 0;
//	for(int i = 0; i < num_states; i++) {
//		norm3 += forward[0][i] * backward[0][i];
//	}
//	printf("norm1 = %f\nnorm2 = %f\nnorm3=%f\n", norm1, norm2, norm3);
//
//	printf("POSTERIORS");
//	printf("\n");
//
//	for(int m = 0; m < num_markers; m++) {
//		for(int s = 0; s < num_states; s++) {
//			posteriors[m][s] = forward[m][s] * backward[m][s] / norm1;
//		}
//	}
//
//	for(int m = 0; m < num_markers; m++) {
//		for(int s = 0; s < num_states; s++) {
//			printf("%f " ,posteriors[m][s]);
//		}
//		printf("\n");
//	}
//	printf("\n");
//
//	FreeFloatMatrix(posteriors, num_markers);
//}


void HaplotypePhaser::InitPriorScaledForward(){
	double * emission_probs = new double[num_states];
	double prior = 1.0 / num_states;
	double c1 = 0.0;

	CalcEmissionProbs(0, emission_probs);


	for(int s = 0; s < num_states; s++){
		c1 += emission_probs[s];
	};

	normalizers[0] = prior*c1;
	//	printf("First Normalizer = %e Prior = %e c1 = %e \n", normalizers[0], prior, c1);

	for(int s = 0; s < num_states; s++){
		s_forward[0][s] = (prior*emission_probs[s]) / normalizers[0];
	};

	delete [] emission_probs;
};

void HaplotypePhaser::InitPriorScaledBackward(){

	//	printf("Init prior Backward : \n");
	//	for(int i = 0; i< num_markers; i++){
	//		printf(" marker = %d normalizer = %e \n", i, normalizers[i]);
	//	}


	for(int s = 0; s < num_states; s++){
		s_backward[num_markers-1][s] = 1.0 / normalizers[num_markers-1];
	};
};


void HaplotypePhaser::CalcScaledForward(){

	double * emission_probs = new double[num_states];
//	double * transition_probs = new double[num_states];
	double ** transition_probs = AllocateDoubleMatrix(num_states, num_states);

	double sum_normalizer;
	double c;
	double * normalizers_test = new double[num_markers];

	InitPriorScaledForward();

	printf("Num procs %d \n ", omp_get_num_procs());

	for(int m = 1; m < num_markers; m++){
//		printf("Entering parallel section \n");
		CalcEmissionProbs(m, emission_probs);
		CalcTransitionProbs(m, transition_probs);
		c = 0.0;
#pragma omp parallel for schedule(dynamic,32)
		for(int s = 0; s < num_states; s++){
			double sum = 0;

//		    auto cpu = sched_getcpu();
//		    std::ostringstream os;
//		        os<<"\nThread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<std::endl;
//		        std::cout<<os.str()<<std::flush;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){
				sum += s_forward[m-1][j] * transition_probs[s][j];
//				s_forward[m][s] += s_forward[m-1][j] * transition_probs[s][j];

			}
			s_forward[m][s] =  emission_probs[s] * sum;
//			c += emission_probs[s]*sum;

//			s_forward[m][s] =  emission_probs[s] * s_forward[m][s];

//wrong
//#pragma omp atomic
//			c += emission_probs[s]*sum;
//
		}

//		//test : correct
		for(int s = 0; s < num_states; s++){
			c+= s_forward[m][s];
		}

		normalizers[m] = c;

		for(int s = 0; s < num_states; s++){
			s_forward[m][s] = s_forward[m][s] / normalizers[m];

		}
	}

	FreeDoubleMatrix(transition_probs, num_states);
	delete [] emission_probs;
	delete [] normalizers_test;

}


void HaplotypePhaser::CalcScaledBackward(){
	double * emission_probs = new double[num_states];
//	double * transition_probs = new double[num_states];
	double ** transition_probs = AllocateDoubleMatrix(num_states, num_states);


	InitPriorScaledBackward();

	//	printf("In Backward First: \n");
	//	for(int i = 0; i< num_states; i++){
	//		printf(" Backward = %e \n", s_backward[num_markers-1][i]);
	//	}

	for(int m = num_markers-2; m >= 0; m--){

		CalcEmissionProbs(m+1, emission_probs);
		CalcTransitionProbs(m+1, transition_probs);
#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
//			s_backward[m][s] = 0;
			double sum = 0.0;
#pragma GCC ivdep
			for(int i = 0; i < num_states; i++){

				sum += s_backward[m+1][i] * transition_probs[s][i] * emission_probs[i];
//				s_backward[m][s] += s_backward[m+1][i] * transition_probs[s][i] * emission_probs[i];
			}
			s_backward[m][s] = sum / normalizers[m];
//			s_backward[m][s] = s_backward[m][s] / normalizers[m];

		}
	}


	//	printf("In Backward: \n");
	//	for(int i = 0; i< num_states; i++){
	//		printf(" Backward = %e \n", s_backward[1][i]);
	//	}

	FreeDoubleMatrix(transition_probs, num_states);
	delete [] emission_probs;
}


/**
 * Calculate posterior probabilities of states given observations.
 *
 * P(Q_m = s_i | O_1 ... O_num_markers) for all markers m and all states s_i.
 *
 * Using scaled forward and backward values.
 *
 */
void HaplotypePhaser::CalcPosterior(){
	double * emission_probs = new double[num_states];
	double * transition_probs = new double[num_states];
	double sum;
	double norm;

	for(int s = 0; s < num_states; s++) {
		if(s_forward[0][s] <= 0.0001 || s_backward[0][s] <= 0.0001) {
			double forw = s_forward[0][s];
			double backw = s_backward[0][s];

		}
	}

	for(int m = 0; m < num_markers; m++) {
		double norm = 0.0;
		for(int i = 0; i < num_states; i++) {
			norm += s_forward[m][i] * s_backward[m][i];
		}

//#pragma GCC ivdep
		for(int s = 0; s < num_states; s++) {
//			if(m==0) {
//				double oldOne = s_forward[m][s];
//				double newOne = s_forward[m][s] * s_backward[m][s] / norm;
//				double backward = s_backward[m][s];
				//				printf("Forward = %e Backward = %e Norm = %e New = %e \n", oldOne, backward, norm, newOne);
				//
				//				printf("\n");
//			}
			s_forward[m][s] = s_forward[m][s] * s_backward[m][s] / norm;
		}
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
 */
void HaplotypePhaser::SampleHaplotypes(int * ml_states){

	double posterior;
	double max_posterior;
	int ml_state;

	for(int m = 0; m < num_markers; m++) {
		max_posterior = 0.0;
		ml_state = -1;
		double norm = 0.0;

		for(int i = 0; i < num_states; i++) {
			norm += s_forward[m][i] * s_backward[m][i];
		}

		// TODO norm should be 1/norm and multiplication instead of division
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

void HaplotypePhaser::SampleHaplotypesNew(int * ml_states){

	double posterior;
	double max_posterior;
	int ml_state;

	for(int m = 0; m < num_markers; m++) {
		max_posterior = 0.0;
		ml_state = -1;

		if(m==0) {
			for(int s = 0; s < num_states; s++) {
				printf("%e, \n",s_forward[m][s]);

			}
		}



		for(int s = 0; s < num_states; s++) {
			posterior = s_forward[m][s];


			if(posterior > max_posterior) {
				ml_state = s;
				max_posterior = posterior;
			}
		}
		ml_states[m] = ml_state;
	}
}

void HaplotypePhaser::InferHaplotypes() {
	int * ml_states = new int[num_markers];
	char * sampled_haplotypes = new char[2*num_markers];

	int num_h = 2*num_inds - 2;
	int c1;
	int c2;

	SampleHaplotypesNew(ml_states);

	printf("C1\tC2\n");
	for(int m = 0; m < num_markers; m++) {
		c1 = ml_states[m] / num_h;
		c2 = ml_states[m] % num_h;

		//TODO check if implied genotype matches observation
		//edit if not
		sampled_haplotypes[m] = haplotypes[c1][m];
		sampled_haplotypes[num_markers + m] = haplotypes[c2][m];
		printf(" %d\t%d",sampled_haplotypes[m], sampled_haplotypes[num_markers + m]);

	}
	delete [] ml_states;
	delete [] sampled_haplotypes;
}


/**
 * Translate maximum likelihood states to haplotypes.
 *
 * Haplotypes printed to out_file and returned.
 *
 */
HaplotypePair HaplotypePhaser::PrintHaplotypesToFile(int * ml_states, const char * out_file){
	std::vector<String> h1;
	std::vector<String> h2;

	int num_h = 2*num_inds - 2;
	int ref_hap1;
	int ref_hap2;

	int tester1 =  ml_states[0];
	int tester2 =  ml_states[1];
	int tester3 =  ml_states[2];


	for(int m = 0; m < num_markers; m++) {
		ref_hap1 = ml_states[m] / num_h;
		ref_hap2 = ml_states[m] % num_h;

		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));

	}

	HaplotypePair hp(h1,h2);
	hp.printToFile(out_file);
	return hp;
	//	Print to file
}


/**
 * Translate maximum likelihood states to haplotypes.
 *
 * Haplotypes returned.
 * Reference haplotypes at each position printed to stdout.
 *
 */
HaplotypePair HaplotypePhaser::PrintReferenceHaplotypes(int * ml_states){


	std::vector<String> h1;
	std::vector<String> h2;

	std::vector<int> ref1;
	std::vector<int> ref2;


	int num_h = 2*num_inds - 2;
	int ref_hap1;
	int ref_hap2;


	std::vector<char> codes;

	for(int i = 65; i < 65+num_h; i++) {
		codes.push_back(i);
	}



	int tester1 =  ml_states[0];

	for(int m = 0; m < num_markers; m++) {
		ref_hap1 = ml_states[m] / num_h;
		ref_hap2 = ml_states[m] % num_h;

		ref1.push_back(ref_hap1);
		ref2.push_back(ref_hap2);


		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));

	}

	for(auto h : ref1){
		printf("%c",codes[h]);
	}
	printf("\n");

	for(auto h : ref2){
		printf("%c",codes[h]);
	}
	printf("\n");

	HaplotypePair hp(h1,h2);
	return hp;
	//	Print to file
}
