
#include <math.h>
#include <string.h>

#include "HaplotypePhaser.h"
#include "MemoryAllocators.h"


HaplotypePhaser::~HaplotypePhaser(){
	FreeCharMatrix(haplotypes, ped.count*2);
	FreeCharMatrix(genotypes, ped.count);
	delete [] phred_probs;
	FreeFloatMatrix(forward, num_markers);
	FreeFloatMatrix(backward, num_markers);
	FreeFloatMatrix(s_forward, num_markers);
	FreeFloatMatrix(s_backward, num_markers);

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
	genotypes = AllocateCharMatrix(num_inds, num_markers*3);

	forward = AllocateFloatMatrix(num_markers, num_states);
	backward = AllocateFloatMatrix(num_markers, num_states);

	s_forward = AllocateFloatMatrix(num_markers, num_states);
	s_backward = AllocateFloatMatrix(num_markers, num_states);
	normalizers = new float[num_markers];

};

void HaplotypePhaser::DeAllocateMemory(){
	FreeCharMatrix(haplotypes, ped.count*2);
	FreeCharMatrix(genotypes, ped.count);
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

	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, genotypes, sample_index);
//	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
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

	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, genotypes, sample_index);
	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
};



/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * forall hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 *
 */

void HaplotypePhaser::CalcEmissionProbs(int marker, float * probs) {

	for (int state = 0; state < num_states; state++) {

		int num_h = 2*num_inds - 2;
		int h1 = haplotypes[state / num_h][marker];
		int h2 = haplotypes[state % num_h][marker];
		//	printf("h1 = %d  h2 = %d \n", h1, h2);


		float sum = 0.0;

		//TODO maybe use log here
		//	printf(" GetEmissionProb: state = %d marker = %d \n", state, marker);
		//	int ph1 = (unsigned char) genotypes[num_inds-1][marker * 3];
		//	int ph2 = (unsigned char) genotypes[num_inds-1][marker * 3 + 1];
		//	int ph3 = (unsigned char) genotypes[num_inds-1][marker * 3 + 2];
		//
		//	printf("indices : %d %d %d \n", ph1, ph2, ph3);


		// g = 00
		if(h1 == 0 and h2 == 0){
			sum += pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
			//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
		}
		else {
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){

				sum += (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
				//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

			}
			else{
				sum +=  pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
				//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

			}
		}

		// g = 01
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
			sum += (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]];
			//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

		}
		else{
			sum += 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]];
			//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

		}

		// g = 11
		if(h1 == 1 and h2 == 1){
			sum += pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
			//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

		} else{
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum += (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
				//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			}
			else{
				sum += pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
				//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

			}
		}
		probs[state] = sum;
	}

}


/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * for the given hidden state
 *
 * P(GLs(marker)|state)
 */
double HaplotypePhaser::GetEmissionProb(int state, int marker){

	int num_h = 2*num_inds - 2;
	int h1 = haplotypes[state / num_h][marker];
	int h2 = haplotypes[state % num_h][marker];
	//	printf("h1 = %d  h2 = %d \n", h1, h2);


	float sum = 0.0;

	//TODO maybe use log here
	//	printf(" GetEmissionProb: state = %d marker = %d \n", state, marker);
	//	int ph1 = (unsigned char) genotypes[num_inds-1][marker * 3];
	//	int ph2 = (unsigned char) genotypes[num_inds-1][marker * 3 + 1];
	//	int ph3 = (unsigned char) genotypes[num_inds-1][marker * 3 + 2];
	//
	//	printf("indices : %d %d %d \n", ph1, ph2, ph3);


	// g = 00
	if(h1 == 0 and h2 == 0){
		sum += pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
		//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);
	}
	else {
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){

			sum += (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
			//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

		}
		else{
			sum +=  pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]];
			//						printf("adding = %f \n", errors[marker] * 2 * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3]]);

		}
	}

	// g = 01
	if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
		sum += (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]];
		//				printf("adding1 = %f \n", (pow(1 - errors[marker], 2) + pow(errors[marker], 2)) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

	}
	else{
		sum += 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]];
		//				printf("adding2 = %f \n", 2 * (1 - errors[marker]) * errors[marker] * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 1]]);

	}

	// g = 11
	if(h1 == 1 and h2 == 1){
		sum += pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
		//				printf("adding = %f \n", pow(1 - errors[marker], 2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

	} else{
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
			sum += (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
			//						printf("adding = %f \n", (1 - errors[marker]) * error * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

		}
		else{
			sum += pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]];
			//						printf("adding = %f \n", pow(errors[marker],2) * phred_probs[(unsigned char) genotypes[num_inds-1][marker * 3 + 2]]);

		}
	}


	/*
	for(int i = 0; i < 3; i++){
		sum += geno_emission_probs[state*3 + i]*genotypes[num_inds-1][marker*3 + i];
	}
	 */
	return sum;
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
 * TODO use log probs
 *
 */
void HaplotypePhaser::CalcTransitionProbs(int marker, int marker_state, double * probs){
	int num_h = 2*num_inds - 2;

	int marker_c1 = marker_state / num_h;
	int marker_c2 = marker_state % num_h;
	int other_c1;
	int other_c2;

//	double scaled_dist = 1-exp(-distances[marker]/num_h);
//	double scaled_dist = 1-exp(-distances[marker]);
	double scaled_dist = 0.01;


	int occurrences_counter_p1 = 0;
	int occurrences_counter_p2 = 0;
	int occurrences_counter_p3 = 0;
	double prob_sum = 0;
	for(int i = 0; i < num_states; i++){
		other_c1 = i / num_h;
		other_c2 = i % num_h;

		//no switch
		if(other_c1 - marker_c1 == 0 and other_c2 - marker_c2 == 0) {
			probs[i] = pow(1 - scaled_dist, 2) + ((2 * (1 - scaled_dist) * scaled_dist) / num_h) + (pow(scaled_dist,2) / pow(num_h,2));
			prob_sum += probs[i];

			occurrences_counter_p3 ++;
		}
		else {
			// both switch
			if(other_c1-marker_c1 != 0 and other_c2 - marker_c2 != 0) {
				probs[i] = pow(scaled_dist/num_h, 2);
				occurrences_counter_p1 ++;
				prob_sum += probs[i];
			}
			// one switch
			else{
				probs[i] = (((1 - scaled_dist) * scaled_dist) / num_h) + pow(scaled_dist/num_h, 2);
				occurrences_counter_p2 ++;
				prob_sum += probs[i];
			}
		}
	}
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

void HaplotypePhaser::InitPriorForward(){
	double emission_prob;
	float prior = 1.0 / num_states;
	for(int s = 0; s < num_states; s++){
		emission_prob = GetEmissionProb(s,0);
		printf("s = %d emission prob = %f \n", s,emission_prob);
		forward[0][s] = prior*emission_prob;
	};
};

void HaplotypePhaser::InitPriorBackward(){

	for(int s = 0; s < num_states; s++){
		backward[num_markers-1][s] = 1.0;
	};
};


void HaplotypePhaser::CalcForward(){
	double emission_prob;
	// TODO is this an efficient way to handle transition probs??
	double * transition_probs = new double[num_states];
	float sum;
	InitPriorForward();

	for(int m = 1; m < num_markers; m++){
		for(int s = 0; s < num_states; s++){
			emission_prob = GetEmissionProb(s,m);
			printf("marker = %d state = %d emissionprob = %e\n", m,s,emission_prob);
			CalcTransitionProbs(m, s, transition_probs);
			sum = 0.0;
			for(int i = 0; i < num_states; i++){
				sum += forward[m-1][i] * transition_probs[i] * emission_prob;
			}
			forward[m][s] = sum;
		}
	}
	delete [] transition_probs;
}


// TODO do the smoothing at the same time as backward
void HaplotypePhaser::CalcBackward(){
	double emission_prob;

	double * transition_probs = new double[num_states];
	float sum;
	InitPriorBackward();

	for(int m = num_markers-2; m >= 0; m--){
		for(int s = 0; s < num_states; s++){
			CalcTransitionProbs(m+1, s, transition_probs);
			sum = 0.0;
			for(int i = 0; i < num_states; i++){
				emission_prob = GetEmissionProb(i,m+1);
				sum += backward[m+1][i] * transition_probs[i] * emission_prob;
			}
			backward[m][s] = sum;
			printf("backward %d %d = %f", m,s,sum);
		}
	}
	delete [] transition_probs;
}

void HaplotypePhaser::CalcPosteriorOld(){
	float ** posteriors = AllocateFloatMatrix(num_markers, num_states);

	float norm1 = 0;
	for(int i = 0; i < num_states; i++) {
		norm1 += forward[num_markers-1][i];
	}

	float norm2 = 0;
	for(int i = 0; i < num_states; i++) {
		norm2 += backward[0][i];
	}

	float norm3 = 0;
	for(int i = 0; i < num_states; i++) {
		norm3 += forward[0][i] * backward[0][i];
	}
	printf("norm1 = %f\nnorm2 = %f\nnorm3=%f\n", norm1, norm2, norm3);

	printf("POSTERIORS");
	printf("\n");

	for(int m = 0; m < num_markers; m++) {
		for(int s = 0; s < num_states; s++) {
			posteriors[m][s] = forward[m][s] * backward[m][s] / norm1;
		}
	}

	for(int m = 0; m < num_markers; m++) {
		for(int s = 0; s < num_states; s++) {
			printf("%f " ,posteriors[m][s]);
		}
		printf("\n");
	}
	printf("\n");

	FreeFloatMatrix(posteriors, num_markers);
}


void HaplotypePhaser::InitPriorScaledForward(){
	float * emission_probs = new float[num_states];
	float prior = 1.0 / num_states;
	float c1 = 0.0;

	CalcEmissionProbs(0, emission_probs);

	for(int s = 0; s < num_states; s++){
		c1 += emission_probs[s];
	};

	normalizers[0] = prior*c1;

	for(int s = 0; s < num_states; s++){
		s_forward[0][s] = prior*emission_probs[s] / c1;
	};

	delete [] emission_probs;
};

void HaplotypePhaser::InitPriorScaledBackward(){

	for(int s = 0; s < num_states; s++){
		s_backward[num_markers-1][s] = 1.0/normalizers[num_markers-1];
	};
};


void HaplotypePhaser::CalcScaledForward(){
	float min_test = 0.0;
	float max_test = 0.5;


	std::clock_t start;
	double duration;

	float * emission_probs = new float[num_states];
	double * transition_probs = new double[num_states];
	float sum;
	float sum_normalizer;
	float c;
	float c_test;
	float * normalizers_test = new float[num_markers];

	InitPriorScaledForward();

	for(int m = 1; m < num_markers; m++){
		start = std::clock();
//		printf("outer loop: %d of %d \n", m, num_markers);

		CalcEmissionProbs(m, emission_probs);

		c = 0.0;
		c_test = 0.0;
		for(int s = 0; s < num_states; s++){
			CalcTransitionProbs(m, s, transition_probs);
			sum = 0.0;
			sum_normalizer = 0.0;
			for(int j = 0; j < num_states; j++){
				sum += s_forward[m-1][j] * transition_probs[j];
				sum_normalizer += s_forward[m-1][j]*transition_probs[j];
			}
			s_forward[m][s] =  emission_probs[s] * sum;
			c += emission_probs[s]*sum;
			c_test += emission_probs[s]*sum;
		}
		normalizers[m] = c;
		normalizers_test[m] = 1.0/c;


		for(int s = 0; s < num_states; s++){

			float test1 = s_forward[m][s] * normalizers[m];
			float test2 = s_forward[m][s] * normalizers_test[m];

			s_forward[m][s] = s_forward[m][s] / normalizers[m];
			//			s_forward[m][s] = s_forward[m][s] * normalizers_test[m];




//			if(s_forward[m][s] > max_test || s_forward[m][s] < min_test) {
//				printf("OUT OF BOUNDS %d %d \n", m, s);
//			}
//
//			if(m==5 && s==0) {
//				float test3 = s_forward[m][s];
//				printf(" Forward ms = %f \n", s_forward[m][s]);
//				printf("\n");
//			}
//			if(abs(test1 - test2) != 0.0) {
//				printf("Not equal \n");
//			}

		}

		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

//		std::cout<<"duration: "<< duration <<'\n';
	}

	//	for(int m = 0; m < num_markers; m++) {
	//		printf("Normalizer %d %f %f \n", m, normalizers[m], 1/normalizers[m]);
	//	}

	float min = 10000.0;
	float max = 0.0;

	for(int m = 0; m < num_markers; m++) {
		for(int s = 0; s < num_states; s++) {
			//			printf("Normalized alpha %d %f \n", m, s_forward[m][s]);
			if(s_forward[m][s] < min) {
				min = s_forward[m][s];
			}
			if(s_forward[m][s] > max) {
				max = s_forward[m][s];
			}
		}
	}
//	printf("Min = %f Max = %f \n", min, max);

	delete [] transition_probs;
	delete [] emission_probs;
	delete [] normalizers_test;

	//Working, done tests on
	//	float * emission_probs = new float[num_states];
	//	double * transition_probs = new double[num_states];
	//	float sum;
	//	float c;
	//
	//	InitPriorScaledForward();
	//
	//	for(int m = 1; m < num_markers; m++){
	//		start = std::clock();
	//		printf("outer loop: %d of %d \n", m, num_markers);
	//
	//		CalcEmissionProbs(m, emission_probs);
	//
	//		c = 0.0;
	//
	//		for(int k = 0; k < num_states; k++){
	//
	//			CalcTransitionProbs(m, k, transition_probs);
	//			sum = 0.0;
	//			for(int j = 0; j < num_states; j++) {
	//
	//				sum += s_forward[m-1][j]*transition_probs[j];
	//			}
	//
	//			c += emission_probs[k]*sum;
	//		}
	//
	//		normalizers[m] = c;
	//
	//		for(int s = 0; s < num_states; s++){
	//			CalcTransitionProbs(m, s, transition_probs);
	//			sum = 0.0;
	//			for(int j = 0; j < num_states; j++){
	//				sum += s_forward[m-1][j] * transition_probs[j];
	//			}
	//			s_forward[m][s] =  emission_probs[s] * sum / c;
	//		}
	//
	//		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//
	//		std::cout<<"printf: "<< duration <<'\n';
	//	}
	//	delete [] transition_probs;
	//	delete [] emission_probs;

}


void HaplotypePhaser::CalcScaledBackward(){
	float * emission_probs = new float[num_states];
	double * transition_probs = new double[num_states];
	float sum;
	InitPriorScaledBackward();

	for(int m = num_markers-2; m >= 0; m--){

		CalcEmissionProbs(m+1, emission_probs);
		for(int s = 0; s < num_states; s++){

			CalcTransitionProbs(m+1, s, transition_probs);
			sum = 0.0;
			for(int i = 0; i < num_states; i++){

				sum += s_backward[m+1][i] * transition_probs[i] * emission_probs[i];
			}
			s_backward[m][s] = sum / normalizers[m];
		}
	}
	delete [] transition_probs;
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
	float * emission_probs = new float[num_states];
	double * transition_probs = new double[num_states];
	float sum;
	float norm;

	for(int s = 0; s < num_states; s++) {
		if(s_forward[0][s] <= 0.0001 || s_backward[0][s] <= 0.0001) {
			float forw = s_forward[0][s];
			float backw = s_backward[0][s];

		}
	}


	for(int m = 0; m < num_markers; m++) {
		float norm = 0.0;
		for(int i = 0; i < num_states; i++) {
			norm += s_forward[m][i] * s_backward[m][i];
		}
		// TODO use vectors,
		// Only one for alphas, betas and posteriors
		// when calculatig betas, store alpha[m][s] * beta[m][s] in vector[m][s]
		// Use vector sum, accumulate or whatever
		for(int s = 0; s < num_states; s++) {
			if(s== 0 && m==0) {
				float oldOne = s_forward[m][s];
				float newOne = s_forward[m][s] * s_backward[m][s] / norm;
				float backward = s_backward[m][s];

//				printf("\n");
			}
			s_forward[m][s] = s_forward[m][s] * s_backward[m][s] / norm;
		}
	}

	//		for(int m = 0; m < num_markers; m++) {
	//			for(int s = 0; s < num_states; s++) {
	//				printf("%f " ,s_forward[m][s]);
	//			}
	//			printf("\n");
	//		}
	//		printf("\n");


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

	float posterior;
	float max_posterior;
	int ml_state;

	for(int m = 0; m < num_markers; m++) {
		max_posterior = 0.0;
		ml_state = -1;
		float norm = 0.0;

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

	float posterior;
	float max_posterior;
	int ml_state;

	for(int m = 0; m < num_markers; m++) {
		max_posterior = 0.0;
		ml_state = -1;

		for(int s = 0; s < num_states; s++) {
			posterior = s_forward[m][s];

			if(m==0) {
				printf("\n");
			}
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

	for(int m = 0; m < num_markers; m++) {
		ref_hap1 = ml_states[m] / num_h;
		ref_hap2 = ml_states[m] % num_h;

		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));
//		printf(" yy \n");

	}

	HaplotypePair hp(h1,h2);
	hp.printToFile(out_file);
	return hp;
//	Print to file
}
