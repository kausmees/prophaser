
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
void HaplotypePhaser::setDistanceCode(int c) {
	distance_code = c;
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
//		VcfUtils::LoadGeneticMap("data/chr20.OMNI.interpolated_genetic_map", ped, distances);


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
//		VcfUtils::LoadGeneticMap("data/chr20.OMNI.interpolated_genetic_map", ped, distances);

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
	}

}

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

	double scaled_dist;

	if(distance_code == 1) {
		scaled_dist = 1-exp(-distances[marker]);
	}

	if(distance_code == 2) {
		scaled_dist = 1-exp(-distances[marker]/num_h);
	}

	if(distance_code == 3) {
		scaled_dist = 1-exp(-(distances[marker]*4*11000)/num_h);
	}

	if(distance_code == 4) {
		scaled_dist = 0.01;
	}


//		double scaled_dist = 1-exp(-distances[marker]);

//	double scaled_dist = 1-exp(-(distances[marker]*4*11000)/num_h);

//	double
//	double scaled_dist = 0.01;

	if(marker==15) {
		printf("Scaled dist = %e \n", scaled_dist);
	}
#pragma omp parallel for
	for(int marker_state = 0; marker_state < num_states; marker_state++) {

		int marker_c1 = marker_state / num_h;
		int marker_c2 = marker_state % num_h;

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
//void HaplotypePhaser::CalcTransitionProbs(int marker, int marker_state, double * probs){
//	int num_h = 2*num_inds - 2;
//
//	int marker_c1 = marker_state / num_h;
//	int marker_c2 = marker_state % num_h;
//	int other_c1;
//	int other_c2;
//
////	double scaled_dist = 1-exp(-distances[marker]/num_h);
////	double scaled_dist = 1-exp(-distances[marker]);
//	double scaled_dist = 0.01;
//
//	double prob_sum = 0;
//	double max = 0;
//	double min = 100;
//
//#pragma omp parallel for
//	for(int i = 0; i < num_states; i++){
//		other_c1 = i / num_h;
//		other_c2 = i % num_h;
//
//		//no switch
//		if(other_c1 - marker_c1 == 0 and other_c2 - marker_c2 == 0) {
//			probs[i] = pow(1 - scaled_dist, 2) + ((2 * (1 - scaled_dist) * scaled_dist) / num_h) + (pow(scaled_dist,2) / pow(num_h,2));
//			prob_sum += probs[i];
//			//						printf("case 3 : %e \n", probs[i]);
//			//
//			//						occurrences_counter_p3 ++;
//		}
//		else {
//			// both switch
//			if(other_c1-marker_c1 != 0 and other_c2 - marker_c2 != 0) {
//				probs[i] = pow(scaled_dist/num_h, 2);
//				//								occurrences_counter_p1 ++;
//				prob_sum += probs[i];
//				//								printf("case 1 : %e \n", probs[i]);
//
//			}
//			// one switch
//			else{
//				probs[i] = (((1 - scaled_dist) * scaled_dist) / num_h) + pow(scaled_dist/num_h, 2);
//				//								occurrences_counter_p2 ++;
//				prob_sum += probs[i];
//				//								printf("case 2 : %e \n", probs[i]);
//
//			}
//		}
//
////		if(probs[i] > max) {
////			max = probs[i];
////		}
////		if(probs[i] < min) {
////			min = probs[i];
////		}
//
//
//
//	}
//
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
	double ** transition_probs = AllocateDoubleMatrix(num_states, num_states);

	double c;

	InitPriorScaledForward();

	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbs(m, emission_probs);
		CalcTransitionProbs(m, transition_probs);
		c = 0.0;
#pragma omp parallel for schedule(dynamic,32)
		for(int s = 0; s < num_states; s++){
			if(m==1 && s==0) {
				printf("Num threads = %d \n", omp_get_num_threads());
			}
			double sum = 0;
			//			double * transition_probs = new double[num_states];
			//			CalcTransitionProbs(m, s,transition_probs);


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){
				sum += s_forward[m-1][j] * transition_probs[s][j];
				//				sum += s_forward[m-1][j] * transition_probs[j];
			}
			s_forward[m][s] =  emission_probs[s] * sum;
			//			c += emission_probs[s]*sum;

			//			delete [] transition_probs;
		}

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

}


void HaplotypePhaser::CalcScaledBackward(){
	double * emission_probs = new double[num_states];
	//	double * transition_probs = new double[num_states];
	double ** transition_probs = AllocateDoubleMatrix(num_states, num_states);


	InitPriorScaledBackward();


	for(int m = num_markers-2; m >= 0; m--){

		CalcEmissionProbs(m+1, emission_probs);
		CalcTransitionProbs(m+1, transition_probs);
#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;
#pragma GCC ivdep
			for(int i = 0; i < num_states; i++){

				sum += s_backward[m+1][i] * transition_probs[s][i] * emission_probs[i];
			}
			s_backward[m][s] = sum / normalizers[m];
		}
	}

	FreeDoubleMatrix(transition_probs, num_states);
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
void HaplotypePhaser::GetMLHaplotypes(int * ml_states){


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
 */
vector<vector<double>>  HaplotypePhaser::GetPosteriorStats(const char * filename){

	vector<vector<double>> stats;

	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(41,-1.0);
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
			posteriors[s] = s_forward[m][s] * s_backward[m][s] / norm;
			sum += posteriors[s];
		}

		vector<size_t> res = sort_indexes(posteriors);

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

	}
	writeVectorToCSV(filename, stats, "w");
	return stats;
}

/**
 * Read stats from filena,e
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
HaplotypePair HaplotypePhaser::PrintReferenceHaplotypes(int * ml_states, const char * out_file){


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


	for(int m = 0; m < num_markers; m++) {
		ref_hap1 = ml_states[m] / num_h;
		ref_hap2 = ml_states[m] % num_h;

		ref1.push_back(ref_hap1);
		ref2.push_back(ref_hap2);


		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));

	}


	FILE * hapout = fopen(out_file, "w");

	for(auto h : ref1) {
		fprintf(hapout,"%d\n",h);
	}
	putc('\n', hapout);
	for(auto h : ref2) {
		fprintf(hapout,"%d\n",h);
	}

	fclose(hapout);


	//	for(auto h : ref1) {
	//		putc(codes[h], hapout);
	//	}
	//	putc('\n', hapout);
	//	for(auto h : ref2) {
	//		putc(codes[h], hapout);
	//	}
	//
	//	fclose(hapout);


//	for(auto h : ref1){
//		printf("%c",codes[h]);
//	}
//	printf("\n");
//
//	for(auto h : ref2){
//		printf("%c",codes[h]);
//	}
//	printf("\n");

	HaplotypePair hp(h1,h2);
	return hp;
	//	Print to file
}


