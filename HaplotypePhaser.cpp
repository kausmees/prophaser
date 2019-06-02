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

	// Number of reference inds.
	num_ref_inds = ped.count;
	// TODO check usage
	num_inds = num_ref_inds +1;

	printf("Num inds tot: %d \n", num_inds);


	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	num_haps = num_ref_inds*2;

	num_states = pow(num_haps,2);

	for(int i = 0; i < num_haps; i++) {
		for(int j = 0; j < num_haps; j++) {
			states.push_back(ChromosomePair(i,j));
		}
	}

	distances.resize(num_markers,0.01);


	// haplotypes(h,m) = 0/1 representing allele of haplotype h at marker m
	haplotypes = MatrixXc(num_haps+2, num_markers);


	sample_gls.resize(num_markers*3);

	s_forward = AllocateDoubleMatrix(num_markers, num_states);
	s_backward = AllocateDoubleMatrix(num_markers, num_states);

	normalizers = new double[num_markers];

};

/**
 * Load reference data from specified file.
 *
 * Load genetic distances from map file.
 *
 */
void HaplotypePhaser::LoadReferenceData(const String &ref_file, const char * map_file){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadReferenceIndividuals(ped,ref_file);
	AllocateMemory();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);
	VcfUtils::LoadGeneticMap(map_file, ped, distances);

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


	for (int state = 0; state < num_states; state++) {

		ChromosomePair chrom_state = states[state];

		// Reference hapotype at chromosome 1 - fixed (0: REF 1: ALT)
		h1 = haplotypes(chrom_state.first, marker);

		// Reference hapotype at chromosome 2 (0: REF 1: ALT)
		h2 = haplotypes(chrom_state.second,marker);

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

void HaplotypePhaser::InitPriorScaledBackward(){

	for(int s = 0; s < num_states; s++){
		s_backward[num_markers-1][s] = normalizers[num_markers-1];
	};
};



void HaplotypePhaser::CalcScaledForward(){
	double * emission_probs = new double[num_states];
	double c;
	double probs[3];
	double scaled_dist;
	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledForward();

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
//			int marker_c1 = s / num_h;
//			int marker_c2 = s % num_h;


			double prob_test = 0.0;
#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){


				int chrom_case = (states[s]).TransitionCase(states[j]);

				sum += s_forward[m-1][j] * probs[chrom_case];
				prob_test +=probs[chrom_case];


//
//
//				int other_c1 = j / num_h;
//				int other_c2 = j % num_h;
//
//				int b = marker_c2-other_c2;
//				int a = marker_c1-other_c1;
//
//				// probs[!a + !b] is p_trans j->s
//				sum += s_forward[m-1][j] * probs[!a+!b];
//
//				prob_test += probs[!a+!b];

				//				int diff = s-j;
				//				int newa = diff / num_h;
				//				int newb = diff-newa*num_h;
				//
				//				if(a!=newa or b!= newb){
				//					int diff = j-s;
				//					int newa = diff / num_h;
				//					int newb = diff-newa*num_h;
				//				}
				//
				//				if(a!=newa or b!= newb){
				//					printf("NOT SAME \n");
				//				}

			}
			s_forward[m][s] =  emission_probs[s] * sum;

			if(abs(prob_test-1.0) > 0.0001) {
				printf("!! Probtest = %f !! \n", prob_test);
			}



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

void HaplotypePhaser::CalcScaledBackward(){
	double * emission_probs = new double[num_states];
	double probs[3];
	double scaled_dist;
	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledBackward();

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
//			int marker_c1 = s / num_haps;
//			int marker_c2 = s % num_haps;

#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

				int chrom_case = (states[s]).TransitionCase(states[j]);

				sum +=  s_backward[m+1][j]* probs[chrom_case] * emission_probs[j];

//
//				int other_c1 = j / num_h;
//				int other_c2 = j % num_h;
//
//				int b = marker_c2-other_c2;
//				int a = marker_c1-other_c1;
//
//				sum += s_backward[m+1][j] * probs[!a + !b] * emission_probs[j];
			}
			s_backward[m][s] = sum * normalizers[m];
		}
	}
	delete [] emission_probs;
}

//
//// ORIGINAL
////void HaplotypePhaser::CalcScaledForward(){
////
////	double * emission_probs = new double[num_states];
////	double ** transition_probs = AllocateDoubleMatrix(num_states, num_states);
////
////	double c;
////	double c2;
////
////
////	InitPriorScaledForward();
////
////	for(int m = 1; m < num_markers; m++){
////		CalcEmissionProbs(m, emission_probs);
////		CalcTransitionProbs(m, transition_probs);
////		c = 0.0;
////		c2 = 0.0;
////#pragma omp parallel for schedule(dynamic,32)
////		for(int s = 0; s < num_states; s++){
////			double sum = 0.0;
////			//			double sum2 = 0.0;
////#pragma GCC ivdep
////			for(int j = 0; j < num_states; j++){
////				sum += s_forward[m-1][j] * transition_probs[s][j];
////				//				sum2 += s_forward2[m-1][j] * transition_probs[s][j];
////
////			}
////			s_forward[m][s] =  emission_probs[s] * sum;
////			//			s_forward2[m][s] =  emission_probs[s] * sum2;
////
////			//			c += emission_probs[s]*sum;
////
////
////		}
////
////		for(int s = 0; s < num_states; s++){
////			c+= s_forward[m][s];
////		}
////
////		//#pragma omp parallel for ordered reduction(+:c2)
////		//		for(int s = 0; s < num_states; s++){
////		//			c2+= s_forward[m][s];
////		//		}
////		//
////		//		if(c != c2) {
////		//			printf("DIFF C %e at %d \n", abs(c-c2), m);
////		//		}
////
////
////		//		normalizers[m] = c;
////		//		normalizers2[m] = 1.0/c2;
////		normalizers[m] = 1.0/c;
////
////
////
////		for(int s = 0; s < num_states; s++){
////			//			s_forward[m][s] = s_forward[m][s] / normalizers[m];
////			//			s_forward2[m][s] = s_forward2[m][s] * normalizers2[m];
////			s_forward[m][s] = s_forward[m][s] * normalizers[m];
////
////
////		}
////	}
////
////	FreeDoubleMatrix(transition_probs, num_states);
////
////	delete [] emission_probs;
////
////}
////
////
////
////
//
//
////// ORIGINAL
////void HaplotypePhaser::CalcScaledBackward(){
////	double * emission_probs = new double[num_states];
////	//	double * transition_probs = new double[num_states];
////	double ** transition_probs = AllocateDoubleMatrix(num_states, num_states);
////
////
////	InitPriorScaledBackward();
////
////
////	for(int m = num_markers-2; m >= 0; m--){
////
////		CalcEmissionProbs(m+1, emission_probs);
////		CalcTransitionProbs(m+1, transition_probs);
////#pragma omp parallel for
////		for(int s = 0; s < num_states; s++){
////			double sum = 0.0;
////			//			double sum2 = 0.0;
////
////#pragma GCC ivdep
////			for(int i = 0; i < num_states; i++){
////
////				sum += s_backward[m+1][i] * transition_probs[s][i] * emission_probs[i];
////				//				sum2 += s_backward2[m+1][i] * transition_probs[s][i] * emission_probs[i];
////
////			}
////			//			s_backward[m][s] = sum / normalizers[m];
////			//			s_backward2[m][s] = sum2 * normalizers2[m];
////			s_backward[m][s] = sum * normalizers[m];
////
////
////		}
////	}
////
////	FreeDoubleMatrix(transition_probs, num_states);
////	delete [] emission_probs;
////}
//
//
//
///**
// * For every marker m, get the state with highest posterior probability at that location, given the entire observation sequence
// * (not most likely sequence of states)
// *
// * Calculate array ml_states s.t.
// * ml_states[m] = s_i that maximises P(Q_m = s_i | O_1 ... O_num_markers)
// *
// * for m in {0 ... num_markers-1}
// *
// */
//void HaplotypePhaser::GetMLHaplotypes(int * ml_states){
//
//
//#pragma omp parallel for
//	for(int m = 0; m < num_markers; m++) {
//		double posterior;
//		double max_posterior = 0.0;
//		int ml_state = -1;
//		double norm = 0.0;
//#pragma GCC ivdep
//		for(int i = 0; i < num_states; i++) {
//			norm += s_forward[m][i] * s_backward[m][i];
//		}
//#pragma GCC ivdep
//		for(int s = 0; s < num_states; s++) {
//			posterior = s_forward[m][s] * s_backward[m][s] / norm;
//
//			if(posterior > max_posterior) {
//				ml_state = s;
//				max_posterior = posterior;
//			}
//		}
//		ml_states[m] = ml_state;
//
//	}
//}



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
 * stats:
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
 */
vector<vector<double>>  HaplotypePhaser::GetPosteriorStats(const char * filename, bool print){
	vector<vector<double>> stats;
	vector<vector<double>> geno_probs;

	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(44,-1.0);
	}

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
			posteriors[s] = s_forward[m][s] * s_backward[m][s] / norm;
			sum += posteriors[s];

			//////////genotype probability/////////////////
			int ref_hap1 = states[s].first;
			int ref_hap2 = states[s].second;


			// AGCT allele
			// String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
			// String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

			// 00, 01, 10, 11
			int hapcode1 = haplotypes(ref_hap1,m);
			int hapcode2 = haplotypes(ref_hap2,m);


			int geno_code;

			if(hapcode1 != hapcode2) {
				geno_code = 1;
			}
			else {
				geno_code = (hapcode1 == 0) ? 0 : 2;
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
		stats[m][41] = geno_probs[m][0];
		stats[m][42] = geno_probs[m][1];
		stats[m][43] = geno_probs[m][2];

	}

	if (print) {
		printf("Writing stats to %s \n", filename);
		writeVectorToCSV(filename, stats, "w");

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

//
///**
// * Translate genotype probabilities to most likely genotypes (represented as haplotype pair).
// *
// * Haplotypes printed to out_file and returned.
// *
// */
//HaplotypePair HaplotypePhaser::PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file){
//	std::vector<String> h1;
//	std::vector<String> h2;
//
//	//hapcodes 00, 10, 11  represent geno_codes 0,1,2
//	int hapcode1;
//	int hapcode2;
//	int max_geno_code;
//	float max_geno_prob;
//
//	for(int m = 0; m < num_markers; m++) {
//
//		max_geno_prob = 0.0;
//		for(int i=0; i<3;i++) {
//			if(stats[m][41+i] > max_geno_prob) {
//				max_geno_prob = stats[m][41+i];
//				max_geno_code = i;
//			}
//		}
//
//		if(max_geno_code == 0) {
//			hapcode1 = 0;
//			hapcode2 = 0;
//		}
//		if(max_geno_code == 1) {
//			hapcode1 = 0;
//			hapcode2 = 1;
//		}
//		if(max_geno_code == 2) {
//			hapcode1 = 1;
//			hapcode2 = 1;
//		}
//
//		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(hapcode1+1));
//		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(hapcode2+1));
//
//	}
//
//	HaplotypePair hp(h1,h2);
//	hp.printToFile(out_file);
//	return hp;
//	//	Print to file
//}
//
//
//
///**
// * Translate maximum likelihood states to haplotypes.
// *
// * Haplotypes printed to out_file and returned.
// *
// */
//HaplotypePair HaplotypePhaser::PrintHaplotypesToFile(int * ml_states, const char * out_file){
//	std::vector<String> h1;
//	std::vector<String> h2;
//
//	int num_h = 2*num_inds - 2;
//	int ref_hap1;
//	int ref_hap2;
//
//	int tester1 =  ml_states[0];
//	int tester2 =  ml_states[1];
//	int tester3 =  ml_states[2];
//
//
//	for(int m = 0; m < num_markers; m++) {
//		ref_hap1 = ml_states[m] / num_h;
//		ref_hap2 = ml_states[m] % num_h;
//
//		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
//		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));
//
//	}
//
//	HaplotypePair hp(h1,h2);
//	hp.printToFile(out_file);
//	return hp;
//	//	Print to file
//}
//
//
///**
// * Translate maximum likelihood states to haplotypes.
// *
// * Haplotypes returned.
// * Reference haplotypes at each position printed to stdout.
// *
// */
//HaplotypePair HaplotypePhaser::PrintReferenceHaplotypes(int * ml_states, const char * out_file){
//
//
//	std::vector<String> h1;
//	std::vector<String> h2;
//
//	std::vector<int> ref1;
//	std::vector<int> ref2;
//
//
//	int num_h = 2*num_inds - 2;
//	int ref_hap1;
//	int ref_hap2;
//
//
//	std::vector<char> codes;
//
//	for(int i = 65; i < 65+num_h; i++) {
//		codes.push_back(i);
//	}
//
//
//	for(int m = 0; m < num_markers; m++) {
//		ref_hap1 = ml_states[m] / num_h;
//		ref_hap2 = ml_states[m] % num_h;
//
//		ref1.push_back(ref_hap1);
//		ref2.push_back(ref_hap2);
//
//
//		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
//		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));
//
//	}
//
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
//
//
//	//	for(auto h : ref1) {
//	//		putc(codes[h], hapout);
//	//	}
//	//	putc('\n', hapout);
//	//	for(auto h : ref2) {
//	//		putc(codes[h], hapout);
//	//	}
//	//
//	//	fclose(hapout);
//
//
//	//	for(auto h : ref1){
//	//		printf("%c",codes[h]);
//	//	}
//	//	printf("\n");
//	//
//	//	for(auto h : ref2){
//	//		printf("%c",codes[h]);
//	//	}
//	//	printf("\n");
//
//	HaplotypePair hp(h1,h2);
//	return hp;
//	//	Print to file
//}




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
//	reader.close();

	printf("Reading from template %s \n", vcf_template);


//	reader.open(sample_file, header_read);


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
 * Print the phased haplotypes indicated by the maximum likelihood states ml_states to vcf.
 *
 * ml_states: n_samples x n_markers vector of most likely states
 * out_file:
 *
 */
void HaplotypePhaser::PrintHaplotypesToVCF(vector<vector<int>> & ml_states, const char * out_file, const char * sample_file, const char * vcf_template){

	std::vector<String> h1;
	std::vector<String> h2;

	int ref_hap1;
	int ref_hap2;
	int prev_ref_hap1;
	int prev_ref_hap2;


	VcfRecord record_template;
	VcfFileReader reader;
	VcfHeader header_read;

	reader.open(vcf_template, header_read);
	reader.readRecord(record_template);
//	reader.close();

	printf("Reading from template %s \n", vcf_template);


//	reader.open(sample_file, header_read);


	int num_samples = header_read.getNumSamples();

	if (num_samples != record_template.getNumSamples()) {
		printf("!!!!! INCONSISTENT NUMBER OF SAMPLES when writing to VCF !!!!!\nSomething is probably wrong with the vcf template \n");
		printf("num samples from header:  %d \n", num_samples);
		printf("num samples in record:  %d \n", record_template.getNumSamples());
	}
	VcfFileWriter writer;
	writer.open((string(out_file) + ".vcf.gz").c_str(), header_read, InputFile::BGZF);

	std::vector<vector<string>> hapstrings;

	for(int sample = 0; sample < num_samples; sample++) {
		hapstrings.push_back({});

		h1.clear();
		h2.clear();

		// First marker. order does not matter
		ref_hap1 = states[ml_states[sample][0]].first;
		ref_hap2 = states[ml_states[sample][0]].second;


		h1.push_back(Pedigree::GetMarkerInfo(0)->GetAlleleLabel(haplotypes(ref_hap1,0)+1));
		h2.push_back(Pedigree::GetMarkerInfo(0)->GetAlleleLabel(haplotypes(ref_hap2,0)+1));

		prev_ref_hap1 = ref_hap1;
		prev_ref_hap2 = ref_hap2;


		std::stringstream ss;
		ss << to_string(int(haplotypes(ref_hap1,0))) << "|" << to_string(int(haplotypes(ref_hap2,0)));
		string GTstring = ss.str();

		hapstrings[sample].push_back(GTstring);


		for(int m = 1; m < num_markers; m++) {

			ref_hap1 = states[ml_states[sample][m]].first;
			ref_hap2 = states[ml_states[sample][m]].second;

//			// No switch
//			if(states[ml_states[sample][m]].NumEquals(states[ml_states[sample][m-1]]) == 2) {
//				ref_hap1 = prev_ref_hap1;
//				ref_hap2 = prev_ref_hap2;
//			}
//			else {
//				// One switch
//				if(states[ml_states[sample][m]].NumEquals(states[ml_states[sample][m-1]]) == 1) {
//
//					if(first == prev_ref_hap1) {
//						ref_hap1 = first;
//						ref_hap2 = second;
//
//					}
//					else {
//						if(first == prev_ref_hap2) {
//							ref_hap1 = second;
//							ref_hap2 = first;
//
//						}
//						else {
//							if(second == prev_ref_hap1) {
//								ref_hap1 = second;
//								ref_hap2 = first;
//
//							}
//							// here we know s == prev_ref_hap2
//							else {
//								ref_hap1 = first;
//								ref_hap2 = second;
//
//							}
//						}
//					}
//				}
//				// if we have two switches then we can use the original order - should not matter?
//				else {
//					ref_hap1 = first;
//					ref_hap2 = second;
//				}
//
//			}
//
//			prev_ref_hap1 = ref_hap1;
//			prev_ref_hap2 = ref_hap2;
//
//
			h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes(ref_hap1,m)+1));
			h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes(ref_hap2,m)+1));

			////////////////////////

			ss.str("");
			ss << to_string(int(haplotypes(ref_hap1,m))) << "|" << to_string(int(haplotypes(ref_hap2,m)));
			GTstring = ss.str();
			hapstrings[sample].push_back(GTstring);

		}
	}


		for(int m = 0; m < num_markers; m++) {

			for(int sample = 0; sample < num_samples; sample++) {

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


				int succ;

				succ = record_template.getGenotypeInfo().setString("GT",sample, hapstrings[sample][m]);

				if (!succ) {
					printf("ERROR IN WRITING TO VCF for marker %d for ind %d \n", m, sample);
				}
			}

			writer.writeRecord(record_template);

		}

		printf("WROTE PHASED TO %s \n", out_file);

		reader.close();
		writer.close();

	}

