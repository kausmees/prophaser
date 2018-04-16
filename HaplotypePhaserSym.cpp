
#include <math.h>
#include <string.h>
#include <omp.h>

#include "HaplotypePhaserSym.h"
#include "MemoryAllocators.h"




HaplotypePhaserSym::~HaplotypePhaserSym(){
	FreeCharMatrix(haplotypes, ped.count*2);
	//	FreeCharMatrix(genotypes, ped.count);
	delete [] phred_probs;
	FreeDoubleMatrix(s_forward, num_markers);
	FreeDoubleMatrix(s_backward, num_markers);

	//TODO need to delete ped?

}


void HaplotypePhaserSym::AllocateMemory(){
	String::caseSensitive = false;

	num_haps = (ped.count-1)*2;
//	num_haps = 4;


	num_states = ((num_haps)*(num_haps+1)) / 2;



	for(int i = 0; i < num_haps; i++) {
		for(int j = i; j < num_haps; j++) {
			states.push_back(ChromosomePair(i,j));
		}
	}

	num_markers = Pedigree::markerCount;
	num_inds = ped.count;

	//	error = 0.01;
	//	theta = 0.01;

	distances.resize(num_markers,0.01);
	//	errors.resize(num_markers,0.01);

	phred_probs = new float[256];

	for(int i = 0; i < 256; i++){
		phred_probs[i] = pow(10.0, -i * 0.1);
	};

	haplotypes = AllocateCharMatrix(num_inds*2, num_markers);
	sample_gls.resize(num_markers*3);


	s_forward = AllocateDoubleMatrix(num_markers, num_states);
	s_backward = AllocateDoubleMatrix(num_markers, num_states);
	normalizers = new double[num_markers];


	pdf_cases = MatrixXi::Zero(num_states, num_states);
	pdf_cases_rowmaj = rowmajdyn::Zero(num_states,num_states);

	//	s_forward2 = AllocateDoubleMatrix(num_markers, num_states);
	//	s_backward2 = AllocateDoubleMatrix(num_markers, num_states);
	//	normalizers2 = new double[num_markers];

};

void HaplotypePhaserSym::DeAllocateMemory(){
	FreeCharMatrix(haplotypes, ped.count*2);
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
	VcfUtils::LoadIndividuals(ped, ref_file, sample_file, sample_index);
	AllocateMemory();
	CalcCases();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);
	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);
	VcfUtils::LoadGeneticMap(map_file, ped, distances);

};

/**
 * Load vcf data from files. Only reference data.
 *
 */
void HaplotypePhaserSym::LoadReferenceData(const String &ref_file, const String &sample_file, int sample_index){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadIndividuals(ped,ref_file, sample_file, sample_index);
	AllocateMemory();
	VcfUtils::LoadHaplotypes(ref_file, ped, haplotypes);
};

/**
 * Load vcf data from files. Only reference data.
 *
 */
void HaplotypePhaserSym::LoadSampleData(const String &ref_file, const String &sample_file, int sample_index){
	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);
	VcfUtils::LoadGeneticMap("/home/kristiina/Projects/Data/1KGData/maps/chr20.OMNI.interpolated_genetic_map", ped, distances);
	//		VcfUtils::LoadGeneticMap("data/chr20.OMNI.interpolated_genetic_map", ped, distances);

};


void HaplotypePhaserSym::CalcCases(){
	for(int j = 0; j < num_states; j ++) {
		for(int s = 0; s < num_states; s ++) {
			int chrom_case = (states[s]).NumEquals2(states[j]);
			if(states[j].first == states[j].second) {
				chrom_case += 3;
			}
			pdf_cases(j,s) = chrom_case;
		}
	}

	pdf_cases_rowmaj = pdf_cases;

//	std::cout << "Pdf cases: m:\n" << pdf_cases << std::endl;

};


/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * forall hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 *
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



	// OPT1
	for (int state = 0; state < num_states; state++) {
		//	for (ChromosomePair state  : states) {

		ChromosomePair chrom_state = states[state];

		// First reference hapotype at state (0: REF 1: ALT)
		h1 = haplotypes[chrom_state.first][marker];
		// Second reference hapotype at state (0: REF 1: ALT)
		h2 = haplotypes[chrom_state.second][marker];

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



void HaplotypePhaserSym::InitPriorScaledForward(){
	double * emission_probs = new double[num_states];
	double prior = 1.0 / num_states;
	double c1 = 0.0;


	CalcEmissionProbs(0, emission_probs);

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


void HaplotypePhaserSym::CalcScaledForward(){
	double * emission_probs = new double[num_states];
	int num_h = 2 * num_inds - 2;
	double c;
	double case_probs[3];
	double probs[6];

	double scaled_dist;
	double c1,c2;

	printf("Num states = %d \n", num_states);
	printf("Num h = %d \n", num_h);

	double pop_const = (4.0 * Ne) / 100.0;

	InitPriorScaledForward();

	double diff_0 = ((num_h-1) * (num_h-2) / 2);
	double diff_1 = 2*(num_h-1);

	double same_0 = (num_h*(num_h-1)/2);
	double same_1 = (num_h-1);


	for(int m = 1; m < num_markers; m++){
		CalcEmissionProbs(m, emission_probs);
		c = 0.0;

		scaled_dist = 1-exp(-(distances[m] * pop_const)/num_h);

		c1 = scaled_dist/num_h;
		c2 = 1 - scaled_dist;

		//both_switch
		case_probs[0] = pow(c1, 2);
		//one switch
		case_probs[1] =  (c2*c1) + pow(c1, 2);
		// no switch
		case_probs[2] = pow(c2, 2) + (2*c2*c1) + pow(c1, 2);


		double norm_const_case_diff = (diff_0 * case_probs[0]) + (diff_1 * case_probs[1]) + case_probs[2];
		double norm_const_case_same = (same_0 * case_probs[0]) + (same_1 * case_probs[1]) + case_probs[2];

		// diff: (j.first != j.second)
		probs[0] = case_probs[0] / norm_const_case_diff;
		//one switch
		probs[1] =  case_probs[1] / norm_const_case_diff;
		// no switch
		probs[2] = case_probs[2] / norm_const_case_diff;

		// same: (j.first == j.second)
		probs[3] = case_probs[0] / norm_const_case_same;
		//one switch
		probs[4] =  case_probs[1] / norm_const_case_same;
		// no switch
		probs[5] = case_probs[2] / norm_const_case_same;




//		// test that the number of cases
//		// over s is correct
//		int case_counter[3];
//		case_counter[0] = 0;
//		case_counter[1] = 0;
//		case_counter[2] = 0;
//		double prob_test = 0;
//		// the fixed j over which to test the sum over s - this j defines a pdf
//		int testj = 90;

#pragma omp parallel for schedule(dynamic,32)
		for(int s = 0; s < num_states; s++){

			double sum = 0.0;
			//			int marker_c1 = s / num_h;
			//			int marker_c2 = s % num_h;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){


				sum += s_forward[m-1][j] * probs[pdf_cases(j,s)];




//				 prob here is p_trans j -> s
//				 conditioning on j
				// this is the one we expect to sum to 1
//				// over s
//				int chrom_case = pdf_cases(j,s) % 3;
//				if (j==testj) {
//					if(states[j].first == states[j].second) {
//#pragma omp atomic
//						prob_test += probs[pdf_cases(j,s)];
//#pragma omp atomic
//						case_counter[chrom_case] += 1;
//					}
//					else{
//#pragma omp atomic
//
//						prob_test += probs[pdf_cases(j,s)];
//#pragma omp atomic
//						case_counter[chrom_case] += 1;
//					}
//				}




			}
			s_forward[m][s] =  emission_probs[s] * sum;
		}

////		 if case counts are wrong
////		 over s
//		if(states[testj].first == states[testj].second) {
//			if(case_counter[0] != num_h*(num_h-1)/2 || case_counter[1] != (num_h-1) || case_counter[2] != 1) {
//				printf("j= %d : (%d,%d) case counts: 0:%d 1:%d 2:%d \n", testj, states[testj].first, states[testj].second,case_counter[0], case_counter[1], case_counter[2]);
//			}
//		}
//		else{
//			if(case_counter[0] != (num_h-2)*(num_h-1)/2 || case_counter[1] != 2*(num_h-1) || case_counter[2] != 1) {
//				printf("j= %d : (%d,%d) case counts: 0:%d 1:%d 2:%d \n", testj, states[testj].first, states[testj].second,case_counter[0], case_counter[1], case_counter[2]);
//			}
//		}
//		// if prob sum is wrong
//		// over s
//		if(abs(prob_test-1.0) > 0.00001) {
//			printf(" m = % d !! Probtest = %f !! \n", m, prob_test);
//		}





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

void HaplotypePhaserSym::CalcScaledBackward(){
	double * emission_probs = new double[num_states];
	int num_h = 2*num_inds - 2;
	double probs[6];
	double case_probs[3];
	vector<double> probs_diff {0.0, 0.0, 0.0};
	vector<double> probs_same {0.0, 0.0, 0.0};


	double scaled_dist;
	double c1;
	double c2;

	double diff_0 = ((num_h-1) * (num_h-2) / 2);
	double diff_1 = 2*(num_h-1);

	double same_0 = (num_h*(num_h-1)/2);
	double same_1 = (num_h-1);


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


		double norm_const_case_diff = (diff_0 * case_probs[0]) + (diff_1 * case_probs[1]) + case_probs[2];
		double norm_const_case_same = (same_0 * case_probs[0]) + (same_1 * case_probs[1]) + case_probs[2];



		// diff: (j.first != j.second)
		probs[0] = case_probs[0] / norm_const_case_diff;
		//one switch
		probs[1] =  case_probs[1] / norm_const_case_diff;
		// no switch
		probs[2] = case_probs[2] / norm_const_case_diff;

		// same: (j.first == j.second)
		probs[3] = case_probs[0] / norm_const_case_same;
		//one switch
		probs[4] =  case_probs[1] / norm_const_case_same;
		// no switch
		probs[5] = case_probs[2] / norm_const_case_same;

		CalcEmissionProbs(m+1, emission_probs);




#pragma omp parallel for
		for(int s = 0; s < num_states; s++){
			double sum = 0.0;
			//			int marker_c1 = s / num_h;
			//			int marker_c2 = s % num_h;


////			// test that the number of cases
////			// over j is correct
//			int case_counter[3];
//			case_counter[0] = 0;
//			case_counter[1] = 0;
//			case_counter[2] = 0;
//			double prob_test = 0;


#pragma GCC ivdep
			for(int j = 0; j < num_states; j++){

//				chrom_case = (states[s]).NumEquals2(states[j]);
//				// use the pdf that s defines
//				sum += s_backward[m+1][j] * probs_this[chrom_case] * emission_probs[j];

				// works - but inefficiewnt with-column major pdf_cases
				sum +=  s_backward[m+1][j]* probs[pdf_cases_rowmaj(s,j)] * emission_probs[j];


//				// prob here is p_trans s -> j
//				// conditioning on s
//				// this is the one we expect to sum to 1
//				// over j
//				int chrom_case = pdf_cases(s,j) % 3;
//				printf("Chrom case = %d \n", chrom_case);

//					if(states[s].first == states[s].second) {
//#pragma omp atomic
//						prob_test += probs[pdf_cases(s,j)];
//#pragma omp atomic
//						case_counter[chrom_case] += 1;
//					}
//					else{
//#pragma omp atomic
//
//						prob_test += probs[pdf_cases(s,j)];
//#pragma omp atomic
//						case_counter[chrom_case] += 1;
//					}


			}

			s_backward[m][s] = sum * normalizers[m];

//////					 if case counts are wrong
//////					 over s
//			if(states[s].first == states[s].second) {
//				if(case_counter[0] != num_h*(num_h-1)/2 || case_counter[1] != (num_h-1) || case_counter[2] != 1) {
//					printf("bw s= %d : (%d,%d) case counts: 0:%d 1:%d 2:%d \n", s, states[s].first, states[s].second,case_counter[0], case_counter[1], case_counter[2]);
//				}
//			}
//			else{
//				if(case_counter[0] != (num_h-2)*(num_h-1)/2 || case_counter[1] != 2*(num_h-1) || case_counter[2] != 1) {
//					printf("bw s= %d : (%d,%d) case counts: 0:%d 1:%d 2:%d \n", s, states[s].first, states[s].second,case_counter[0], case_counter[1], case_counter[2]);
//				}
//			}
//			// if prob sum is wrong
//			// over s
//			if(abs(prob_test-1.0) > 0.00001) {
//				printf("bw m = % d !! Probtest = %f !! \n", m, prob_test);
//			}



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
 */
vector<vector<double>>  HaplotypePhaserSym::GetPosteriorStats(const char * filename){
	int num_h = 2*num_inds - 2;
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

//			cout << "Wtf3 " <<  posteriors[s] << " " <<  s_forward[m][s] << " " << s_backward[m][s] << "\n";

			//////////genotype probability/////////////////
			int ref_hap1 = states[s].first;
			int ref_hap2 = states[s].second;


			// AGCT allele
			//			String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
			//			String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

			// 00, 01, 10, 11
			int hapcode1 = haplotypes[ref_hap1][m];
			int hapcode2 = haplotypes[ref_hap2][m];

			int geno_code;

			if(hapcode1 != hapcode2) {
				geno_code = 1;
			}
			else {
				geno_code = (hapcode1 == 0) ? 0 : 2;
			}

			geno_probs[m][geno_code] += posteriors[s];

		}

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


	writeVectorToCSV(filename, stats, "w");


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
HaplotypePair HaplotypePhaserSym::PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file){
	std::vector<String> h1;
	std::vector<String> h2;

	//hapcodes 00, 10, 11  represent geno_codes 0,1,2
	int hapcode1;
	int hapcode2;
	int max_geno_code;
	float max_geno_prob;

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

	}

	HaplotypePair hp(h1,h2);
	hp.printToFile(out_file);
	return hp;
	//	Print to file
}



/**
 * Translate maximum likelihood states to haplotypes.
 *
 * Haplotypes printed to out_file and returned.
 * Also printed to out_file.vcf.gz
 *
 */
HaplotypePair HaplotypePhaserSym::PrintHaplotypesToFile(int * ml_states, const char * out_file, const char * sample_file){
	std::vector<String> h1;
	std::vector<String> h2;

	int ref_hap1;
	int ref_hap2;
	int prev_ref_hap1;
	int prev_ref_hap2;

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
	ref_hap1 = states[ml_states[0]].first;
	ref_hap2 = states[ml_states[0]].second;

	h1.push_back(Pedigree::GetMarkerInfo(0)->GetAlleleLabel(haplotypes[ref_hap1][0]+1));
	h2.push_back(Pedigree::GetMarkerInfo(0)->GetAlleleLabel(haplotypes[ref_hap2][0]+1));

	prev_ref_hap1 = ref_hap1;
	prev_ref_hap2 = ref_hap2;


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
	ss << to_string(int(haplotypes[ref_hap1][0])) << "|" << to_string(int(haplotypes[ref_hap2][0]));
	string GTstring = ss.str();

	int succ = record_template.getGenotypeInfo().setString("GT",0, GTstring.c_str());
//	int succ = record_template.getGenotypeInfo().setString("GT",0, "1|1");

	//	record_template.getGenotypeInfo().addStoreField("GT");

	if (succ == 0) {
		printf("ERROR IN WRITING TO VCF %s \n", GTstring.c_str());
	}
	else {
		printf("Wrote %s \n", GTstring.c_str());
		writer.writeRecord(record_template);
	}


	for(int m = 1; m < num_markers; m++) {


		int first = states[ml_states[m]].first;
		int second = states[ml_states[m]].second;

		// No switch
		if(states[ml_states[m]].NumEquals(states[ml_states[m-1]]) == 2) {
			ref_hap1 = prev_ref_hap1;
			ref_hap2 = prev_ref_hap2;
		}
		else {
			// One switch
			if(states[ml_states[m]].NumEquals(states[ml_states[m-1]]) == 1) {

				if(first == prev_ref_hap1) {
					ref_hap1 = first;
					ref_hap2 = second;

				}
				else {
					if(first == prev_ref_hap2) {
						ref_hap1 = second;
						ref_hap2 = first;

					}
					else {
						if(second == prev_ref_hap1) {
							ref_hap1 = second;
							ref_hap2 = first;

						}
						// here we know s == prev_ref_hap2
						else {
							ref_hap1 = first;
							ref_hap2 = second;

						}
					}
				}
			}
			// if we have two switches then we can use the original order - should not matter?
			else {
				ref_hap1 = first;
				ref_hap2 = second;
			}

		}

		prev_ref_hap1 = ref_hap1;
		prev_ref_hap2 = ref_hap2;

		h1.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1));
		h2.push_back(Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1));

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
		ss << to_string(int(haplotypes[ref_hap1][m])) << "|" << to_string(int(haplotypes[ref_hap2][m]));
		GTstring = ss.str();

		int succ = record_template.getGenotypeInfo().setString("GT",0, GTstring.c_str());
//		int succ = record_template.getGenotypeInfo().setString("GT",0, "0|1");

		if (succ == 0 && m==1) {
			printf("ERROR IN WRITING TO VCF %s %d \n", GTstring.c_str(), haplotypes[ref_hap2][m]);

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


	std::vector<int> ref1;
	std::vector<int> ref2;


	int num_h = 2*num_inds - 2;
	int ref_hap1;
	int ref_hap2;
	int prev_ref_hap1;
	int prev_ref_hap2;


	std::vector<char> codes;

	for(int i = 65; i < 65+num_h; i++) {
		codes.push_back(i);
	}

	ref_hap1 = states[ml_states[0]].first;
	ref_hap2 = states[ml_states[0]].second;

	prev_ref_hap1 = ref_hap1;
	prev_ref_hap2 = ref_hap2;
	ref1.push_back(ref_hap1);
	ref2.push_back(ref_hap2);
	for(int m = 1; m < num_markers; m++) {

		int first = states[ml_states[m]].first;
		int second = states[ml_states[m]].second;

		// No switch
		if(states[ml_states[m]].NumEquals(states[ml_states[m-1]]) == 2) {
			ref_hap1 = prev_ref_hap1;
			ref_hap2 = prev_ref_hap2;
		}
		else {
			// One switch
			if(states[ml_states[m]].NumEquals(states[ml_states[m-1]]) == 1) {

				if(first == prev_ref_hap1) {
					ref_hap1 = first;
					ref_hap2 = second;

				}
				else {
					if(first == prev_ref_hap2) {
						ref_hap1 = second;
						ref_hap2 = first;

					}
					else {
						if(second == prev_ref_hap1) {
							ref_hap1 = second;
							ref_hap2 = first;

						}
						// here we know s == prev_ref_hap2
						else {
							ref_hap1 = first;
							ref_hap2 = second;

						}
					}
				}
			}
			// if we have two switches then we can use the original order - should not matter?
			else {
				ref_hap1 = first;
				ref_hap2 = second;
			}

		}
		prev_ref_hap1 = ref_hap1;
		prev_ref_hap2 = ref_hap2;

		ref1.push_back(ref_hap1);
		ref2.push_back(ref_hap2);


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

//	return refs;
	//	Print to file
}


