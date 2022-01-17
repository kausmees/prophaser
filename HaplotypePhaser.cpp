#include <math.h>
#include <string.h>
#include <omp.h>

#include "HaplotypePhaser.h"
#include "MemoryAllocators.h"

#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/device_vector.h>


HaplotypePhaser::~HaplotypePhaser(){
	delete [] normalizersf;
	delete [] normalizersb;
	delete [] emission_probs;
}

struct HapSummer
{
	const int num_haps;
	const int num_states;	
	vector<phaserreal> hapSums[2];

	void reset() {
		for (int k = 0; k < 2; k++) {
			auto& hapSum = hapSums[k];
		#ifdef TARGET_GPU
			hapSum.resize(num_haps);
			#pragma omp target teams distribute parallel for //reduction(+:sums0[:5000], sums1[:5000])			
			for (int i = 0; i < num_haps; i++)
			{
				hapSum[i] = 0;
			}
#else			
			hapSum.assign(num_haps, 0);
#endif		
		}
	}

	HapSummer(int num_haps, int num_states) :
		num_haps(num_haps), num_states(num_states) {
		reset();
	}

	void sum(const phaserreal* __restrict table, const vector<ChromosomePair>& states, const phaserreal* __restrict doEmissions) {		
	#ifdef TARGET_GPU
		int nh2 = num_haps;
		auto key = thrust::make_transform_iterator(thrust::make_counting_iterator(0), [num_haps = nh2] (int j)
		{
			return j / num_haps;
		});

		auto input0 = thrust::make_transform_iterator(thrust::make_counting_iterator(0), [doEmissions, table] (int j)
		{
			return doEmissions ? table[j] * doEmissions[j] : table[j];
		});

		auto input1 = thrust::make_transform_iterator(thrust::make_counting_iterator(0), [doEmissions, table, num_haps = nh2] (int j)
		{
			int major = j / num_haps;
			int minor = j - major * num_haps;
			int j2 = minor * num_haps + major;
			return doEmissions ? table[j2] * doEmissions[j2] : table[j2];
		});

		thrust::reduce_by_key(thrust::device, key, key + num_states, input0, thrust::make_discard_iterator(), hapSums[0].begin());
		thrust::reduce_by_key(thrust::device, key, key + num_states, input1, thrust::make_discard_iterator(), hapSums[1].begin());
#else

		reset();
		auto* sums0 = &hapSums[0][0];
		auto* sums1 = &hapSums[1][0];
#ifdef TARGET_GPU		
		#pragma omp target teams distribute parallel for
#else		
		#pragma omp parallel for schedule(dynamic,131072) reduction(+:sums0[:num_haps], sums1[:num_haps])
#endif
		for (int j = 0; j < num_states; j++) {
			const ChromosomePair cp = states[j];
			const phaserreal val = doEmissions ? table[j] * doEmissions[j] : table[j];
			//phaserreal& val0 = sums0[cp.first];
			//phaserreal& val1 = sums1[cp.second];
#pragma omp atomic
			sums0[cp.first] += val;			
			//atomicAdd(&val0, val);
#pragma omp atomic
			sums1[cp.second] += val;
			//atomicAdd(&val1, val);
}
#endif		

	}
	const array<phaserreal, 3> caseProbs(const phaserreal* __restrict table, const phaserreal* __restrict doEmissions, int s, const ChromosomePair cp) {
		const phaserreal diagonal = doEmissions ? table[s] * doEmissions[s] : table[s];
		const phaserreal halfmatch = hapSums[0][cp.first] + hapSums[1][cp.second] - 2 * diagonal;
		// Nothing special really happens when cp.first == cp.second... RIGHT???
		const phaserreal therest = 1.0 - diagonal - halfmatch;

		return { therest, halfmatch, diagonal };
	}

	const vector<phaserreal>& operator[](int index) {
		return hapSums[index];
	}
};

/*
 * Set the number of haplotypes to use when phasing a sample.
 * The first num_haps haplotypes will be used.
 *
 */
void HaplotypePhaser::SetNumHaps(int new_num_haps){

	num_haps = new_num_haps;
	num_states = pow(num_haps, 2);

	printf("num states: %d\n", num_states);
	emission_probs = new phaserreal[num_states];

	FillStates();

	hapSum = new HapSummer(num_haps, num_states);

	// <decltype(CalcSingleScaledForwardObj)>
	new(&s_forward) StepMemoizer<HaplotypePhaser>(this, num_markers, num_states, true, 64, &HaplotypePhaser::CalcSingleScaledForward);
	new(&s_backward) StepMemoizer<HaplotypePhaser>(this, num_markers, num_states, false, 64, &HaplotypePhaser::CalcSingleScaledBackward);

}

/*
 * Create the states vector of ChromospomePairs corresponding to all pairs of reference haps,
 * when using num_haps haplotypes.
 *
 * if the value of num_haps is changed, this must be called.
 *
 *
 */
void HaplotypePhaser::FillStates(){
	states.clear();

	for(int i = 0; i < num_haps; i++) {
		for(int j = 0; j < num_haps; j++) {
			states.push_back(ChromosomePair(i,j));
		}
	}
}

void HaplotypePhaser::AllocateMemory(){
	String::caseSensitive = false;

	pop_const = (4.0 * Ne) / 100.0;

	distances.resize(num_markers,0.01);


	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	SetNumHaps((num_inds - 1) * 2);


	// haplotypes(h,m) = 0/1 representing allele of haplotype h at marker m
	haplotypes = MatrixXc(num_haps+2, num_markers);

	sample_gls = MatrixXpreal(num_sample_inds, num_markers * 3);

	normalizersf = new phaserreal[num_markers];
	normalizersb = new phaserreal[num_markers];

	for (int i = 0; i < num_markers; i++)
	{
		normalizersf[i] = 0;
		normalizersb[i] = 0;
	}

	// temporary: only consider reference inds ever
	//SetNumHaps(num_ref_inds * 2);

//	hapSum = new HapSummer(num_haps, num_states);

	case_1 = (pow(1 - error, 2) + pow(error, 2));
	case_2 = 2 * (1 - error) * error;
	case_3 = pow(1 - error, 2);
	case_4 = (1 - error) * error;
	case_5 = pow(error,2);

	// <decltype(CalcSingleScaledForwardObj)>
//	new(&s_forward) StepMemoizer<HaplotypePhaser>(this, num_markers, num_states, true, 32, &HaplotypePhaser::CalcSingleScaledForward);
//	new(&s_backward) StepMemoizer<HaplotypePhaser>(this, num_markers, num_states, false, 32, &HaplotypePhaser::CalcSingleScaledBackward);
};


/**
 * Load information about reference and sample individuals, and reference markers
 * into pedigree from the specified files.
 *
 *
 * Load genetic distances from map file if specified.
 *
 */
void HaplotypePhaser::LoadData(const String &ref_file, const String &sample_file, const String &map_file){
	VcfUtils::LoadReferenceMarkers(ref_file);
	VcfUtils::LoadIndividuals(ped,ref_file);
	num_ref_inds = ped.count;
	VcfUtils::LoadIndividuals(ped, sample_file);
	num_sample_inds = ped.count - num_ref_inds;

	num_inds = ped.count;
	num_markers = Pedigree::markerCount;

	printf("Num inds : %d  \nNum markers: %d\n", num_inds, num_markers);

	AllocateMemory();

	VcfUtils::LoadPhasedHaplotypes(ref_file, ped, haplotypes);
	VcfUtils::InitSampleHaplotypes(sample_file, ped, haplotypes, num_ref_inds * 2);
	VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls);

	if(!map_file.IsEmpty()) {
		VcfUtils::LoadGeneticMap(map_file.c_str(), ped, distances);
	};
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
	//VcfUtils::LoadSampleIndividual(ped, sample_file, sample_index);
	//VcfUtils::LoadGenotypeLikelihoods(sample_file, ped, sample_gls, sample_index);

};


/**
 * Translate a state representing a chromosome to an index in the haploytpes matrix.
 *
 * If its a sample that comes after the current sample in haplotypes,
 * we shift the indices 2 steps towards the end of matrix.
 * Because the n_haps first haplotypes in the array are used, but we don't want to
 * use the current sample's own haplotypes as templates.
 *
 */
int HaplotypePhaser::ChromStateToHaplotypeIndex(int chrom_state){

	int newone = chrom_state + 2*(chrom_state >= current_sample_hap_index);

//	if (newone != chrom_state) {
//		printf("!!!!!!!!! chrom_state  %d -> %d  whe ncurrent sample = % d \n", chrom_state, newone, current_sample);
//	}

	return newone;


};


phaserreal HaplotypePhaser::CalcEmissionProb(const ChromosomePair chrom_state, int marker)
{
//		// Reference hapotype at chromosome 1 - fixed (0: REF 1: ALT)
//		int h1 = haplotypes(chrom_state.first, marker);
//
//		// Reference hapotype at chromosome 2 (0: REF 1: ALT)
//		int h2 = haplotypes(chrom_state.second,marker);

		// Reference hapotype at chromosome 1 - fixed (0: REF 1: ALT)
		int h1 = haplotypes(ChromStateToHaplotypeIndex(chrom_state.first), marker);

		// Reference hapotype at chromosome 2 (0: REF 1: ALT)
		int h2 = haplotypes(ChromStateToHaplotypeIndex(chrom_state.second), marker);


		phaserreal sum = 0.0f;
		// case1: g = 0

		if(h1 == 0 and h2 == 0){
			sum += case_3 *  sample_gls(current_sample,marker * 3);
		}
		else {
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){

				sum += case_4 * sample_gls(current_sample,marker * 3);

			}
			else{
				sum +=  case_5 *  sample_gls(current_sample,marker * 3);

			}
		}

		// case2: g = 1
		if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
			sum += case_1 *  sample_gls(current_sample,marker * 3 + 1);

		}
		else{
			sum +=case_2 * sample_gls(current_sample,marker * 3 + 1);

		}

		// case3: g = 2
		if(h1 == 1 and h2 == 1){
			sum += case_3 * sample_gls(current_sample,marker * 3 + 2);

		} else{
			if((h1 == 0 and h2 == 1) or (h1 == 1 and h2 == 0)){
				sum += case_4 * sample_gls(current_sample,marker * 3 + 2);

			}
			else{
				sum += case_5 * sample_gls(current_sample,marker * 3 + 2);

			}
		}

	return sum;
}



/**
 * Get the emission probability of the observed genotype likelihoods at marker
 * for all hidden states
 *
 * P(GLs(marker)|s) for all s in 0...num_states-1
 * for fixed reference haplotype at chromosome 1
 */

void HaplotypePhaser::CalcEmissionProbs(int marker, phaserreal * probs) {

#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for
#else
	#pragma omp parallel for schedule(dynamic,131072)
#endif
	for (int state = 0; state < num_states; state++) {
		probs[state] = CalcEmissionProb(states[state], marker);
	}
}



void HaplotypePhaser::InitPriorScaledForward(){
	phaserreal * emission_probs = new phaserreal[num_states];
	phaserreal prior = 1.0 / num_states;
	phaserreal c1 = 0.0;

	for (int i = 0; i < num_markers; i++)
	{
		normalizersf[i] = 0;
	}


	CalcEmissionProbs(0, emission_probs);

	for(int s = 0; s < num_states; s++){
		c1 += emission_probs[s];
	};

	normalizersf[0] = 1.0/(prior*c1);

	for(int s = 0; s < num_states; s++){
		s_forward[0][s] = (prior*emission_probs[s]) * normalizersf[0];
	}

	delete [] emission_probs;
};

void HaplotypePhaser::InitPriorScaledBackward(){
	phaserreal prior = 1.0 / num_states;
	for (int i = 0; i < num_markers; i++)
	{
		normalizersb[i] = 0;
	}
	normalizersb[num_markers - 1] = 1.0 / prior;
	for(int s = 0; s < num_states; s++){
		s_backward[num_markers-1][s] = normalizersb[num_markers-1];
	}
};

void HaplotypePhaser::CalcSingleScaledForward(int m, const phaserreal* __restrict prev, phaserreal* __restrict now) {
	phaserreal c, c1, c2;
	phaserreal probs[3];
	phaserreal scaled_dist;
	HapSummer& hapSum = *this->hapSum;

	hapSum.sum(prev, states, nullptr);
	c = 0.0;

	scaled_dist = 1 - exp(-(distances[m] * pop_const) / num_haps);

	c1 = scaled_dist / num_haps;
	c2 = 1 - scaled_dist;

	//both_switch
	probs[0] = pow(c1, 2);
	//one switch
	probs[1] = (c2 * c1) + pow(c1, 2);
	// no switch
	probs[2] = pow(c2, 2) + (2 * c2 * c1) + pow(c1, 2);

	// Based on normalization scheme, sum over all previous forwards is always 1
	// Let's precalc sum over all halves in first and second half of pair
	//if (m % 1000 == 0) fprintf(stderr, "%d\n", m);
	phaserreal oldnormalizer = normalizersf[m];

#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for reduction(+ : c)
#else
#pragma omp parallel for schedule(dynamic,131072) reduction(+ : c)
#endif
	for (int s = 0; s < num_states; s++) {
		const ChromosomePair cp = states[s];
		phaserreal sum = 0.0;

		const array<phaserreal, 3> allcases = hapSum.caseProbs(prev, 0, s, cp);
		for (int chrom_case = 0; chrom_case < 3; chrom_case++) {
			sum += allcases[chrom_case] * probs[chrom_case];
		}
		now[s] = CalcEmissionProb(cp, m) * sum;
		if (oldnormalizer) now[s] *= oldnormalizer;
		else
		c += now[s];
	}

	if (!oldnormalizer)
	{
	normalizersf[m] = 1.0 / c;

#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic,131072)
#endif
	for (int s = 0; s < num_states; s++) {
		now[s] = now[s] * normalizersf[m];
	}
	}
}


void HaplotypePhaser::CalcScaledForward(){
	InitPriorScaledForward();
	s_forward.fillAllButFirst();	
}


void HaplotypePhaser::CalcSingleScaledBackward(int m, const phaserreal* __restrict prev, phaserreal* __restrict now) {
	phaserreal probs[3];
	phaserreal scaled_dist;
	phaserreal c, c1, c2;
	vector<ChromosomePair>& states{this->states};
	HapSummer hapSum(num_haps, num_states);
	c = 0;
	phaserreal* emission_probs = this->emission_probs;
	CalcEmissionProbs(m + 1, emission_probs);
	hapSum.sum(prev, states, emission_probs);
	scaled_dist = 1 - exp(-(distances[m + 1] * pop_const) / num_haps);

	// TODO: Shared between both, refactor into helper
	c1 = scaled_dist / num_haps;
	c2 = 1 - scaled_dist;

	//both_switch
	probs[0] = pow(c1, 2);
	//one switch
	probs[1] = (c2 * c1) + pow(c1, 2);
	// no switch
	probs[2] = pow(c2, 2) + (2 * c2 * c1) + pow(c1, 2);
	//if (m % 1000 == 0) fprintf(stderr, "%d\n", m);
	phaserreal oldnormalizer = normalizersb[m];


	if (!oldnormalizer)
	{
	#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for reduction(+ : c)
	#else
	#pragma omp parallel for schedule(dynamic,131072) reduction(+ : c)
	#endif
		for (int s = 0; s < num_states; s++) {
			// Lambda of this not supported in current nvc++, keep in synch with loop below
			const ChromosomePair& cp = states[s];
			phaserreal sum = 0.0f;

			const array<phaserreal, 3> allcases = hapSum.caseProbs(prev, emission_probs, s, cp);

			for (int chrom_case = 0; chrom_case < 3; chrom_case++) {

				sum += allcases[chrom_case] * probs[chrom_case];
			}

			now[s] = sum;
			if (oldnormalizer) now[s] *= oldnormalizer;
			else
			c += now[s];
		}
		normalizersb[m] = 1.0 / c;	

	#ifdef TARGET_GPU
	//#pragma omp target teams distribute parallel for
	#pragma omp target teams num_teams(512) thread_limit(512) distribute parallel for
	#else
	#pragma omp parallel for schedule(dynamic,131072)
	#endif
		for (int s = 0; s < num_states; s++) {
			now[s] = now[s] * normalizersb[m];
		}
	}
	else
	{
		#ifdef TARGET_GPU
		#pragma omp target teams num_teams(512) thread_limit(512) distribute parallel for
		//#pragma omp target loop
		#else
		#pragma omp parallel for schedule(dynamic,131072)
		#endif
		for (int s = 0; s < num_states; s++) {
			// Lambda of this not supported in current nvc++, keep in synch with loop above
			const ChromosomePair& cp = states[s];
			phaserreal sum = 0.0f;

			const array<phaserreal, 3> allcases = hapSum.caseProbs(prev, emission_probs, s, cp);

			for (int chrom_case = 0; chrom_case < 3; chrom_case++) {

				sum += allcases[chrom_case] * probs[chrom_case];
			}

			now[s] = sum;
			if (oldnormalizer) now[s] *= oldnormalizer;
			else
			c += now[s];
		}		
	}
}

void HaplotypePhaser::CalcScaledBackward() {
	InitPriorScaledBackward();
	s_backward.fillAllButFirst();

}





/**
 *
 * Calculate posterior probabilities based on the current values of s_forward and s_backward.
 *
 * For every marker m, get the state with highest posterior probability at that location, given the entire observation sequence
 * (not most likely sequence of states)
 *
 * Calculate array ml_states s.t.
 * ml_states[m] = s_i that maximises P(Q_m = s_i | O_1 ... O_num_markers)
 *
 *
 */
vector<int>  HaplotypePhaser::GetMLStates(){

	vector<int> ml_states(num_markers, -1);

	for(int m = 0; m < num_markers; m++) {
#ifdef TARGET_GPU
		thrust::device_vector<phaserreal> posteriors;
#else
vector<phaserreal> posteriors;
#endif
		posteriors.resize(num_states, -1);

		phaserreal norm = 0.f;
		phaserreal* forward = &s_forward[m][0];
		phaserreal* backward = &s_backward[m][0];
#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for reduction(+ : norm)
#else
#pragma omp parallel for schedule(dynamic,131072) reduction(+ : norm)
#endif
		for(int i = 0; i < num_states; i++) {
			norm += forward[i] * backward[i];
		}

		phaserreal sum = 0.f;
		phaserreal geno0 = 0.f;
		phaserreal geno1 = 0.f;
		phaserreal geno2 = 0.f;
		norm = 1.f / norm;
#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for reduction(+ : sum)
	//
#else
	#pragma omp parallel for schedule(dynamic,131072) reduction(+ : sum)
#endif
		for(int s = 0; s < num_states; s++) {
			phaserreal posterior = forward[s] * backward[s] * norm;

			posteriors[s] = posterior;
			sum += posterior;

		}

#ifdef TARGET_GPU
		thrust::device_vector<phaserreal>::iterator iter = thrust::max_element(posteriors.begin(), posteriors.end());
		unsigned int maxPosteriorIndex = iter - posteriors.begin();
		ml_states[m] = maxPosteriorIndex;

#else
		int maxPosteriorIndex = *max_element(posteriors.rbegin(), posteriors.rend());
		ml_states[m] = maxPosteriorIndex;
#endif

	}
	return ml_states;

}


/**
 * Update the current estimate of a sample's haploytpes according to the most likely states for this sample.
 *
 *
 */
void HaplotypePhaser::UpdateSampleHaplotypes(vector<int> & ml_states) {

	cout << "In UpdateSampleHaplotypes" << endl << endl;
	int template_hap1;
	int template_hap2;
	int hapcode1;
	int hapcode2;


	for (int m = 0;  m < num_markers; m++) {

		if (current_sample == 10) {
			cout << "m = " << m << "ml state = " << ml_states[m] << endl << endl;
			cout << "first =" << states[ml_states[m]].first << "second = " << states[ml_states[m]].second << endl << endl;
		}

		//////////genotype probability/////////////////
		template_hap1 = ChromStateToHaplotypeIndex(states[ml_states[m]].first);
		template_hap2 = ChromStateToHaplotypeIndex(states[ml_states[m]].second);

		if (current_sample == 10) {

		cout << "first translated = " << template_hap1 << "second translated " << template_hap2 << endl << endl;
		}

		// AGCT allele
		// String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
		// String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

		// 00, 01, 10, 11
		hapcode1 = haplotypes(template_hap1, m);
		hapcode2 = haplotypes(template_hap2, m);
		if (current_sample == 10) {

		cout << "hapcode 1 = " << hapcode1 << "hapcode 2" << hapcode2 << endl << endl;
		}
		haplotypes(current_sample_hap_index, m) = hapcode1;
		haplotypes(current_sample_hap_index + 1, m) = hapcode2;

	}
};


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
vector<vector<phaserreal>>  HaplotypePhaser::GetPosteriorStats(const char * filename, bool print){
	vector<vector<phaserreal>> stats;
	vector<vector<phaserreal>> geno_probs;

	for(int m = 0; m < num_markers; m++) {
		stats.push_back({});
		stats[m].resize(44,-1.0);
	}

	for(int m = 0; m < num_markers; m++) {
		geno_probs.push_back({});
		geno_probs[m].resize(3,0.f);
	}


	for(int m = 0; m < num_markers; m++) {
		vector<phaserreal> posteriors;
		posteriors.resize(num_states, -1);

		phaserreal norm = 0.f;
		phaserreal* forward = &s_forward[m][0];
		phaserreal* backward = &s_backward[m][0];
#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for reduction(+ : norm)
#else
#pragma omp parallel for schedule(dynamic,131072) reduction(+ : norm)
#endif
		for(int i = 0; i < num_states; i++) {
			norm += forward[i] * backward[i];
		}

		phaserreal sum = 0.f;
		phaserreal* geno_probs_m = &geno_probs[m][0];
		phaserreal geno0 = 0.f;
		phaserreal geno1 = 0.f;
		phaserreal geno2 = 0.f;
		norm = 1.f / norm;
#ifdef TARGET_GPU
	#pragma omp target teams distribute parallel for reduction(+ : sum) reduction(+ : geno0, geno1, geno2)
	//   
#else
	#pragma omp parallel for schedule(dynamic,131072) reduction(+ : sum) reduction(+ : geno0, geno1, geno2)
#endif
		for(int s = 0; s < num_states; s++) {
			phaserreal posterior = forward[s] * backward[s] * norm;
			// TODO: map back to posteriors[s]?
			sum += posterior;

			//////////genotype probability/////////////////
			int template_hap1 = ChromStateToHaplotypeIndex(states[s].first);
			int template_hap2 = ChromStateToHaplotypeIndex(states[s].second);
				

			// AGCT allele
			// String allele1 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap1][m]+1);
			// String allele2 = Pedigree::GetMarkerInfo(m)->GetAlleleLabel(haplotypes[ref_hap2][m]+1);

			// 00, 01, 10, 11
			int hapcode1 = haplotypes(template_hap1,m);
			int hapcode2 = haplotypes(template_hap2,m);


			int geno_code;

			if(hapcode1 != hapcode2) {
				geno_code = 1;
			}
			else {
				geno_code = (hapcode1 == 0) ? 0 : 2;
			}


			switch (geno_code)
			{
				case 0:
				geno0 += posterior;
				break;
				case 1:
				geno1 += posterior;
				break;
				case 2:
				geno2 += posterior;
				break;
			}
			//geno_probs_m[geno_code] += posteriors[s];

		}
		geno_probs_m[0] = geno0;
		geno_probs_m[1] = geno1;
		geno_probs_m[2] = geno2;

		// check that the sum of genotype probabilities adds up to 1
		float check_sum = 0.f;
		for(int i = 0; i < 3 ; i++) {
			check_sum += geno_probs[m][i];
		}
		if(abs(check_sum - 1.0) > 0.0001 || !isfinite(check_sum)) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!Sum of all geno probs is %f at marker %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1\n ", check_sum, m);
		}


		/*		vector<size_t> res = VcfUtils::sort_indexes(posteriors);

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
			}*/

		stats[m][40] = sum / posteriors.size();
		stats[m][41] = geno_probs[m][0];
		stats[m][42] = geno_probs[m][1];
		stats[m][43] = geno_probs[m][2];

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
vector<vector<phaserreal>>  HaplotypePhaser::ReadPosteriorStats(const char * filename){

	vector<vector<phaserreal>> stats;


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
 * Print the best-guess genotypes and the corresponding posterior genotypes probvabilities to a VCF.
 *
 */
void HaplotypePhaser::PrintPostGenotypesToVCF(vector<vector<int>> & genotypes, vector<vector<vector<phaserreal>>> & postprobs, const char * out_file, const char * sample_file, const char * vcf_template){
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
            string s_postprobs;

            s_postprobs = to_string(postprobs[sample][m][0]) + "," + to_string(postprobs[sample][m][1]) + "," + to_string(postprobs[sample][m][2]);
//            if (m % 1000 == 0) {
//                printf("Writing GP %s for samnple %d at marker %d\n", s_postprobs, sample, m);
//            }
            succ = record_template.getGenotypeInfo().setString("GP", sample, s_postprobs);

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

//	std::vector<String> h1;
//	std::vector<String> h2;

	int template_hap1;
	int template_hap2;


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

				template_hap1 = ChromStateToHaplotypeIndex(states[ml_states[sample][m]].first);
				template_hap2 = ChromStateToHaplotypeIndex(states[ml_states[sample][m]].second);

				std::stringstream ss;
				ss << to_string(int(haplotypes(template_hap1,m))) << "|" << to_string(int(haplotypes(template_hap2,m)));
				string GTstring = ss.str();

				succ = record_template.getGenotypeInfo().setString("GT",sample, GTstring);

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

