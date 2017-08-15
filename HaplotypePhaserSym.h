#ifndef __HAPLOPHASERSYM_H__
#define __HAPLOPHASERSYM_H__

#include "Pedigree.h"
#include "VcfUtils.h"
#include <vector>
#include <chrono>
#include <algorithm>


/**
 * Represents the reference haplotypes at a loci as an unordered pair.
 *
 */
struct ChromosomePair {
    int first;
    int second;

    ChromosomePair(int a, int b) {
    	first = std::min(a,b);
    	second = std::max(a,b);
    }

    int NumEquals(ChromosomePair other) {
    	if (first == other.first ) {
    		 if (second == other.second) {
    			 return 2;
    		 }
    		 else {
    			 return 1;
    		 }
    	}
    	else {
    		if(first == other.second) {
    			return 1;
    		}
    		if(second == other.first) {
    			return 1;
    		}
    		return 0;
    	}

    }
};

/**

Main haplotype phasing and imputation functionality.

*/
class HaplotypePhaserSym {

public:
	char ** haplotypes;
//	char ** genotypes;
	vector <double> sample_gls;
	Pedigree ped;
	float theta;
	float error;
	int Ne;

	std::vector<double> distances;

	int distance_code;
	int num_haps;
	vector<ChromosomePair> states;


	// phred_probs[i] contains linear-scale genotype likelihood
	float * phred_probs;

	~HaplotypePhaserSym();
	void LoadData(const String &ref_file, const String &sample_file, int sample_index);
	void LoadReferenceData(const String &ref_file, const String &sample_file, int sample_index);
	void LoadSampleData(const String &ref_file, const String &sample_file, int sample_index);
	void setDistanceCode(int c);




//private:

	int num_states;
	int num_markers;

	//TODO
	int num_inds;

	double ** s_backward;
	double ** s_forward;
	double * normalizers;


//	double ** s_backward2;
//	double ** s_forward2;
//	double * normalizers2;


	void AllocateMemory();
	void DeAllocateMemory();

	void CalcTransitionProbs(int marker, int marker_state, double * probs);
	void CalcTransitionProbs(int marker, double ** probs);
	void CalcEmissionProbs(int marker, double * probs);

	void InitPriorScaledForward();
	void InitPriorScaledBackward();
	void CalcScaledForward();
	void CalcScaledBackward();
	void GetMLHaplotypes(int * ml_states);
	vector<vector<double>> GetPosteriorStats(const char * filename);
	vector<vector<double>>  ReadPosteriorStats(const char * filename);

	HaplotypePair PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file);
	HaplotypePair PrintHaplotypesToFile(int * states, const char * out_file);
	HaplotypePair PrintReferenceHaplotypes(int * ml_states, const char * out_file);
};


#endif
