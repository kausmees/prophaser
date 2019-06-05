#ifndef __HAPLOPHASER_H__
#define __HAPLOPHASER_H__

#include "Pedigree.h"
#include "VcfUtils.h"
#include <vector>
#include <chrono>
#include <algorithm>
#include <Eigen/Dense>
#include <math.h>
#include <iomanip>




using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;
using Eigen::RowMajor;
using Eigen::Dynamic;
using Eigen::Matrix;

using namespace Eigen;

typedef Eigen::Matrix<int, Dynamic,Dynamic, RowMajor> rowmajdyn;

/**
 * Represents the reference haplotypes at a loci as an ordered pair.
 *
 */
struct ChromosomePair {
	int first;
	int second;

	ChromosomePair(int a, int b) {
		first = a;
		second = b;

	}

	/**
	 * What of the three transition cases it is to move from this state to the given other state.
	 *
	 * cases:
	 *
	 * 0: both chromosomes switch reference
	 * 1: one switches
	 * 2: no switch
	 */
	int TransitionCase(ChromosomePair other) {

		int num = 0;

		if (first == other.first) {
			num += 1;
		}
		if (second == other.second) {
			num += 1;

		}

		return num;
	}


};

/**

Main haplotype phasing and imputation functionality.

 */
class HaplotypePhaser {

public:

	int prob_precision = 10;

	// the current chromosome that forward-backward is performed along. 0 or 1
	int curr_chrom;
	MatrixXc haplotypes;


	vector <double> sample_gls;
	Pedigree ped;
	float error;
	float Ne;
	double pop_const;
	std::vector<double> distances;

	//	vector<ChromosomePair> states;

	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	int num_haps;


	~HaplotypePhaser();
	void LoadReferenceData(const String &ref_file, String &map_file);
	void LoadSampleData(const String &sample_file,  int sample_index);




	MatrixXd L;
	VectorXd u;
	VectorXd x;

	VectorXd coefs_p0;
	VectorXd coefs_p1;
	double p0;
	double p1;

	//private:

	int num_states;
	int num_markers;

	int num_ref_inds;
	int num_inds;

	double ** s_backward;
	double ** s_forward;
	double * normalizers;


	VectorXi allele_counts;
	int * ml_states_c1;
	int * ml_states_c2;

	int * ml_alleles_c1;
	int * ml_alleles_c2;

	MatrixXd all_posteriors;

	void AllocateMemory();

	void CalcAlleleCounts();
	void CalcEmissionProbs(int marker, double * probs);

	void CalcEmissionProbsMarginalized(int marker, double * probs);
	void CalcEmissionProbsIntegrated(int marker, double * probs);
	void CalcEmissionProbsMla(int marker, double * probs);


	void InitPriorScaledForward();
	void InitPriorScaledForwardMarginalized();

	void InitPriorScaledBackward();

	void CalcScaledForwardFixed();
	void CalcScaledForwardSeparate();
	void CalcScaledForwardIntegrated();
	void CalcScaledForwardMarginalized();

	void CalcMarginalStateProbs(int marker);

	void CalcScaledBackwardFixed();
	void CalcScaledBackwardIntegrated();
	void CalcScaledBackwardMarginalized();
	void CalcScaledBackwardSeparate();

	void GetMLHaplotypes(int * ml_states);
	vector<vector<double>> GetPosteriorStats(const char * filename, bool ml_state, bool print);
	vector<vector<double>> GetPosteriorStatsMarginalized(const char * filename, bool print);

	vector<vector<double>>  ReadPosteriorStats(const char * filename);

	void PrintGenotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);
	void PrintMLSHaplotypesToVCF(vector<vector<int>> & ml_states_h1,vector<vector<int>> & ml_states_h2, const char * out_file, const char * sample_file, const char * vcf_template);
	void PrintMLAHaplotypesToVCF( vector<vector<int>> & ml_alleles_h1, vector<vector<int>> & ml_alleles_h2, const char * out_file, const char * sample_file, const char * vcf_template);




};


#endif
