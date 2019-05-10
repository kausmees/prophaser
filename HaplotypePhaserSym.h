#ifndef __HAPLOPHASERSYM_H__
#define __HAPLOPHASERSYM_H__

#include "Pedigree.h"
#include "VcfUtils.h"
#include <vector>
#include <chrono>
#include <algorithm>
#include <Eigen/Dense>
#include <math.h>
//#include <sys/stat.h>
//#include <iostream>
//#include <fstream>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;
using Eigen::RowMajor;
using Eigen::Dynamic;
using Eigen::Matrix;

using namespace Eigen;

typedef Eigen::Matrix<int, Dynamic,Dynamic, RowMajor> rowmajdyn;

//using namespace Eigen;
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

		int num = 0;

		if(first == second) {
			if(first == other.first) {
				num += 1;
			}
			if(first == other.second) {
				num += 1;
			}
			return num;
		}

		else {
			if (first == other.first) {
				num += 1;

			}
			else {
				if(first == other.second) {
					num += 1;
				}
			}
			if (second == other.first) {
				num += 1;

			}
			else {
				if(second == other.second) {
					num += 1;
				}
			}
		}

		return num;
	};

	int NumEquals2(ChromosomePair other) {

		int num = 0;

		if(! (first == second)) {

			if (first == other.first) {
				num += 1;

			}
			else {
				num += first == other.second;

			}

			if (second == other.first) {
				num += 1;

			}
			else {
				num += second == other.second;

			}
		}

		else{
			return (first == other.first) + (first == other.second);
		}

		return num;
	}



	//
	// About 2 times longer time than NumEquals2
	int NumEquals3(ChromosomePair other) {


		return (first==second) * ((first==other.first) + (first == other.second)) +
			   (first!=second) * ((first==other.first) + ((first!=other.first)*(first==other.second)) + (second==other.first) + ((second!=other.first)*(second==other.second)));
	};
};

/**

Main haplotype phasing and imputation functionality.

 */
class HaplotypePhaserSym {

public:

//	char ** haplotypes;

	int prob_precision = 10;
	int curr_hap;
	MatrixXc haplotypes;

//	MatrixXi haplotypes;
	// vector of count of alt allele in reference for every marker
	VectorXi allele_counts;

	//	char ** genotypes;
	vector <double> sample_gls;
	Pedigree ped;
	float theta;
	float error;
	float Ne;

	std::vector<double> distances;

	vector<ChromosomePair> states;

	int distance_code;

	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	int num_haps;

	MatrixXd L;
	VectorXd u;
	VectorXd x;

	VectorXd coefs_p0;
	VectorXd coefs_p1;
	double p0;
	double p1;

	MatrixXd all_posteriors;


//	MatrixXi pdf_cases;
//	rowmajdyn pdf_cases_rowmaj;


	// phred_probs[i] contains linear-scale genotype likelihood
	float * phred_probs;

	~HaplotypePhaserSym();
	void LoadData(const String &ref_file, const String &sample_file, int sample_index, const String &map_file);
	void LoadReferenceData(const String &ref_file);
	void LoadSampleData(const String &sample_file, int sample_index);
	void setDistanceCode(int c);




	//private:

	int num_states;
	int num_markers;

	//TODO
	int num_inds;

	double ** s_backward;
	double ** s_forward;
	double * normalizers;

	int * ml_states_h1;
	int * ml_states_h2;

	int * ml_alleles_h1;
	int * ml_alleles_h2;


	//	double ** s_backward2;
	//	double ** s_forward2;
	//	double * normalizers2;


	void AllocateMemory();
	void DeAllocateMemory();

	void CalcAlleleCounts();

	void CalcTransitionProbs(int marker, int marker_state, double * probs);
	void CalcTransitionProbs(int marker, double ** probs);
	void CalcEmissionProbs(int marker, double * probs);
	void CalcMarginalStateProbs(int marker);
	void CalcEmissionProbsMarginalized(int marker, double * probs);
	void CalcEmissionProbsIntegrated(int marker, double * probs);
	void CalcEmissionProbsMla(int marker, double * probs);



//	void CalcCases();

	void InitPriorScaledForward();
	void InitPriorScaledForwardMarginalized();
	void InitPriorScaledBackward();

	void CalcScaledForward();
	void CalcScaledForwardSeparate();
	void CalcScaledForwardIntegrated();



	void CalcScaledForwardMarginalized();

	void CalcScaledBackwardSeparate();
	void CalcScaledBackwardIntegrated();

	void CalcScaledBackward();
	void CalcScaledBackwardMarginalized();
	void GetMLHaplotypes(int * ml_states);
	vector<vector<double>> GetPosteriorStats(const char * filename);
	vector<vector<double>> GetPosteriorStatsSpecial(const char * filename);
	vector<vector<double>> GetPosteriorStatsMarginalized(const char * filename);
	vector<vector<double>>  ReadPosteriorStats(const char * filename);

	HaplotypePair PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file, const char * sample_file);
	void PrintGenotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file);

	HaplotypePair PrintHaplotypesToFile(const char * out_file,  const char * sample_file);

	void PrintReferenceHaplotypes(int * ml_states, const char * out_file);
	void PrintReferenceHaplotypesMarginalized(const char * out_file);

};


#endif
