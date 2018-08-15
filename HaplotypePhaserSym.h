#ifndef __HAPLOPHASERSYM_H__
#define __HAPLOPHASERSYM_H__

#include "Pedigree.h"
#include "VcfUtils.h"
#include <vector>
#include <chrono>
#include <algorithm>
#include <Eigen/Dense>
#include <math.h>


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

	MatrixXi haplotypes;
	// vector of count of reference allele in reference for every marker
	VectorXi allele_counts;

	//	char ** genotypes;
	vector <double> sample_gls;
	Pedigree ped;
	float theta;
	float error;
	float Ne;

	std::vector<double> distances;

	int distance_code;
	int num_haps;
	vector<ChromosomePair> states;

//	MatrixXi pdf_cases;
//	rowmajdyn pdf_cases_rowmaj;


	// phred_probs[i] contains linear-scale genotype likelihood
	float * phred_probs;

	~HaplotypePhaserSym();
	void LoadData(const String &ref_file, const String &sample_file, int sample_index, const String &map_file);
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
	void CalcEmissionProbsMarginalized(int marker, double * probs);



//	void CalcCases();

	void InitPriorScaledForward();
	void InitPriorScaledBackward();
	void CalcScaledForward();
	void CalcScaledForwardMarginalized();
	void CalcScaledBackward();
	void GetMLHaplotypes(int * ml_states);
	vector<vector<double>> GetPosteriorStats(const char * filename);
	vector<vector<double>>  ReadPosteriorStats(const char * filename);

	HaplotypePair PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file, const char * sample_file);
	HaplotypePair PrintHaplotypesToFile(int * states, const char * out_file,  const char * sample_file);
	void PrintReferenceHaplotypes(int * ml_states, const char * out_file);
};


#endif
