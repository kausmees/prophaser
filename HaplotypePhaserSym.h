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


	/**
	 *
	 * Number of values that exist in oth pairs (0,1 or 2).
	 *
	 * Represents the transition cases when states are unordered pairs.
	 * Assumed that if the same haplotype exists in both pairs, then there is no switch.
	 *
	 */
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
//	MatrixXi haplotypes;

	int prob_precision = 10;
	int curr_hap;
	MatrixXc haplotypes;

	//	char ** genotypes;
	vector <double> sample_gls;
	Pedigree ped;
	float theta;
	float error;
	float Ne;

	std::vector<double> distances;

	vector<ChromosomePair> states;

	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	int num_haps;

	MatrixXd all_posteriors;


	~HaplotypePhaserSym();
	void LoadReferenceData(const String &ref_file, const char * map_file);
	void LoadSampleData(const String &sample_file,  int sample_index);


	//private:

	int num_states;
	int num_markers;

	int num_ref_inds;
	int num_inds;

	double ** s_backward;
	double ** s_forward;
	double * normalizers;

	void AllocateMemory();

	void CalcEmissionProbs(int marker, double * probs);



	void InitPriorScaledForward();
	void InitPriorScaledBackward();

	void CalcScaledForward();
	void CalcScaledBackward();
	void GetMLHaplotypes(int * ml_states);
	vector<vector<double>> GetPosteriorStats(const char * filename, bool print);
	vector<vector<double>>  ReadPosteriorStats(const char * filename);

//	HaplotypePair PrintGenotypesToFile(vector<vector<double>> & stats, const char * out_file, const char * sample_file);
	void PrintGenotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);

	void PrintHaplotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);

//	void PrintReferenceHaplotypes(int * ml_states, const char * out_file);

};


#endif
