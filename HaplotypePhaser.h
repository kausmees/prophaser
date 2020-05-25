#ifndef __HAPLOPHASER_H__
#define __HAPLOPHASER_H__

#include "Pedigree.h"
#include "VcfUtils.h"
#include <vector>
#include <chrono>
#include <algorithm>

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
//using Eigen::Matrix;

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

template<class T>
class StepMemoizer
{
	vector<double*> mainTable;
	vector<double*> auxTable;
	T* parent;
	int num_markers;
	bool dirup;
	int step;
	int auxAnchor = -1;	
	int totalAnchor;
	using gent = void (T::*)(int, const double*, double*);
	gent generator;
	

public:
	StepMemoizer() = default;

	// TODO: Copy assignment is dangerous.
	StepMemoizer(T* parent, int num_markers, int num_states, bool dirup, int step, gent generator) : parent(parent), num_markers(num_markers), dirup(dirup), step(step), generator(generator) {
		auxTable.resize(step - 1);
		for (double*& t : auxTable) {
			t = new double[num_states];
		}

		int numElems = (num_markers - 1) / step + 1;
		mainTable.resize(numElems);
		for (double*& t : mainTable) {
			t = new double[num_states];
		}

		totalAnchor = dirup ? 0 : num_markers - 1;
	}

	double* operator [](int i) {
		const int dir = dirup ? 1 : -1;
		const int i2 = totalAnchor + i * dir;

		const int anchor = i2 / step;
		const int auxIndex = (i2 % step) - 1;
		if (auxIndex == -1) {
			return mainTable[anchor];
		}

		if (anchor != auxAnchor) {
			double* prev = mainTable[anchor];
			for (int j = 0; j < auxTable.size(); j++) {
				int jorig = totalAnchor + (anchor * step + j + 1) * dir;
				if (jorig < 0 || jorig >= num_markers) break;
				(parent->*generator)(jorig, prev, auxTable[j]);
				prev = auxTable[j];
			}
			auxAnchor = anchor;
		}

		return auxTable[auxIndex];
	}

	void fillAllButFirst() {
		int dir = dirup ? 1 : -1;
		for (int j = 1; j < mainTable.size(); j++) {
			int jorig = totalAnchor + (j * step) * dir;			

			// Will trigger auxTable filling as needewd
			(parent->*generator)(jorig, (*this)[jorig - dir], mainTable[j]);
		}
	}

	~StepMemoizer() {
		for (double* t : mainTable) {
			delete[] t;
		}
		for (double* t : auxTable) {
			delete[] t;
		}
	}
};

/**

Main haplotype phasing and imputation functionality.

 */
class HaplotypePhaser {

public:

//	char ** haplotypes;
//	char ** genotypes;

	int prob_precision = 10;
	int curr_hap;

	MatrixXc haplotypes;
	MatrixXf sample_gls;



	Pedigree ped;
	float error;
	float Ne;
	double pop_const;

	std::vector<double> distances;

	vector<ChromosomePair> states;

	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	int num_haps;


	~HaplotypePhaser();
	void LoadData(const String &ref_file, const String &sample_file, const String &map_file);

	//private:

	int num_states;
	int num_markers;

	int num_ref_inds;
	int num_sample_inds;

	int num_inds;


	// index of the sample we are curently updating haplotype estimates for
	int current_sample;
	// index of current_sample's first haplotype in haplotypes matrix
	int current_sample_hap_index;

	int ChromStateToHaplotypeIndex(int chrom_state);
	void UpdateCurrentSample(int sample);
	void SetNumHaps(int num_haps);
	void FillStates();

	void CalcSingleScaledForward(int marker, const double* prev, double* now);
	void CalcSingleScaledBackward(int marker, const double* prev, double* now);

	StepMemoizer<HaplotypePhaser> s_forward;
	StepMemoizer<HaplotypePhaser> s_backward;
	double * normalizersf = 0;
	double * normalizersb = 0;


	void AllocateMemory();

	void CalcEmissionProbs(int marker, double * probs);

	void InitPriorScaledForward();
	void InitPriorScaledBackward();

	void CalcScaledForward();
	void CalcScaledBackward();
	void GetMLHaplotypes(int * ml_states);
	vector<int>  GetMLStates();
	vector<vector<double>> GetPosteriorStats(const char * filename, bool print);
	vector<vector<double>>  ReadPosteriorStats(const char * filename);

	void UpdateSampleHaplotypes(vector<int> ml_states);

	void PrintGenotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);
	void PrintHaplotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);

};


#endif
