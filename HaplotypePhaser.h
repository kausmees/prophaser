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
	__uint16_t first;
	__uint16_t second;

	ChromosomePair(__uint16_t a, __uint16_t b) {
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
	vector<double*> auxTable[2];
	T* parent;
	int num_markers;
	bool dirup;
	int step;
	int auxAnchor[2] = {-1, -1};
	int auxCeil[2] = {0, 0};	
	int totalAnchor;
	using gent = void (T::*)(int, const double*, double*);
	gent generator;
	

public:
	StepMemoizer() = default;

	// TODO: Copy assignment is dangerous.
	StepMemoizer(T* parent, int num_markers, int num_states, bool dirup, int step, gent generator) : parent(parent), num_markers(num_markers), dirup(dirup), step(step), generator(generator) {
		for (auto& subv : auxTable)
		{
			subv.resize(step - 1);
			for (double*& t : subv) {
			t = new double[num_states];
		}
		}		

		int numElems = (num_markers - 1) / step / step + 1;
		mainTable.resize(numElems);
		for (double*& t : mainTable) {
			t = new double[num_states];
		}

		totalAnchor = dirup ? 0 : num_markers - 1;
	}

	double* operator [](int i) {
		const int dir = dirup ? 1 : -1;
		const int i2 = totalAnchor + i * dir;

		
		const int auxIndex0 = (i2 % (step * step)) / step - 1;
		const int auxIndex1 = (i2 % (step)) - 1;
		const int midAnchor = i2 / step / step;
		double* prev = mainTable[midAnchor];
		if (auxIndex0 == -1 && auxIndex1 == -1) {
			return prev;
		}

		const int anchor = i2 / step;
		if (anchor != auxAnchor[1] || auxIndex1 >= auxCeil[1]) {
			if (auxIndex0 != -1)
			{
				if (midAnchor != auxAnchor[0])
				{	
					auxAnchor[0] = midAnchor;			
					auxCeil[0] = 0;
				}
					for (int& j = auxCeil[0]; j <= auxIndex0; j++) {
						int jorig = totalAnchor + (auxAnchor[0] * step * step + (j + 1) * step) * dir;
						if (jorig < 0 || jorig >= num_markers) break;
						(parent->*generator)(jorig, (*this)[jorig - dir], auxTable[0][j]);
					}
			
				prev = auxTable[0][auxIndex0];
		}

			if (auxIndex1 != -1)
			{
				if (auxAnchor[1] != anchor)
				{
					auxAnchor[1] = anchor;
					auxCeil[1] = 0;
				}
				if (auxCeil[1] != 0)
				{
					prev = auxTable[1][auxCeil[1] - 1];
				}
				for (int& j = auxCeil[1]; j <= auxIndex1; j++) {
					int jorig = totalAnchor + (auxAnchor[1] * step + j + 1) * dir;
				if (jorig < 0 || jorig >= num_markers) break;
					(parent->*generator)(jorig, prev, auxTable[1][j]);
					prev = auxTable[1][j];
			}
			}
		}

		if (auxIndex1 == -1) {
			return prev;
		}

		return auxTable[1][auxIndex1];
	}

	void fillAllButFirst() {
		int dir = dirup ? 1 : -1;
		for (int j = 1; j < mainTable.size(); j++) {
			int jorig = totalAnchor + (j * step * step) * dir;			

			// Will trigger auxTable filling as needewd
			(parent->*generator)(jorig, (*this)[jorig - dir], mainTable[j]);
		}
	}

	~StepMemoizer() {
		for (double* t : mainTable) {
			delete[] t;
		}
		for (auto& subv : auxTable)
		{
			for (double* t : subv) {
			delete[] t;
		}
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


	vector <double> sample_gls;
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
	void LoadReferenceData(const String &ref_file, String &map_file);
	void LoadSampleData(const String &sample_file,  int sample_index);


	//private:

	int num_states;
	int num_markers;

	int num_ref_inds;
	int num_inds;
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
	vector<vector<double>> GetPosteriorStats(const char * filename, bool print);
	vector<vector<double>>  ReadPosteriorStats(const char * filename);

	void PrintGenotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);
    void PrintPostGenotypesToVCF(vector<vector<int>> & ml_genotypes, vector<vector<vector<double>>> & ml_postprobs, const char * out_file, const char * sample_file, const char * vcf_template);
    void PrintHaplotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);

};


#endif
