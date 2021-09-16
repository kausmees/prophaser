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

using phaserreal = float;

typedef Eigen::Matrix<phaserreal,Dynamic, Dynamic> MatrixXpreal;


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
	vector<phaserreal*> mainTable;
	vector<phaserreal*> auxTable[2];
	T* parent;
	int num_markers;
	int num_states;
	bool dirup;
	int step;
	int auxAnchor[2] = {-1, -1};
	int auxCeil[2] = {0, 0};	
	int totalAnchor;
	using gent = void (T::*)(int, const phaserreal*, phaserreal*);
	gent generator;
	

public:
	StepMemoizer() = default;

	// TODO: Copy assignment is dangerous.
	StepMemoizer(T* parent, int num_markers, int num_states, bool dirup, int step, gent generator) : parent(parent), num_markers(num_markers), dirup(dirup), step(step), generator(generator), num_states(num_states) {		
		for (auto& subv : auxTable)
		{
			subv.resize(step - 1);
			for (phaserreal*& t : subv) {
				t = new phaserreal[num_states];
		}
		}		

		int numElems = (num_markers - 1) / step / step + 1;
		mainTable.resize(numElems);
		for (phaserreal*& t : mainTable) {
			t = new phaserreal[num_states];
		}

		totalAnchor = dirup ? 0 : num_markers - 1;
	}

	phaserreal* rawindexer(int i) {
		const int dir = dirup ? 1 : -1;
		const int i2 = totalAnchor + i * dir;

		
		const int auxIndex0 = (i2 % (step * step)) / step - 1;
		const int auxIndex1 = (i2 % (step)) - 1;
		const int midAnchor = i2 / step / step;
		phaserreal* prev = mainTable[midAnchor];
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

	phaserreal* operator [](int i) {
		phaserreal* data = rawindexer(i);
		return data;
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
		for (phaserreal* t : mainTable) {
			delete[] t;
		}
		for (auto& subv : auxTable)
		{
			for (phaserreal* t : subv) {
			delete[] t;
		}
	}
	}
};

struct HapSummer;

/**

Main haplotype phasing and imputation functionality.

 */
class HaplotypePhaser {
	phaserreal* emission_probs;
	HapSummer* hapSum;
	phaserreal case_1, case_2, case_3, case_4, case_5;
public:

//	char ** haplotypes;
//	char ** genotypes;

	int prob_precision = 10;
	int curr_hap;
	MatrixXc haplotypes;


	MatrixXpreal sample_gls;
	Pedigree ped;
	float error;
	float Ne;
	phaserreal pop_const;

	std::vector<phaserreal> distances;

	vector<ChromosomePair> states;

	// Number of reference haplotypes that will be considered when handling one sample
	// The num_haps first haplotypes in the matrix haplotypes will be used.
	int num_haps;


	~HaplotypePhaser();

	void LoadData(const String &ref_file, const String &sample_file, const String &map_file);

	void LoadReferenceData(const String &ref_file, String &map_file);
	void LoadSampleData(const String &sample_file,  int sample_index);


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

	void CalcSingleScaledForward(int marker, const phaserreal* prev, phaserreal* now);
	void CalcSingleScaledBackward(int marker, const phaserreal* prev, phaserreal* now);

	StepMemoizer<HaplotypePhaser> s_forward;
	StepMemoizer<HaplotypePhaser> s_backward;
	phaserreal* normalizersf = 0;
	phaserreal* normalizersb = 0;


	void AllocateMemory();

	phaserreal CalcEmissionProb(const ChromosomePair cp, int marker);
	void CalcEmissionProbs(int marker, phaserreal* probs);

	void InitPriorScaledForward();
	void InitPriorScaledBackward();

	void CalcScaledForward();
	void CalcScaledBackward();
	void GetMLHaplotypes(int * ml_states);
	vector<vector<phaserreal>> GetPosteriorStats(const char * filename, bool print);
	vector<vector<phaserreal>>  ReadPosteriorStats(const char * filename);

	void PrintGenotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);
    void PrintPostGenotypesToVCF(vector<vector<int>> & ml_genotypes, vector<vector<vector<phaserreal>>> & ml_postprobs, const char * out_file, const char * sample_file, const char * vcf_template);
    void PrintHaplotypesToVCF(vector<vector<int>> & ml_genotypes, const char * out_file, const char * sample_file, const char * vcf_template);

};


#endif
