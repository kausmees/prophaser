#include "Pedigree.h"
#include "VcfUtils.h"
#include <vector>
#include <chrono>

/**

Main haplotype phasing and imputation functionality.

*/
class HaplotypePhaser {

public:
	char ** haplotypes;
	char ** genotypes;
	Pedigree ped;
	float theta;
	float error;
	int Ne;

	std::vector<double> errors;
	std::vector<double> distances;

	// phred_probs[i] contains linear-scale genotype likelihood
	float * phred_probs;

	~HaplotypePhaser();
	void LoadData(const String &ref_file, const String &sample_file, int sample_index);
	void LoadReferenceData(const String &ref_file, const String &sample_file, int sample_index);
	void LoadSampleData(const String &ref_file, const String &sample_file, int sample_index);



//private:

	int num_states;
	int num_markers;

	//TODO
	int num_inds;

	float ** transition_probs;


	// remove
	float ** forward;
	float ** backward;
	float ** s_backward;

	float ** s_forward;
	float * normalizers;

	void AllocateMemory();
	//TODO remove, for testing
	void DeAllocateMemory();
	void CalcGenoEmissions();
	double GetEmissionProb(int state, int marker);
	void CalcTransitionProbs(int marker, int marker_state, double * probs);
	void CalcEmissionProbs(int marker, float * probs);
	void InitPriorForward();
	void InitPriorBackward();

	//Old
	void CalcForward();
	void CalcBackward();
	void CalcPosteriorOld();


	void InitPriorScaledForward();
	void InitPriorScaledBackward();

	void CalcScaledForward();
	void CalcPosterior();

	//Old
	void CalcScaledBackward();

	void SampleHaplotypes(int * ml_states);
	void SampleHaplotypesNew(int * ml_states);
	void InferHaplotypes();

	HaplotypePair PrintHaplotypesToFile(int * states, const char * out_file);
};
