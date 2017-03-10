#include <fstream>

/**
 *
 * Utility methods for handling VCF files
 *
 */
namespace VcfUtils {

extern int max_pl;
extern std::vector<int> unphased_marker_subset;

struct HaplotypePair {
	std::vector<String> h1;
	std::vector<String> h2;

	void print(){
		for (auto & t : h1) {
			printf("%s", t.ToUpper().c_str());

		}
		printf("\n");
		for (auto & t : h2) {
			printf("%s", t.ToUpper().c_str());
		}
		printf("\n");
	}

	bool isEqual(HaplotypePair other) {
		return ((h1 == other.h1 && h2 == other.h2) || (h1 == other.h2 && h2 == other.h1));
	};

};

void LoadReferenceMarkers(const String &file_name);
void LoadIndividuals(Pedigree &ped, const String &ref_file, const String & sample_file, int sample_id);
void LoadHaplotypes(const String &file_name, const Pedigree &ped, char** haplotypes);
void LoadGenotypeLikelihoods(const String &file_name, const Pedigree &ped, char** genotypes, int sample_ind_file);
int GetMarkerPos(int marker_id);
void LoadGeneticMap(const char * file_name, const Pedigree &ped, vector<double> &thetas);
HaplotypePair loadHaplotypesFromVCF(const char * file_name, int sample_ind);
HaplotypePair loadHaplotypesFromMACH(const char * file_name);
void CompareHaplotypes();
}
