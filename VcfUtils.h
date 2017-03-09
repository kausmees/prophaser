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
};

void LoadReferenceMarkers(const String &file_name);
void LoadIndividuals(Pedigree &ped, const String &ref_file, const String & sample_file, int sample_id);
void LoadHaplotypes(const String &file_name, const Pedigree &ped, char** haplotypes);
void LoadGenotypeLikelihoods(const String &file_name, const Pedigree &ped, char** genotypes, int sample_ind_file);
void LoadGeneticMap(const char * file_name, const Pedigree &ped, vector<double> &thetas);
int GetMarkerPos(int marker_id);
void CompareHaplotypes();
}
