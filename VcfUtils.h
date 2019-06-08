#ifndef __VCFUTILS_H__
#define __VCFUTILS_H__

#include <fstream>
#include "VcfRecordGenotype.h"
#include "VcfFileReader.h"
#include "Pedigree.h"
//#include "GenoUtils.h"
#include <Eigen/Dense>
#include <numeric>

#include <algorithm>
#include "VcfFileWriter.h"



using Eigen::Dynamic;
typedef Eigen::Matrix<int,Dynamic, Dynamic> MatrixXc;

using namespace Eigen;
using namespace std;
/**
 *
 * Utility methods for handling VCF files
 *
 */
namespace VcfUtils {

extern int max_pl;
extern std::vector<int> unphased_marker_subset;

//typedef std::set<std::vector<String>>

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
// initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  stable_sort(idx.begin(), idx.end(),
		         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});


  return idx;
};


void LoadReferenceMarkers(const String &file_name);
void LoadReferenceIndividuals(Pedigree &ped, const String &ref_file);
void LoadSampleIndividual(Pedigree &ped, const String & sample_file, int sample_id);


void LoadHaplotypes(const String &file_name, const Pedigree &ped, char** haplotypes);
void LoadHaplotypes(const String &file_name, const Pedigree &ped, MatrixXc & haplotypes );
void FillSampleHaplotypes(const Pedigree &ped, MatrixXc & haplotypes, int index);
void LoadGenotypeLikelihoods(const String &file_name, const Pedigree &ped, vector<double> & genotypes, int sample_ind_file);
int GetMarkerPos(int marker_id);
void LoadGeneticMap(const char * file_name, const Pedigree &ped, vector<double> &thetas);
void writeVectorToCSV(const char* file_name, const std::vector<vector<double>>& v, const char* opentype);
vector<double> get_GL(VcfHeader& header, VcfRecord& record, int sample);
std::string get_likelihood_format(VcfHeader& header);
}


#endif
