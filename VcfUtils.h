#ifndef __VCFUTILS_H__
#define __VCFUTILS_H__

#include <fstream>
#include "VcfRecordGenotype.h"
#include "VcfFileReader.h"
#include "Pedigree.h"
#include "GenoUtils.h"
#include <Eigen/Dense>

using Eigen::Dynamic;
typedef Eigen::Matrix<int,Dynamic, Dynamic> MatrixXc;

using namespace Eigen;

/**
 *
 * Utility methods for handling VCF files
 *
 */
namespace VcfUtils {

extern int max_pl;
extern std::vector<int> unphased_marker_subset;

//typedef std::set<std::vector<String>>




void LoadReferenceMarkers(const String &file_name);
void LoadReferenceIndividuals(Pedigree &ped, const String &ref_file);
void LoadSampleIndividual(Pedigree &ped, const String & sample_file, int sample_id);


void LoadHaplotypes(const String &file_name, const Pedigree &ped, char** haplotypes);
void LoadHaplotypes(const String &file_name, const Pedigree &ped, MatrixXc & haplotypes );
void FillSampleHaplotypes(const Pedigree &ped, MatrixXc & haplotypes, int index);
void LoadGenotypeLikelihoods(const String &file_name, const Pedigree &ped, vector<double> & genotypes, int sample_ind_file);
int GetMarkerPos(int marker_id);
void LoadGeneticMap(const char * file_name, const Pedigree &ped, vector<double> &thetas);

void writeVectorToCSV(const char* file_name, std::vector<vector<double>> v, const char* opentype);
}


#endif
