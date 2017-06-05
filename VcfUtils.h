#ifndef __VCFUTILS_H__
#define __VCFUTILS_H__

#include <fstream>
#include "VcfRecordGenotype.h"
#include "VcfFileReader.h"
#include "Pedigree.h"
#include "GenoUtils.h"

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
void LoadIndividuals(Pedigree &ped, const String &ref_file, const String & sample_file, int sample_id);
void LoadHaplotypes(const String &file_name, const Pedigree &ped, char** haplotypes);
void LoadGenotypeLikelihoods(const String &file_name, const Pedigree &ped, vector<double> & genotypes, int sample_ind_file);
int GetMarkerPos(int marker_id);
void LoadGeneticMap(const char * file_name, const Pedigree &ped, vector<double> &thetas);


}


#endif
