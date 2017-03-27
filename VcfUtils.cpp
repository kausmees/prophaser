#include "VcfRecordGenotype.h"
#include "VcfFileReader.h"
#include "Pedigree.h"
#include "VcfUtils.h"


#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>


namespace VcfUtils{

HaplotypePair::HaplotypePair(std::vector<String> in1, std::vector<String> in2) {
	if(in1.size() != in2.size())
		throw std::invalid_argument("Haplotypes not of same length.");

	length = in1.size();
	h1=in1;
	h2=in2;
}

void HaplotypePair::print(){
	for (auto & t : h1) {
		printf("%s", t.ToUpper().c_str());

	}
	printf("\n");
	for (auto & t : h2) {
		printf("%s", t.ToUpper().c_str());
	}
	printf("\n");
}

bool HaplotypePair::isEqual(HaplotypePair other) {
	return HaplotypePair::pairEqual(h1,h2,other.h1,other.h2);
};

double HaplotypePair::switchError(HaplotypePair other) {
	double count = 0.0;

	if(length != other.length)
		return -1.0;

	std::vector<String>::iterator this_1 = h1.begin()+1;
	std::vector<String>::iterator this_2 = h2.begin()+1;
	std::vector<String>::iterator other_1 = other.h1.begin()+1;
	std::vector<String>::iterator other_2 = other.h2.begin()+1;

	while(this_1 != h1.end()) {
		if(!pairEqual(subVector(this_1-1,2), subVector(this_2-1,2), subVector(other_1-1,2),subVector(other_2-1,2))) {
			count++;
		}
		this_1++;
		this_2++;
		other_1++;
		other_2++;
	}
	return count / (h1.size() - 1);
}

/**
 * Compare genotype implied by h1 and h2 to that of haplotype pair h_true
 * at marker indices specified in positions.
 *
 * Return percentage of incorrect markers.
 *
 */
double HaplotypePair::imputeError(HaplotypePair h_true, std::vector<int> positions){
	double count = 0.0;
	for(int & i : positions) {
		if(!pairEqual(h1[i], h2[i], h_true.h1[i], h_true.h2[i])) {
			count++;
		}
	}
	return count / positions.size();
}




int max_pl = 255;
std::vector<int> unphased_marker_subset;

/**
 * Loads all correct markers in file_name into Pedigree.
 *
 * A correct marker :
 * Is a SNP with 1 alternative allele and with no missing base in the reference or alt allele.
 * Has no missing or unphased genotypes.
 *
 *
 * Throws error if unphased or missing genotypes present in any sample.
 */
void LoadReferenceMarkers(const String &file_name){

	printf("Loading markers from file %s \n", file_name.c_str());

	VcfFileReader reader;
	VcfHeader header;
	reader.open(file_name, header);

	const char * marker_name;

	VcfRecord record;
	while(reader.readRecord(record)){
		//A marker's name in the pedigree is its chromosome and position.
		//So they can be identified as same when originating from different files

		std::stringstream ss;
		ss << record.getChromStr() << ":" << record.get1BasedPosition();
		marker_name = ss.str().c_str();
		ss.clear();
		//TODO should use sprintf instead, faster? know chrom and pos wont be non-null terminated
		// or large to cause buffer overflow issues?
		//sprintf(marker_id, "%s:%d",record.getChromStr(), record.get1BasedPosition());

		if(!record.allPhased()) {
			//TODO throw error.
			printf("ERROR: NOT ALL PHASED!\n");
			printf("all phased = %d \n", record.allPhased());
			printf("all genotypes = %d \n", record.hasAllGenotypeAlleles());
		}

		// if marker is a monomorphic SNP with no missing bases
		if(record.getNumRefBases() == 1 && record.getNumAlts() == 1 && strcmp(record.getRefStr(),"N") != 0 && strlen(record.getRefStr()) == 1 && strlen(record.getAltStr()) == 1) {

			int marker_id = Pedigree::GetMarkerID(marker_name);
			printf("added marker (%d, %s) to pedigree \n", marker_id, marker_name);
			String ref_base = String(record.getRefStr()[0]);
			String alt_base = String(record.getAltStr()[0]);

			printf("ref base = %s \n", ref_base.c_str());
			printf("alt base = %s \n", alt_base.c_str());


			//TODO perform check on ref and alt values?
			int ref = Pedigree::LoadAllele(marker_id, ref_base);
			int alt = Pedigree::LoadAllele(marker_id, alt_base);

			//printf("ref str = %s  num = %d \n", ref_base.c_str(), ref);
			//printf("alt str = %s  num = %d \n", alt_base.c_str(), alt);
		}
	}
	unphased_marker_subset.resize(Pedigree::markerCount,0);
	reader.close();
};

/**
 * Add individuals from file_name to the pedigree ped.
 * TODO possibly no need to have individuals in pedigree.

 */
void LoadIndividuals(Pedigree &ped, const String &ref_file, const String &sample_file, int sample_index) {

	VcfFileReader reader;
	VcfHeader header;
	reader.open(ref_file, header);

	int num_samples = header.getNumSamples();
	printf("Num reference inds : %d \n", num_samples);

	if(num_samples == 0) {
		//TODO add exceptions
		printf("ERROR: No reference individuals in file");
	}

	for(int i = 0; i < num_samples; i++) {
		ped.AddPerson(header.getSampleName(i), header.getSampleName(i), "0", "0", 0, 1);
	}

	reader.close();
	printf("Num reference inds in ped: %d \n", ped.count);

	reader.open(sample_file, header);
	ped.AddPerson(header.getSampleName(sample_index), header.getSampleName(sample_index), "0", "0", 0, 1);
	reader.close();

};

/**
 * Load all alleles from file_name at markers present in pedigree into haplotypes.
 *
 * Assumes the file only contains samples that are phased at all markers, have no missing reference of alternative alleles,
 * and that the pedigree only contains monomorphic SNPs.
 *
 */
void LoadHaplotypes(const String &file_name, const Pedigree &ped, char** haplotypes) {

	printf("Loading haplotypes from file %s \n", file_name.c_str());

	VcfFileReader reader;
	VcfHeader header;
	reader.open(file_name, header);

	const char * marker_name;
	VcfRecord record;
	while(reader.readRecord(record)){
		std::stringstream ss;
		ss << record.getChromStr() << ":" << record.get1BasedPosition();
		marker_name = ss.str().c_str();
		int marker_id = Pedigree::LookupMarker(marker_name);
		if(marker_id >= 0){

			for (int ind = 0; ind < ped.count; ind++) {
				int i0 = record.getGT(ind,0);
				int i1 = record.getGT(ind,1);

				haplotypes[ind*2][marker_id] = i0;
				haplotypes[ind*2 + 1][marker_id] = i1;

			}
		}
	}
	reader.close();
};

/**
 * Load genotypes of individual sample_index_file from file_name into last row of genotypes.
 * Markers present in pedigree but not in file_name are given genotype phreds 0.
 *
 * temproary: fill all other genotypes with 0
 */
void LoadGenotypeLikelihoods(const String &file_name, const Pedigree &ped, char** genotypes, int sample_index_file) {
	printf("Loading genotype likelihoods from file %s \n", file_name.c_str());

	int sample_index_ped = ped.count-1;
	std::string lformat;
	const std::string *likelihood;
	const char * marker_name;
	int pl_00, pl_01, pl_11;
	int num_common_markers = 0;
	int num_total_markers = 0;

	std::string sub00, sub01, sub11;

	VcfFileReader reader;
	VcfHeader header;
	reader.open(file_name, header);
	VcfRecord record;


	while(reader.readRecord(record)){
		num_total_markers += 1;

		std::stringstream ss;
		ss << record.getChromStr() << ":" << record.get1BasedPosition();
		marker_name = ss.str().c_str();

		int marker_id = Pedigree::LookupMarker(marker_name);

		if(marker_id >= 0){
			printf("DOING MARKER ID %d\n", marker_id);
			// if marker is a monomorphic SNP with no missing bases
			if(record.getNumRefBases() != 1 || record.getNumAlts() != 1) {
				//TODO throw error
				printf("ERROR: Sample marker is not monomorphic SNP\n");
				//TODO check how this affects imputation. this marker will behave as one we have no info on in sample,
				//do we want it to be imputed as if it was missing in sample?
				continue;

			}
			String sample_ref_base = String(record.getRefStr()[0]);
			String sample_alt_base = String(record.getAltStr()[0]);

			MarkerInfo * mi = Pedigree::GetMarkerInfo(marker_id);

			String phased_ref_base = mi->GetAlleleLabel(1);
			String phased_alt_base = mi->GetAlleleLabel(2);

			if(sample_ref_base != phased_ref_base || sample_alt_base != phased_alt_base){

				//TODO throw error or warning
				printf("ERROR: Sample alleles do not match phased reference's alles\n");
				//TODO check how this affects imputation. this marker will behave as one we have no info on in sample,
				//do we want it to be imputed as if it was missing in sample?
				continue;


			}

			printf("PED REF = %s  SAMPLE REF = %s \n",  phased_ref_base.c_str(), sample_ref_base.c_str());
			printf("PED ALT = %s  SAMPLE ALT = %s \n",  phased_alt_base.c_str(), sample_alt_base.c_str());


			lformat = "GL";
			//check if this record has GL or PL
			likelihood = record.getGenotypeInfo().getString(lformat,0);
			if(likelihood == NULL) {
				lformat = "PL";
				likelihood = record.getGenotypeInfo().getString(lformat,0);

				if(likelihood == NULL) {
					//TODO throw error
					printf("ERROR: RECORD WITH NO PL OR GL FOUND\n");
					//TODO check how this affects imputation. this marker will behave as one we have no info on in sample,
					//do we want it to be imputed as if it was missing in sample?
					continue;
				}
			}


			// if VCF record has GL/PL but sample has missing field with :    then likelihood = "."
			// if VCF record has GL/PL but sample has missing field with nothing    then likelihood = ""
			// if VCF record does not have GL/PL then likelihood is NULL pointer
			likelihood = record.getGenotypeInfo().getString(lformat,sample_index_file);

			if(*likelihood == "" || *likelihood == ".") {
				//TODO throw error
				printf("ERROR: INDIVIDUAL WITH NO LIKELIHOOD FOUND \n");
			}


			std::istringstream iss(*likelihood);

			std::getline(iss, sub00, ',');
			std::getline(iss, sub01, ',');
			std::getline(iss, sub11, ',');


			pl_00 = (lformat == "PL") ? atoi(sub00.c_str()) : static_cast<int>(-10.0 * atof(sub00.c_str()));
			pl_01 = (lformat == "PL") ? atoi(sub01.c_str()) : static_cast<int>(-10.0 * atof(sub01.c_str()));
			pl_11 = (lformat == "PL") ? atoi(sub11.c_str()) : static_cast<int>(-10.0 * atof(sub11.c_str()));

			if(pl_00 < 0 || pl_01 < 0 || pl_11 < 0) {
				//TODO throw error
				printf("ERROR: NEGATIVE PL \n");
				//TODO check how this affects imputation. this marker will behave as one we have no info on in sample,
				//do we want it to be imputed as if it was missing in sample?
				continue;
			}

			pl_00 = std::min(pl_00, max_pl);
			pl_01 = std::min(pl_01, max_pl);
			pl_11 = std::min(pl_11, max_pl);

			genotypes[sample_index_ped][marker_id*3] = pl_00;
			genotypes[sample_index_ped][marker_id*3+1] = pl_01;
			genotypes[sample_index_ped][marker_id*3+2] = pl_11;

			unphased_marker_subset[marker_id] = 1;
			num_common_markers += 1;
		}
	}
	reader.close();

	for(int marker_id = 0; marker_id < ped.markerCount; marker_id++){
		if(!unphased_marker_subset[marker_id]){
			genotypes[sample_index_ped][marker_id*3] = 0;
			genotypes[sample_index_ped][marker_id*3+1] = 0;
			genotypes[sample_index_ped][marker_id*3+2] = 0;
		}

		//TODO temproary, fill in reference genotypes as 0
		for(int i = 0; i < ped.count -1 ; i++){
			genotypes[i][marker_id*3] = 0;
			genotypes[i][marker_id*3+1] = 0;
			genotypes[i][marker_id*3+2] = 0;

		}
	}
	printf("Num reference markers : %d, Num sample markers: %d, %d of which are also in reference. \n", ped.markerCount, num_total_markers, num_common_markers);
}


int GetMarkerPos(int marker_id){
	std::string marker_name = Pedigree::GetMarkerInfo(marker_id)->name.c_str();
	std::size_t delim = marker_name.find(":");
	return std::stoi(marker_name.substr(delim+1));
}


/**
 * Fill thetas with crossover probabilities so that
 * thetas[m] = prob of crossover happening in interval between m-1 and m.
 *
 *
 */
void LoadGeneticMap(const char *file_name, const Pedigree &ped, vector<double> &distances) {
	int ped_pos;
	int map_pos;
	double dist;
	double prev_dist;
	int c = 0;
	char marker_name [256];
	int result;
	//TODO dont hardcode Ne
	double pop_const = 4*11418;
	FILE * mapstream = fopen(file_name, "r");

	result = fscanf(mapstream, "%255s %d %lf", marker_name, &map_pos, &prev_dist);

	for (int i = 1; i < Pedigree::markerCount && result == 3; i++) {
		ped_pos = GetMarkerPos(i);
		while(map_pos < ped_pos && result == 3) {
			result = fscanf(mapstream, "%255s %d %lf", marker_name, &map_pos, &dist);
			c +=1;
		}
		if(map_pos == ped_pos) {

			//TODO what should it be set to if it is zero??
			distances[i] = (dist - prev_dist > 0.00000001) ? (dist - prev_dist)*pop_const : 0.001*pop_const;
			prev_dist = dist;
		}
		else{
			distances[i] = distances[i-1];
			//	TODO throw warning
			printf("WARNING: Marker %d is not present in crossover file. \n", GetMarkerPos(i));
			continue;
		}
	}

	for (int i = 0; i < Pedigree::markerCount; i++) {
		printf("Marker : %d, Theta: %f \n", GetMarkerPos(i), distances[i]);
	}
	printf("%d lines read from map \n", c);

}

/**
 * Load haplotypes of individual sample_ind from file_name at markers
 * specified in pedigree.
 */
HaplotypePair loadHaplotypesFromVCF(const char * file_name, int sample_ind) {
	std::vector<String> h1;
	std::vector<String> h2;

	VcfFileReader reader;
	VcfHeader header;
	reader.open(file_name, header);
	int g1, g2;

	const char * marker_name;
	VcfRecord record;
	while(reader.readRecord(record)){
		std::stringstream ss;
		ss << record.getChromStr() << ":" << record.get1BasedPosition();
		marker_name = ss.str().c_str();
		int marker_id = Pedigree::LookupMarker(marker_name);
		if(marker_id >= 0){

			g1 = record.getGT(sample_ind,0);
			g2 = record.getGT(sample_ind,1);
			if(marker_id == 14 || marker_id ==13 || marker_id == 15){
				printf("%dth: %d %d \n", marker_id, g1, g2);
				printf("Alleles: %s %s \n", Pedigree::GetMarkerInfo(marker_id)->GetAlleleLabel(g1+1).c_str(), Pedigree::GetMarkerInfo(marker_id)->GetAlleleLabel(g2+1).c_str());

			}
			h1.push_back(Pedigree::GetMarkerInfo(marker_id)->GetAlleleLabel(g1+1));
			h2.push_back(Pedigree::GetMarkerInfo(marker_id)->GetAlleleLabel(g2+1));
		}
	}
	reader.close();
	HaplotypePair haps(h1,h2);
	return haps;
}

HaplotypePair loadHaplotypesFromMACH(const char * file_name) {
	std::vector<String> h1;
	std::vector<String> h2;

	FILE * hapstream = fopen(file_name, "r");
	int c;

	while(!feof(hapstream) && c != ' '){
		c = fgetc(hapstream);
	}
	c = fgetc(hapstream);
	while(!feof(hapstream) && c != ' '){
		c = fgetc(hapstream);
	}
	c = fgetc(hapstream);
	while(!feof(hapstream) && c != '\n'){
		h1.push_back(String(char(c)));
		c = fgetc(hapstream);
	}
	while(!feof(hapstream) && c != ' '){
		c = fgetc(hapstream);
	}
	c = fgetc(hapstream);
	while(!feof(hapstream) && c != ' '){
		c = fgetc(hapstream);
	}
	c = fgetc(hapstream);
	while(!feof(hapstream) && c != '\n'){
		h2.push_back(String(char(c)));
		c = fgetc(hapstream);

	}

	//	while(!feof(hapstream) && *fgets(buffer, 1, hapstream) != '\n'){
	//		hp.h2.push_back(String(buffer[0]));
	//	}
	HaplotypePair hp(h1,h2);

	return hp;


}

void CompareHaplotypes(){

	char ha1 [Pedigree::markerCount+1];
	char ha2 [Pedigree::markerCount+1];
	char hb1 [Pedigree::markerCount+1];
	char hb2 [Pedigree::markerCount+1];


	std::string filea = "/home/kristiina/Projects/Data/tests/mine.txt";
	std::string fileb = "/home/kristiina/Projects/Data/tests/mach.txt";

	std::ifstream hapsa (filea);
	std::ifstream hapsb (fileb);

	hapsa.getline(ha1, Pedigree::markerCount+1);
	hapsa.getline(ha2, Pedigree::markerCount+1);
	hapsb.getline(hb1, Pedigree::markerCount+1);
	hapsb.getline(hb2, Pedigree::markerCount+1);

	std::cout << ha1 << "\n";
	std::cout << ha2 << "\n";
	std::cout << hb1 << "\n";
	std::cout << hb2 << "\n";


	for(int i = 0; i < Pedigree::markerCount; i++) {
		if(ha1[i] != hb1[i]){
			printf("ha1 and hb1 differ at position %d \n", i);
		}

	}

	for(int i = 0; i < Pedigree::markerCount; i++) {
		if(ha1[i] != hb2[i]){
			printf("ha1 and hb2 differ at position %d \n", i);
		}
	}

	for(int i = 0; i < Pedigree::markerCount; i++) {
		if(ha2[i] != hb1[i]){
			printf("ha2 and hb1 differ at position %d \n", i);
		}
	}

	for(int i = 0; i < Pedigree::markerCount; i++) {
		if(ha2[i] != hb2[i]){
			printf("ha2 and hb2 differ at position %d \n", i);
		}
	}


}


}
