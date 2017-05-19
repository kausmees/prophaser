#include <chrono>
#include "HaplotypePhaser.h"
#include "Parameters.h"
#include "Pedigree.h"
#include "GenoUtils.h"

using namespace std;


int main(int argc, char ** argv){
//	String sample_file;
//	String ref_file;
//
//	ParameterList parameter_list;
//	LongParamContainer long_parameters;
//	long_parameters.addGroup("Input files");
//	long_parameters.addString("sampleVCF", &sample_file);
//	long_parameters.addString("referenceVCF", &ref_file);
//
//	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
//	parameter_list.Read(argc, argv);
//	parameter_list.Status();




//	String dir = "../../Data/1KGData/vcfs/chrom20/";
//	String data_id = "4";
//	String subset_id = "A";

	String dir = "../../Data/1KGData/vcfs/chrom20/";
	String data_id = "4";
	String subset_id = "B";


//	vector<double> coverages = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	vector<double> coverages = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

	double error = 0.001;
	String file_name_in;
	String file_name_out;

//	String sample_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_snps.vcf";
//	int sample_ind = 624;

	String sample_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_NA18909_snps.vcf";
	int sample_ind = 0;


	String ref_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_10inds_snps.vcf";

	HaplotypePhaser phaser;
	phaser.LoadReferenceData(ref_file, sample_file, sample_ind);

	HaplotypePair true_haps = loadHaplotypesFromVCF(sample_file, sample_ind);


	////////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	file_name_in = dir + subset_id + "/" + data_id + "_" + subset_id + "_" + "NA18909" + "_snps.vcf" ;
	file_name_out = "./Results/" + data_id + "_" + subset_id + "_" + "NA18909" + "_snps_phased_4" ;
	phaser. (ref_file, file_name_in, 0);

	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
		cout << "Starting Forward \n";
	phaser.CalcScaledForward();

		cout << "Starting Backward \n";
	phaser.CalcScaledBackward();

	cout << "Starting Posterior \n";
	phaser.CalcPosterior();

	chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

	cout << "Original: " << endl;
	cout << "Time = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;

	int * ml_states = new int[phaser.num_markers];

	cout << "Sampling " << endl;

	phaser.SampleHaplotypesNew(ml_states);
	cout << "Printing " << endl;

	HaplotypePair mine_haps = phaser.PrintHaplotypesToFile(ml_states, file_name_out);
	mine_haps.print();
	delete [] ml_states;



	for(auto cov : coverages) {
		file_name_in = dir + subset_id + "/" + data_id + "_" + subset_id + "_" + "NA18909" + "_snps_" + (to_string(cov)).c_str() + "_" + (to_string(error)).c_str() + ".vcf" ;
		file_name_out ="./Results/" + data_id + "_" + subset_id + "_" + "NA18909" + "_snps_" + (to_string(cov)).c_str() + "_" + (to_string(error)).c_str() + "_phased_4" ;
		phaser.LoadSampleData(ref_file, file_name_in, 0);

		chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//		cout << "Starting Forward \n";
		phaser.CalcScaledForward();

//		cout << "Starting Backward \n";
		phaser.CalcScaledBackward();

//		cout << "Starting Posterior \n";
		phaser.CalcPosterior();
		chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
		cout << "Coverage: " << cov << endl;
		cout << "Time = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;


		int * ml_states = new int[phaser.num_markers];
		phaser.SampleHaplotypesNew(ml_states);
		HaplotypePair mine_haps = phaser.PrintHaplotypesToFile(ml_states, file_name_out);
		delete [] ml_states;

	}




	///////////////////////////////////////////////////////////////////////////////////////////////////////////

//	String orig_sample_file = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps.vcf";
//	String ref_file_ = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_10inds_snps.vcf";
//
//
//	const char * file4_lowcov1 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.100000_0.001000.vcf";
//	const char * file4_lowcov2 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.200000_0.001000.vcf";
//	const char * file4_lowcov3 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.300000_0.001000.vcf";
//	const char * file4_lowcov4 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.400000_0.001000.vcf";
//	const char * file4_lowcov5 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.500000_0.001000.vcf";
//	const char * file4_lowcov6 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.600000_0.001000.vcf";
//	const char * file4_lowcov7 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.700000_0.001000.vcf";
//	const char * file4_lowcov8 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.800000_0.001000.vcf";
//	const char * file4_lowcov9 = "../../Data/1KGData/vcfs/chrom20/A/First/4_A_snps_0.900000_0.001000.vcf";

//	String orig_sample_file = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps.vcf";
//	String ref_file = "../../Data/1KGData/vcfs/chrom20/B/4_B_10inds_snps.vcf";
//
//
//	const char * file4_lowcov1 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.100000_0.001000.vcf";
//	const char * file4_lowcov2 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.200000_0.001000.vcf";
//	const char * file4_lowcov3 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.300000_0.001000.vcf";
//	const char * file4_lowcov4 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.400000_0.001000.vcf";
//	const char * file4_lowcov5 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.500000_0.001000.vcf";
//	const char * file4_lowcov6 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.600000_0.001000.vcf";
//	const char * file4_lowcov7 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.700000_0.001000.vcf";
//	const char * file4_lowcov8 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.800000_0.001000.vcf";
//	const char * file4_lowcov9 = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps_0.900000_0.001000.vcf";
//
//
//	const char * sample_name = "NA18909";
//
//	PedFile ped_file = PedFile(2);
//
//
//	int sample_ind = get_sample_index(orig_sample_file.c_str(), sample_name);
//	printf("Ind index = %d \n", sample_ind);
//
//
//
//	vector<String> sample_files =  {file4_lowcov1,
//									file4_lowcov2,
//									file4_lowcov3,
//									file4_lowcov4,
//									file4_lowcov5,
//									file4_lowcov6,
//									file4_lowcov7,
//									file4_lowcov8,
//									file4_lowcov9,
//									orig_sample_file};
////
////
//	vector<double> switch_error_mine;
//	vector<double> switch_error_mach;
//
//	vector<double> geno_error_mine;
//	vector<double> geno_error_mach;
//
//	HaplotypePhaser phaser;
//	phaser.LoadReferenceData(ref_file, orig_sample_file, sample_ind);
//
//
//	for(auto sample_file : sample_files) {
////		printf("in for \n");
//
//		phaser.LoadSampleData(ref_file, sample_file, sample_ind);
//		HaplotypePair true_haps = loadHaplotypesFromVCF(orig_sample_file.c_str(), 624);
//
//		chrono::steady_clock::time_point begin = chrono::steady_clock::now();
////		cout << "Starting Forward \n";
//		phaser.CalcScaledForward();
//
////		cout << "Starting Backward \n";
//		phaser.CalcScaledBackward();
//
////		cout << "Starting Posterior \n";
//		phaser.CalcPosterior();
//		chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//		cout << "Time = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
//
//
//		int * ml_states = new int[phaser.num_markers];
//		phaser.SampleHaplotypesNew(ml_states);
//
//		HaplotypePair mine_haps = phaser.PrintHaplotypes(ml_states);
//
//		delete [] ml_states;
//
//		string file_name = sample_file.c_str();
//		string::size_type pos4 = file_name.rfind("/4", file_name.length()) +1;
//		string mach_end = "_NA18909.out";
//		string subst = file_name.substr(pos4);
//		string::size_type posdot = subst.rfind('.', file_name.length());
//		subst.replace(posdot, mach_end.length(), mach_end);
//		string mach_res = "/home/kristiina/Programs/mach.1.0.18.Linux/executables/First/mach1_"+ subst;
//		printf("Sample file %s reading mach file %s \n", sample_file.c_str(), mach_res.c_str());
//
//		HaplotypePair mach_haps = loadHaplotypesFromMACH(mach_res.c_str());
//
//		double mach_switch = mach_haps.switchError(true_haps);
//		double mach_geno = mach_haps.genotypeError(true_haps);
//
//		double mine_switch = mine_haps.switchError(true_haps);
//		double mine_geno = mine_haps.genotypeError(true_haps);
//
//		printf("File: %s \n \n", sample_file.c_str());
//		printf("Mine: %f    %f \n", mine_switch, mine_geno);
//		printf("Mach: %f    %f \n", mach_switch, mach_geno);
//
//
//		switch_error_mine.push_back(mine_switch);
//		geno_error_mine.push_back(mine_geno);
//		switch_error_mach.push_back(mach_switch);
//		geno_error_mach.push_back(mach_geno);
//
//
//	}
//	const char * out_file = "../../Data/1KGData/vcfs/chrom20/A/results_A_macherror0_5.csv";
//
//	writeVectorToCSV(out_file, switch_error_mine);
//	writeVectorToCSV(out_file, switch_error_mach);
//	writeVectorToCSV(out_file, geno_error_mine);
//	writeVectorToCSV(out_file, geno_error_mach);




//
//	HaplotypePhaser phaser;
////
////	const char * sample_name = "NA18909";
////	int sample_ind = VcfUtils::get_sample_index(sample_file, sample_name);
//
//	int sample_ind = 624;
//	phaser.LoadData(ref_file, sample_file, sample_ind);
////
////	VcfUtils::HaplotypePair true_haps = VcfUtils::loadHaplotypesFromVCF(sample_file, 624);
////
////
//////	VcfUtils::HaplotypePair beagle_haps = VcfUtils::loadHaplotypesFromVCF("/home/kristiina/Programs/beagle/4_A_snps_NA18909_beagle.vcf", 0);
////	VcfUtils::HaplotypePair beagle_haps = VcfUtils::loadHaplotypesFromVCF("/home/kristiina/Programs/beagle/4_A_snps_0.500000_0.001000_NA18909_beagle.vcf", 0);
////
////
////
//////	VcfUtils::HaplotypePair mach_haps =
//////			VcfUtils::loadHaplotypesFromMACH("/home/kristiina/Programs/mach.1.0.18.Linux/executables/errorRate0_1/mach1_4_A_snps_NA18909.out");
////	VcfUtils::HaplotypePair mach_haps =
////			VcfUtils::loadHaplotypesFromMACH("/home/kristiina/Programs/mach.1.0.18.Linux/executables/errorRate0_1/mach1_4_A_snps_0.500000_0.001000_NA18909.out");
////
////
//	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//	cout << "Starting Forward \n";
//	phaser.CalcScaledForward();
//
//	cout << "Starting Backward \n";
//	phaser.CalcScaledBackward();
//
//	cout << "Starting Posterior \n";
//	phaser.CalcPosterior();
//	chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//	cout << "Time difference = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
//
//
//	int * ml_states = new int[phaser.num_markers];
//	phaser.SampleHaplotypesNew(ml_states);
//
//	phaser.PrintHaplotypes(ml_states);
//	delete [] ml_states;

//
//
//
//	printf("True: \n");
//	true_haps.print();
//	printf("\n");
//	printf("Inferred: \n");
//	inferred_haps.print();
//	printf("\n");
//	printf("Mach: \n");
//	mach_haps.print();
//	printf("\n");
//	printf("Beagle: \n");
//	beagle_haps.print();
//	printf("\n");
//
//
//////	printf("lengths true %d %d \n", true_haps.h1.size(), true_haps.h2.size());
//////	printf("lengths inferred %d %d \n", inferred_haps.h1.size(), inferred_haps.h2.size());
//////	printf("lengths mach %d %d \n", mach_haps.h1.size(), mach_haps.h2.size());
//
//	printf("Inferred \n");
//	printf("[true], [inferred] equal: %d \n", inferred_haps.isEqual(true_haps));
//	printf("[true], [inferred] switch: %f \n", inferred_haps.switchError(true_haps));
//	printf("[true], [inferred] geno: %f \n", inferred_haps.genotypeError(true_haps));
//	printf("\nMach \n");
//
//	printf("[true], [mach] equal: %d \n", true_haps.isEqual(mach_haps));
//	printf("[true], [mach] switch: %f \n", mach_haps.switchError(true_haps));
//	printf("[true], [mach] geno: %f \n", mach_haps.genotypeError(true_haps));
//
//	printf("\nBeagle \n");
//
//	printf("[true], [beagle] equal: %d \n", true_haps.isEqual(beagle_haps));
//	printf("[true], [beagle] switch: %f \n", beagle_haps.switchError(true_haps));
//	printf("[true], [beagle] geno: %f \n", beagle_haps.genotypeError(true_haps));
//
//
//	printf("\n \n Inferred Geno Errors \n");
//	inferred_haps.printGenoErrors(true_haps);
//
//	printf("\n \n Mach Geno Errors \n");
//	mach_haps.printGenoErrors(true_haps);
//
//	printf("\n \n Beagle Geno Errors \n");
//	beagle_haps.printGenoErrors(true_haps);













//	const char * file4_lowcov1_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.100000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov2_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.200000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov3_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.300000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov4_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.400000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov5_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.500000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov6_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.600000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov7_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.700000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov8_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.800000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov9_beagle = "/home/kristiina/Programs/beagle/4_A_snps_0.900000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov10_beagle = "/home/kristiina/Programs/beagle/4_A_snps_NA18909_beagle.vcf";


//	String orig_sample_file = "../../Data/1KGData/vcfs/chrom20/B/4_B_snps.vcf";
//	String reff_file = "../../Data/1KGData/vcfs/chrom20/B/4_B_10inds_snps.vcf";
//
//	const char * sample_name = "NA18909";
//	int sample_ind = VcfUtils::get_sample_index(sample_file, sample_name);
//	phaser.LoadData(reff_file, orig_sample_file, sample_ind);
//
//
//
//	VcfUtils::HaplotypePair true_haps = VcfUtils::loadHaplotypesFromVCF(orig_sample_file.c_str(), 624);
//
//
//
//	const char * file4_lowcov1_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.100000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov2_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.200000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov3_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.300000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov4_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.400000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov5_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.500000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov6_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.600000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov7_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.700000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov8_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.800000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov9_beagle = "/home/kristiina/Programs/beagle/4_B_snps_0.900000_0.001000_NA18909_beagle.vcf";
//	const char * file4_lowcov10_beagle = "/home/kristiina/Programs/beagle/4_B_snps_NA18909_beagle.vcf";
//
//	vector<const char *> files {file4_lowcov1_beagle,
//								file4_lowcov2_beagle,
//								file4_lowcov3_beagle,
//								file4_lowcov4_beagle,
//								file4_lowcov5_beagle,
//								file4_lowcov6_beagle,
//								file4_lowcov7_beagle,
//								file4_lowcov8_beagle,
//								file4_lowcov9_beagle,
//								file4_lowcov10_beagle};
//
//
//
//
//	vector<double>  switch_error_beagle;
//	vector<double>  geno_error_beagle;
//
//
//
//	for (auto file : files) {
//
//
//		VcfUtils::HaplotypePair beagle_haps = VcfUtils::loadHaplotypesFromVCF(file, 0);
//		double geno = beagle_haps.genotypeError(true_haps);
//		double switcher = beagle_haps.switchError(true_haps);
//		switch_error_beagle.push_back(switcher);
//		geno_error_beagle.push_back(geno);
//
//
//	}
//
//	const char * out_file = "../../Data/1KGData/vcfs/chrom20/B/results_beagle.csv";
//
//	VcfUtils::writeVectorToCSV(out_file, switch_error_beagle);
//	VcfUtils::writeVectorToCSV(out_file, geno_error_beagle);













	//
////	true_haps.printDiffs(mach_haps);
//
//
//	for(int i = 0; i < true_haps.length; i++) {
//		if(i==165 || i==217 || i==223|| i==275 || i==288||i==301 || i==313||i==323 || i == 338 || i == 343|| i == 345  ) {
//			printf("%s ", true_haps.h1[i].c_str());
//
//		}
//	}
//	printf("\n");
//	for(int i = 0; i < true_haps.length; i++) {
//		if(i==165 || i==217 || i==223|| i==275 || i==288||i==301 || i==313||i==323 || i == 338 || i == 343|| i == 345  ) {
//			printf("%s ", true_haps.h2[i].c_str());
//
//		}
//	}
//	printf("\n");
//	printf("\n");
//	for(int i = 0; i < true_haps.length; i++) {
//		if(i==165 || i==217 || i==223|| i==275 || i==288||i==301 || i==313||i==323 || i == 338 || i == 343|| i == 345  ) {
//			printf("%s ", mach_haps.h1[i].c_str());
//
//		}
//	}
//	printf("\n");
//	for(int i = 0; i < true_haps.length; i++) {
//		if(i==165 || i==217 || i==223|| i==275 || i==288||i==301 || i==313||i==323 || i == 338 || i == 343|| i == 345  ) {
//			printf("%s ", mach_haps.h2[i].c_str());
//
//		}
//	}
//	printf("\n");
//
//	//_________________________________________________________________________________________
//
//	printf("\n");
//	for(int i = 0; i < true_haps.length; i++) {
//		if(i==165 || i==217 || i==223|| i==275 || i==288||i==301 || i==313||i==323 || i == 338 || i == 343|| i == 345  ) {
//			printf("%s ", inferred_haps.h1[i].c_str());
//
//		}
//	}
//	printf("\n");
//	for(int i = 0; i < true_haps.length; i++) {
//		if(i==165 || i==217 || i==223|| i==275 || i==288||i==301 || i==313||i==323 || i == 338 || i == 343|| i == 345  ) {
//			printf("%s ", inferred_haps.h2[i].c_str());
//
//		}
//	}
//	printf("\n");
//
//
//
//	printf("[true], [mach gens] genotype error: %f \n", mach_haps.genotypeError(true_haps));
//	printf("[true], [inferred gens] genotype error: %f \n", inferred_haps.genotypeError(true_haps));
//

//	printf("mach inferred equal: %d \n", mach_haps.isEqual(inferred_haps));





//	std::vector<String> t1{ "A", "T", "G", "A" };
//	std::vector<String> t2{ "T", "C", "C", "T" };


//	std::vector<String> i1{ "A", "T", "C", "T" };
//	std::vector<String> i2{ "T", "C", "G", "A" };

//
//	std::vector<String> i1{ "A", "T", "G", "A" };
//	std::vector<String> i2{ "T", "C", "C", "T" };

//	std::vector<String> i1{ "A", "T", "G", "A" };
//	std::vector<String> i2{ "R", "R", "C", "R" };

//
//	VcfUtils::HaplotypePair th(t1,t2);
//	VcfUtils::HaplotypePair ih(i1,i2);
////
//	double e1 = ih.genotypeError(th);
//
//	printf("genotype Error = %f \n", e1);

//	std::vector<int> pos1{0,1,2,3};
//	std::vector<int> pos2{0};
//	std::vector<int> pos3{2};
//	std::vector<int> pos4{3};


//	printf("equal = %d \n", th.isEqual(ih));
//	printf("switch = %f \n", th.switchError(ih));
//	printf("genotype error1 = %f \n", ih.genotypeError(th, pos1));
//	printf("genotype error2 = %f \n", ih.genotypeError(th, pos2));
//	printf("genotype error3 = %f \n", ih.genotypeError(th, pos3));
//	printf("genotype error4 = %f \n", ih.genotypeError(th, pos4));

	/*

	cout << "HAPLOTYPES" << "\n";

	for(int i = 0; i < phaser.ped.count; i++ ) {
		printf("person %d : ", i);
		for(int m = 0; m < phaser.ped.markerCount; m++){
			printf(" %d ",phaser.haplotypes[i*2][m]);
		}
		printf("----");
		for(int m = 0; m < phaser.ped.markerCount; m++){
			printf(" %d ",phaser.haplotypes[i*2+1][m]);
		}

		printf("\n");
	}


	for(int i = 0; i < phaser.ped.markerCount; i++ ) {
		printf("marker id %d unphase marker index %d \n", i, VcfUtils::unphased_marker_subset[i]);
	}


	printf("FORWARD  \n");
	for(int m = 0; m < phaser.num_markers; m++) {
		for(int s = 0; s < phaser.num_states; s++){
			printf(" %f \t ", phaser.forward[m][s]);
		}
		printf("\n");
	}

	printf("SCALED FORWARD  \n");
	for(int m = 0; m < phaser.num_markers; m++) {
		for(int s = 0; s < phaser.num_states; s++){
			printf(" %f \t ", phaser.s_forward[m][s]);
		}
		printf("\n");
	}

	printf("BACKWARD  \n");

	for(int m = 0; m < phaser.num_markers; m++) {
		for(int s = 0; s < phaser.num_states; s++){
			printf(" %f \t ", phaser.backward[m][s]);
		}
		printf("\n");
	}



	printf("SCALED BACKWARD  \n");

	for(int m = 0; m < phaser.num_markers; m++) {
		for(int s = 0; s < phaser.num_states; s++){
			printf(" %f \t ", phaser.s_backward[m][s]);
		}
		printf("\n");
	}

	printf("NORMALIZERS \n");
	for(int m = 0; m < phaser.num_markers; m++) {
		printf(" %f \n ", phaser.normalizers[m]);
	}




	// TODO make test cases out of these
	for(int m = 0; m < phaser.num_markers; m++) {
		float cm = 0.0;
		for(int j = 0; j < phaser.num_states; j++){
			cm += phaser.forward[m][j];
		}
		printf("cm = %f \n", cm);
		for(int s = 0; s < phaser.num_states; s++){
			printf("m=%d s=%d : Scaled forward = %f \t cm*forward = %f \n", m, s, phaser.s_forward[m][s], phaser.forward[m][s]/cm);

		}

		float cmv2 = 1.0;
		for(int j = 0; j <= m; j++){
			cmv2 *= phaser.normalizers[j];
		}
		printf("cmv2 = %f \n", cmv2);
		for(int s = 0; s < phaser.num_states; s++){
			printf("m=%d s=%d : Scaled forward = %f \t cmv2*forward = %f \n", m, s, phaser.s_forward[m][s], phaser.forward[m][s]/cmv2);

		}



		float dm = 1.0;
		for(int j = m; j < phaser.num_markers; j++){
			dm *= phaser.normalizers[j];
		}

		printf("dm = %f \n", dm);

		for(int s = 0; s < phaser.num_states; s++){
			printf("m=%d s=%d : Scaled backward = %f \t dm*backward = %f \n", m, s, phaser.s_backward[m][s], phaser.backward[m][s]/dm);

		}




	}

*/

	//	printf("\n\n");
	//	cout << "GENOTYPES" << "\n";
	//
	//	printf("marker \t");
	//	for(int m = 0; m < phaser.ped.markerCount; m++){
	//		printf(" %d \t \t", m);
	//	}
	//	printf("\n");
	//	for(int i = 0; i < phaser.ped.count; i++ ) {
	//		printf("person %d : ", i);
	//		for(int m = 0; m < phaser.ped.markerCount; m++){
	//			//phaser.genotypes[i][m*3] = phaser.genotypes[i][m*3] +1;
	//			//phaser.genotypes[i][m*3+1] = phaser.genotypes[i][m*3+1] +1;
	//			//phaser.genotypes[i][m*3+2] = phaser.genotypes[i][m*3+2] -1;
	//
	//			cout << (int) (unsigned char) phaser.genotypes[i][m*3] << " ";
	//			cout << (int) (unsigned char) phaser.genotypes[i][m*3+1] << " ";
	//			cout << (int) (unsigned char) phaser.genotypes[i][m*3+2] << " ";
	//
	//
	//			printf("\t");
	//		}
	//
	//		printf("\n");
	//	}
	//	printf("num inds = %d \n", phaser.num_inds);
	//
	//
	//	printf("0: %e \n", phaser.phred_probs[0]);
	//	printf("3: %e \n", phaser.phred_probs[3]);
	//	printf("12: %e \n", phaser.phred_probs[12]);
	//	printf("87: %e \n", phaser.phred_probs[87]);
	//	printf("255: %e \n", phaser.phred_probs[255]);
	return(0);
}
