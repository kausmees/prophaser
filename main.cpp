#include <chrono>
#include "HaplotypePhaser.h"
#include "HaplotypePhaserSym.h"
#include "Parameters.h"
#include "Pedigree.h"
#include "GenoUtils.h"


int main(int argc, char ** argv){
//	String sample_file;
//	String ref_file;
	String dir;


	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Input files");
//	long_parameters.addString("sampleVCF", &sample_file);
//	long_parameters.addString("referenceVCF", &ref_file);
	long_parameters.addString("directory", &dir);


	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
	parameter_list.Read(argc, argv);
	parameter_list.Status();


//	String dir = "../../Data/1KGData/vcfs/chrom20/";
	string data_id = "4";
	string individual = "NA12890";
//	string individual = "NA12717";
//	string individual = "NA12812";



//	string subset_id = "C";
	string subset_id = "B";
//	string subset_id = "A";


//	string ref_set = "CEU_10";
//	string ref_set = "CEU_25";
//	string ref_set = "CEU_50";
	string ref_set = "CEU_ALL_"+individual;

//	string ref_set = "YRI_ALL";





	double error = 0.001;
//	vector<double> coverages = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

// new standard
//	vector<double> coverages = {0.0, 0.02, 0.05, 0.08, 0.1, 0.4, 0.7, 1.0};
	vector<double> coverages = {1.0};

//	vector<double> coverages = {};

	vector<string> parameters = {};
	for (auto c : coverages) {
		parameters.push_back("_" + to_string(c) + "_" + to_string(error));
	}

	parameters.push_back("");


	string fileend = "_s.vcf";
	string filetype = "_s";


	string file_name_in;
	string file_name_out;

//	String sample_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_snps.vcf";
//	int sample_ind = ???;

	int sample_ind = 0;
	string sample_file = string(dir) + subset_id+"/" + data_id +"_" + subset_id + "_" + individual + fileend;
	string true_file = string(dir) + subset_id+"/" + data_id +"_" + subset_id + "_" + individual + "_snps.vcf";
	string ref_file = string(dir) + subset_id+"/" + data_id +"_" + subset_id + "_" + ref_set + fileend;

	string distance_code;

//
//	HaplotypePhaserSym phaser;
////	 distance code new_sym_t20 has Ne 11418
//	distance_code = "new_sym_t20_Ne10x6";
//    phaser.Ne = 1000000.0;


//	HaplotypePhaserSym phaser;
////	 distance code new_sym_t20 has Ne 11418
//	distance_code = "new_sym_t20_Ne11418";
//    phaser.Ne = 11418.0;



//	HaplotypePhaser phaser;
////	 distance code 6_new has Ne 11418
//	distance_code = "6_Ne10x6";
//	phaser.Ne = 1000000.0;

	HaplotypePhaser phaser;
//	 distance code 6_new has Ne 11418
	distance_code = "6_Ne11418";
	phaser.Ne = 11418.0;




//	phaser.setDistanceCode(atoi(distance_code.c_str()));
//
	phaser.LoadReferenceData(ref_file.c_str(), sample_file.c_str(), sample_ind);

	printf("Loading true haps from %s \n", true_file.c_str());
	HaplotypePair true_haps = loadHaplotypesFromVCF(true_file.c_str(), sample_ind);


	//////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////
	printf("start \n");

	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;

	for(auto par : parameters) {
		file_name_in = string(dir) + subset_id + "/" + data_id + "_" + subset_id + "_" + individual + filetype + par + ".vcf" ;

		printf("Phasing file : %s \n", file_name_in.c_str());
		printf("Reference file : %s \n", ref_file.c_str());


		file_name_out ="./Results/" + subset_id+ "/"+ ref_set + "/" + data_id + "_" + subset_id + "_" + individual + filetype + par + ".phased_"+ distance_code;
		printf("Writing to to file : %s \n", file_name_out.c_str());

		phaser.LoadSampleData(ref_file.c_str(), file_name_in.c_str(), 0);
		begin1 = chrono::steady_clock::now();
		begin = chrono::steady_clock::now();
		cout << "Starting Forward \n";
		phaser.CalcScaledForward();
		end= std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		begin = chrono::steady_clock::now();
		cout << "Starting Backward \n";
		phaser.CalcScaledBackward();
		end= std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		begin = chrono::steady_clock::now();
		cout << "Starting Stats \n";

		vector<vector<double>> stats = phaser.GetPosteriorStats((file_name_out+"_stats").c_str());
		int * ml_states = new int[phaser.num_markers];
		for (int i = 0; i < phaser.num_markers; i++) {
			ml_states[i] = stats[i][39];
		}

		FILE * mout = fopen((file_name_out+"_mlstates").c_str(), "w");

		for(int i = 0; i<phaser.num_markers; i++) {
			fprintf(mout,"%d\n", ml_states[i]);
		}

		HaplotypePair mine_haps1 = phaser.PrintHaplotypesToFile(ml_states, file_name_out.c_str());
		mine_haps1.print();
		HaplotypePair mine_haps2 = phaser.PrintReferenceHaplotypes(ml_states,(file_name_out+"_ref_haps").c_str());
		printf("Equal haps= %d \n", mine_haps1.isEqual(mine_haps2));

		HaplotypePair mine_genos = phaser.PrintGenotypesToFile(stats, (file_name_out+"_genos").c_str());
		printf("Equal hap geno= %d \n", mine_haps1.isEqual(mine_genos));

		end= std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

	}

}
