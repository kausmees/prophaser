#include <chrono>
#include "HaplotypePhaser.h"
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
	string subset_id = "B";
	string individual = "NA12890";
	string ref_set = "CEU_25";



	double error = 0.001;
	vector<double> coverages = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
//	vector<double> coverages = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};

	//	vector<double> coverages = {0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
//	vector<double> coverages = {0.01, 0.5};
//	vector<double> coverages = {0.1};
//	vector<double> coverages = {};


	vector<string> parameters = {};

	for (auto c : coverages) {
		parameters.push_back("_" + to_string(c) + "_" + to_string(error));
	}

	parameters.push_back("");





	string file_name_in;
	string file_name_out;

//	String sample_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_snps.vcf";
//	int sample_ind = ???;

	string sample_file = string(dir) + subset_id+"/" + data_id +"_" + subset_id + "_" + individual + "_snps.vcf";
	int sample_ind = 0;

	string ref_file = string(dir) + subset_id+"/" + data_id +"_" + subset_id + "_" + ref_set + "_snps.vcf";

	HaplotypePhaser phaser;

	string distance_code;
	phaser.setDistanceCode(atoi(distance_code.c_str()));

	phaser.LoadReferenceData(ref_file.c_str(), sample_file.c_str(), sample_ind);

	HaplotypePair true_haps = loadHaplotypesFromVCF(sample_file.c_str(), sample_ind);


	//////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;
	for(auto par : parameters) {
		file_name_in = string(dir) + subset_id + "/" + data_id + "_" + subset_id + "_" + individual + "_snps" + par + ".vcf" ;

		////////////////////////////////////////
		distance_code = "2";
		phaser.setDistanceCode(atoi(distance_code.c_str()));
		printf("Distance code %d \n", phaser.distance_code);


		file_name_out ="./Results/" + subset_id+ "/"+ ref_set + "/" + data_id + "_" + subset_id + "_" + individual + "_snps" + par + ".phased_"+ distance_code;
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
		HaplotypePair mine_haps1 = phaser.PrintHaplotypesToFile(ml_states, file_name_out.c_str());
		HaplotypePair mine_haps2 = phaser.PrintReferenceHaplotypes(ml_states,(file_name_out+"_ref_haps").c_str());
		printf("Equal = %d \n", mine_haps1.isEqual(mine_haps2));

		end= std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;;

		///////////////////////////////////////////////

		distance_code = "3";
		phaser.setDistanceCode(atoi(distance_code.c_str()));
		printf("Distance code %d \n", phaser.distance_code);


		file_name_out ="./Results/" + subset_id+ "/"+ ref_set + "/" + data_id + "_" + subset_id + "_" + individual + "_snps" + par + ".phased_"+ distance_code;
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
		vector<vector<double>> stats2 = phaser.GetPosteriorStats((file_name_out+"_stats").c_str());
		int * ml_states2 = new int[phaser.num_markers];
		for (int i = 0; i < phaser.num_markers; i++) {
			ml_states2[i] = stats2[i][39];
		}
		HaplotypePair mine_haps3 = phaser.PrintHaplotypesToFile(ml_states, file_name_out.c_str());
		HaplotypePair mine_haps4 = phaser.PrintReferenceHaplotypes(ml_states,(file_name_out+"_ref_haps").c_str());
		printf("Equal = %d \n", mine_haps3.isEqual(mine_haps4));

		end= std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;;

	}
}
