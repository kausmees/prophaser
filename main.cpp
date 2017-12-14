#include <chrono>
#include "HaplotypePhaser.h"
#include "HaplotypePhaserSym.h"
#include "Parameters.h"
#include "Pedigree.h"
#include "GenoUtils.h"

// eg ./phase --Ne 30000 --individual NA12717 --reference CEU_10 --cov 0.900000


int main(int argc, char ** argv){

//	String sample_file;
//	String ref_file;
	String dir;
	String subset_id;
	String individual;
	String ref_set;
	String map_file;
	String coverage;
	double Ne;
	double error;
	double seq_error;

	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Input files");
	long_parameters.addString("directory", &dir);
	long_parameters.addString("subset", &subset_id);
	long_parameters.addString("individual", &individual);
	long_parameters.addString("reference", &ref_set);
	long_parameters.addString("cov", &coverage);
	long_parameters.addDouble("seq_error", &seq_error);
	long_parameters.addString("map", &map_file);
	long_parameters.addGroup("Parameters");
	long_parameters.addDouble("Ne", &Ne);
	long_parameters.addDouble("error", &error);


	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
	parameter_list.Read(argc, argv);
	parameter_list.Status();

	// Defaults
	Ne = Ne ? Ne : 11418.0;
	seq_error = seq_error ? seq_error : 0.001;
	error = error ? error : 0.01;
	dir = !dir.IsEmpty() ? dir : "../../Data/1KGData/vcfs/chrom20/";
	map_file = !map_file.IsEmpty() ? map_file : "../../Data/1KGData/vcfs/chrom20/maps/5_snps_interpolated_HapMap2_map_20";
	subset_id = !subset_id.IsEmpty() ? subset_id : "B";
	individual = !individual.IsEmpty() ? individual : "NA12890";
	ref_set = !ref_set.IsEmpty() ? ref_set : "CEU_10";
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2) << Ne;


//	ChromosomePair a(1, 2);
//	ChromosomePair b(3, 4);


//	return 0;

////////// Fixed parameters and settings ////////////////

	string data_id = "5";
	string filetype = "_s";

	string par;
	if(coverage.IsEmpty()) {
		par = "";
	}
	else{
		par = "_" + string(coverage.c_str()) + "_" + to_string(seq_error);
	}
	HaplotypePhaserSym phaser;
	string distance_code = "sym_Ne"+ss.str();

//	HaplotypePhaser phaser;
//	string distance_code = "orig_Ne"+ss.str();

//////////////////// Setup //////////////////////////////

	phaser.Ne = Ne;
	phaser.error = error;

	string sample_file=string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(individual) + par+ string(filetype) +".vcf.gz";
	string true_file = string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(individual) + "_snps.vcf.gz";
	string ref_file =  string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(ref_set)  + "_snps.vcf.gz";

	string result_file ="./Results/" + string(subset_id)+ "/"+ string(ref_set) + "/" + string(data_id) + "_" + string(subset_id) + "_" + string(individual) + par  +
			string(filetype)+ ".phased_"+ distance_code;

	printf("Phasing file : %s \n", sample_file.c_str());
	printf("Reference file : %s \n", ref_file.c_str());
	printf("Writing to : %s \n", result_file.c_str());
//
	printf("With: \nNe %f \n", phaser.Ne);
//	printf("error %f \n", phaser.error);

	////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	phaser.LoadData(ref_file.c_str(), sample_file.c_str(), 0, map_file.c_str());

//	return 0;

	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;


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

	vector<vector<double>> stats = phaser.GetPosteriorStats((result_file+"_stats").c_str());
	cout << "Done Stats \n";
	int * ml_states = new int[phaser.num_markers];
	for (int i = 0; i < phaser.num_markers; i++) {
		ml_states[i] = stats[i][39];
	}
	cout << "Done ml states \n";

//	FILE * mout = fopen((result_file+"_mlstates").c_str(), "w");
//	for(int i = 0; i<phaser.num_markers; i++) {
//		fprintf(mout,"%d\n", ml_states[i]);
//	}

	HaplotypePair mine_haps1 = phaser.PrintHaplotypesToFile(ml_states, result_file.c_str());
	HaplotypePair mine_haps2 = phaser.PrintReferenceHaplotypes(ml_states,(result_file+"_ref_haps").c_str());
	HaplotypePair mine_genos = phaser.PrintGenotypesToFile(stats, (result_file+"_genos").c_str());

	cout << "Done Haps \n";

	if(!mine_haps1.isEqual(mine_haps2)){
		printf("Haps not equal \n");
	}
	end= std::chrono::steady_clock::now();
	cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
	cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

}

