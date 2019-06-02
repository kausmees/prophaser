#include <chrono>
//#include "HaplotypePhaserSym.h"
#include "HaplotypePhaser.h"
#include "Parameters.h"
#include "Pedigree.h"


int main(int argc, char ** argv){

	String dir;
	String res_dir;
	String s_file;
	String r_file;
	String map_file;
	double Ne;
	double error;

	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Files");
	long_parameters.addString("directory", &dir);
	long_parameters.addString("results_directory", &res_dir);
	long_parameters.addString("sample_file", &s_file);
	long_parameters.addString("reference_file", &r_file);
	long_parameters.addString("map", &map_file);

	long_parameters.addGroup("Parameters");
	long_parameters.addDouble("Ne", &Ne);
	long_parameters.addDouble("Error", &error);


	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
	parameter_list.Read(argc, argv);
	parameter_list.Status();

	// Defaults
	Ne = Ne ? Ne : 11418.0;
	error = error ? error : 0.001;

	dir = !dir.IsEmpty() ? dir : "../../Data/1KGData/vcfs/chrom20/";
	res_dir = !res_dir.IsEmpty() ? res_dir : "./Results/";
	map_file = !map_file.IsEmpty() ? map_file : "/home/kristiina/Projects/Data/1KGData/vcfs/chrom20/maps/5_snps_interpolated_HapMap2_map_20";


//	HaplotypePhaserSym phaser;
	HaplotypePhaser phaser;


// /////
//	vector<ChromosomePair> states;
//
//	for(int i = 0; i < 5; i++) {
//		for(int j = 0; j < 5; j++) {
//			states.push_back(ChromosomePair(i,j));
//		}
//	}
//
//	printf("Num in vec: %d \n", states.size());
//
//	for (int s1 = 0; s1 < states.size(); s1++) {
//		for (int s2 = 0; s2 < states.size(); s2++) {
//
//			int eq = states[s1].TransitionCase(states[s2]);
//			printf("%d %d | %d %d ==> %d\n", states[s1].first, states[s1].second, states[s2].first, states[s2].second, eq);
//		}
//
//	}
//
//		return 0;
//	////




	phaser.Ne = Ne;
	phaser.error = error;


	string sample_file = string(dir) + string(s_file);
	string sfilebase = string(s_file).substr(0, s_file.Find(".vcf.gz"));
	string ref_file = string(r_file);
	string result_file = string(res_dir) + sfilebase;

	string vcf_template = string(dir) + sfilebase + ".template.vcf";

	printf("Phasing file : %s \n", sample_file.c_str());
	printf("----------------\n");

	printf("Reference file : %s \n", ref_file.c_str());
	printf("Map file : %s \n", map_file.c_str());
	printf("Writing to : %s \n\n", result_file.c_str());
	//
	printf("With: \nNe %f \n", phaser.Ne);
	printf("error %f \n", phaser.error);


	printf("\nTemplate vcf: %s \n\n ", vcf_template.c_str());


	////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	phaser.LoadReferenceData(r_file.c_str(), map_file.c_str());

	VcfFileReader reader;
	VcfHeader header_read;

	reader.open((sample_file).c_str(), header_read);
	int num_samples = header_read.getNumSamples();
	reader.close();

	// n_samples x n_markers array of most likely genotypes
	vector<vector<int>> ml_genotypes;

	// n_samples x n_markers array of most likely states
	vector<vector<int>> ml_states;


	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;

	begin1 = chrono::steady_clock::now();

	for(int sample = 0; sample< num_samples; sample++) {

		cout << "\n\n------------- Starting sample " << sample << "-------------"<< endl;

		phaser.LoadSampleData(sample_file.c_str(), sample);
		printf("Num inds in ped after load sample: %d \n ", phaser.ped.count);


		cout << "Haplotypes:\n" << phaser.haplotypes.block(0,0,4,10) << endl;

		cout << "Sample GLs:\n" << endl;
		for (int m = 0; m < 4 ; m++) {
			printf("--- marker %d : %f %f %f \n" , m, phaser.sample_gls[m*3], phaser.sample_gls[m*3 +1 ], phaser.sample_gls[m*3+2]);
		}

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

		vector<vector<double>> stats = phaser.GetPosteriorStats((result_file+"_"+to_string(sample)+"_stats").c_str(), true);
		cout << "Done Stats \n";

		// push back a vector for this sample
		ml_genotypes.push_back({});
		ml_states.push_back({});

		for(int m = 0; m < phaser.num_markers; m++) {
			float max_geno_prob = 0.0;
			int max_geno_code;
			for(int i = 0; i < 3; i++) {
				if(stats[m][41+i] > max_geno_prob) {
					max_geno_prob = stats[m][41+i];
					max_geno_code = i;
				};
			};
			ml_genotypes[sample].push_back(max_geno_code);
			ml_states[sample].push_back(stats[m][39]);
		};

		end = std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;

	};


//	phaser.PrintGenotypesToVCF(ml_genotypes, (result_file + ".genos").c_str(), sample_file.c_str(), vcf_template.c_str());
//	phaser.PrintHaplotypesToVCF(ml_states, (result_file + ".phased").c_str(), sample_file.c_str(), vcf_template.c_str());

	phaser.PrintGenotypesToVCF(ml_genotypes, (result_file + ".genos.full").c_str(), sample_file.c_str(), vcf_template.c_str());
	phaser.PrintHaplotypesToVCF(ml_states, (result_file + ".phased.full").c_str(), sample_file.c_str(), vcf_template.c_str());



//	for(int sample = 0; sample < num_samples; sample ++) {
//		printf("sample %d \n", sample);
//		for(int marker = 0; marker < 10; marker ++) {
//			printf("%d ", ml_states[sample][marker]);
//		}
//		printf("\n");
//	}


	end = std::chrono::steady_clock::now();
	cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

	return 0;
}

