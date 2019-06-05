#include <chrono>
#include "HaplotypePhaser.h"
#include "Parameters.h"
#include "Pedigree.h"


int main(int argc, char ** argv){

	String dir;
	String res_dir;
	String s_file;
	String r_file;
	String map_file;
	String algo;

	double Ne;
	double error;
	int n_iter;
	bool geno_mls;


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
	long_parameters.addInt("Iterations", &n_iter);
	long_parameters.addString("Algorithm", &algo);
	// use most likely states of other haplotype to calculate posterior genotype probability for current haplotype (false: use most likley allele)
	long_parameters.addBool("Geno_mls", &geno_mls);



	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
	parameter_list.Read(argc, argv);
	parameter_list.Status();

	// Defaults
	Ne = Ne ? Ne : 11418.0;
	error = error ? error : 0.001;
	n_iter = n_iter ? n_iter : 3;
	geno_mls = geno_mls ? geno_mls : false;


	dir = !dir.IsEmpty() ? dir : "./";
	res_dir = !res_dir.IsEmpty() ? res_dir : "./";
	algo = !algo.IsEmpty() ? algo : "fixed";


	if (!(!algo.Compare("fixed") || !algo.Compare("separate") || !algo.Compare("integrated"))) {
		printf("No such algorithm \n ");
		return -1;
	};

	string geno_mlss;
	if (geno_mls) {
		geno_mlss = "gmls";
	}else {
		geno_mlss = "gmla";
	}


	string a = algo.c_str();
	string suffix = ".ls."+a+"."+geno_mlss+"."+to_string(n_iter);


	HaplotypePhaser phaser;
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
	//
	printf("With: \nNe %f \n", phaser.Ne);
	printf("error %f \n", phaser.error);
	printf("algorithm %s \n", algo.c_str());

	printf("\nTemplate vcf: %s \n\n ", vcf_template.c_str());

	printf("----------------\n");
	printf("Writing to : %s \n\n", (result_file+suffix).c_str());


	////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	phaser.LoadReferenceData(r_file.c_str(), map_file);

	VcfFileReader reader;
	VcfHeader header_read;

	reader.open((sample_file).c_str(), header_read);
	int num_samples = header_read.getNumSamples();
	reader.close();

	// n_samples x n_markers array of most likely genotypes
	vector<vector<int>> ml_genotypes;

	// n_samples x n_markers array of most likely states for each chromosome
	vector<vector<int>> ml_states_h1;
	vector<vector<int>> ml_states_h2;

	// n_samples x n_markers array of most likely alleles for each chromosome
	vector<vector<int>> ml_alleles_h1;
	vector<vector<int>> ml_alleles_h2;


	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;

	begin1 = chrono::steady_clock::now();

	for(int sample = 0; sample< num_samples; sample++) {

		cout << "\n\n------------- Starting sample " << sample << "-------------"<< endl;

		phaser.LoadSampleData(sample_file.c_str(), sample);
		printf("Num inds in ped after load sample: %d \n ", phaser.ped.count);


		cout << "Haplotypes:\n" << phaser.haplotypes.block(0,0,4,10) << endl;

		phaser.CalcAlleleCounts();
		cout << "Allele counts:\n" << phaser.allele_counts.segment(0, 10) << endl;


		cout << "Sample GLs:\n" << endl;
		for (int m = 0; m < 4 ; m++) {
			printf("--- marker %d : %f %f %f \n" , m, phaser.sample_gls[m*3], phaser.sample_gls[m*3 +1 ], phaser.sample_gls[m*3+2]);
		}

		begin = chrono::steady_clock::now();
		cout << "Starting Forward haplotype 1\n";
		phaser.CalcScaledForwardMarginalized();

		end= std::chrono::steady_clock::now();

		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;

		begin = chrono::steady_clock::now();
		cout << "Starting Backward haplotype 1 \n";
		//	phaser.CalcScaledBackward();
		phaser.CalcScaledBackwardMarginalized();


		end= std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		begin = chrono::steady_clock::now();
		cout << "Starting Stats haplotype 1\n";

		vector<vector<double>> stats_marginalized = phaser.GetPosteriorStatsMarginalized((result_file+"_stats_h1").c_str(), false);
		cout << "Done Stats h1\n";




		vector<vector<double>> stats;

		for (int n = 1; n < n_iter; n++) {
			cout << "Starting iteration " << n << "\n";

			int hap = (n % 2) + 1;

			phaser.curr_chrom = hap;


			if (!algo.Compare("fixed")) {
				cout << "Starting Forward haplotype " << hap << "\n";
				phaser.CalcScaledForwardFixed();
				end= std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
				begin = chrono::steady_clock::now();
				cout << "Starting Backward haplotype " << hap << "\n";
				phaser.CalcScaledBackwardFixed();
			}

			if (!algo.Compare("separate")) {
				cout << "Starting Forward haplotype " << hap << "\n";
				phaser.CalcScaledForwardSeparate();
				end= std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
				begin = chrono::steady_clock::now();
				cout << "Starting Backward haplotype " << hap << "\n";
				phaser.CalcScaledBackwardSeparate();
			}

			if (!algo.Compare("integrated")) {
				cout << "Starting Forward haplotype " << hap << "\n";
				phaser.CalcScaledForwardIntegrated();
				end= std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
				begin = chrono::steady_clock::now();
				cout << "Starting Backward haplotype " << hap << "\n";
				phaser.CalcScaledBackwardIntegrated();
			}


			end= std::chrono::steady_clock::now();
			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
			begin = chrono::steady_clock::now();
			cout << "Starting Stats haplotype" << hap << "\n";
			stats = phaser.GetPosteriorStats((result_file+"_stats").c_str(), geno_mls, false);
			cout << "Done Stats \n";

		};


		// push back a vector for this sample
		ml_genotypes.push_back({});
		ml_states_h1.push_back({});
		ml_states_h2.push_back({});

		ml_alleles_h1.push_back({});
		ml_alleles_h2.push_back({});



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

			ml_states_h1[sample].push_back(phaser.ml_states_c1[m]);
			ml_states_h2[sample].push_back(phaser.ml_states_c2[m]);

			ml_alleles_h1[sample].push_back(phaser.ml_alleles_c1[m]);
			ml_alleles_h2[sample].push_back(phaser.ml_alleles_c2[m]);


		};

		end = std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;

	};

	phaser.PrintGenotypesToVCF(ml_genotypes, (result_file + suffix + ".genos").c_str(), sample_file.c_str(), vcf_template.c_str());
	phaser.PrintMLSHaplotypesToVCF(ml_states_h1, ml_states_h2, (result_file + suffix + ".phased.MLS" ).c_str(), sample_file.c_str(), vcf_template.c_str());
	phaser.PrintMLAHaplotypesToVCF(ml_alleles_h1, ml_alleles_h2, (result_file + suffix + ".phased.MLA").c_str(), sample_file.c_str(), vcf_template.c_str());


	end = std::chrono::steady_clock::now();
	cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

	return 0;
}

