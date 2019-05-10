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
	String res_dir;
	String subset_id;
	String individual;
	String s_file;
	String r_file;
	String ref_set;
	String map_file;
	String coverage;
	double Ne;
	double error;
	double seq_error;
	int split;
	int niter;

	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Input files");
	//	long_parameters.addString("directory", &dir);
	long_parameters.addString("results_directory", &res_dir);
	long_parameters.addString("sample_file", &s_file);
	long_parameters.addString("reference_file", &r_file);

	//	long_parameters.addString("subset", &subset_id);
	//	long_parameters.addString("individual", &individual);
	//	long_parameters.addString("reference", &ref_set);
	//	long_parameters.addString("cov", &coverage);
	//	long_parameters.addDouble("seq_error", &seq_error);
	long_parameters.addString("map", &map_file);
	long_parameters.addGroup("Parameters");
	long_parameters.addDouble("Ne", &Ne);
	long_parameters.addDouble("error", &error);
	long_parameters.addInt("split", &split);
	long_parameters.addInt("niter", &niter);



	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
	parameter_list.Read(argc, argv);
	parameter_list.Status();

	// Defaults
	Ne = Ne ? Ne : 11418.0;
	seq_error = seq_error ? seq_error : 0.001;
	error = error ? error : 0.001;

	dir = !dir.IsEmpty() ? dir : "../../Data/1KGData/vcfs/chrom20/";
	res_dir = !res_dir.IsEmpty() ? res_dir : "./Results/";
	map_file = !map_file.IsEmpty() ? map_file : "../../Data/1KGData/vcfs/chrom20/maps/5_snps_interpolated_HapMap2_map_20";
	subset_id = !subset_id.IsEmpty() ? subset_id : "B";
	individual = !individual.IsEmpty() ? individual : "NA12890";
	ref_set = !ref_set.IsEmpty() ? ref_set : "CEU_10";
	std::stringstream ss;
	ss << std::fixed << std::setprecision(3) << Ne << "_" << error;



	////////////////////////////////////////////////////////////////

	//	vector<double> posteriors;
	//	posteriors.resize(5);
	//
	//	posteriors[0] = 0.02;
	//	posteriors[1] = 0.02;
	//	posteriors[2] = 0.02;
	//	posteriors[3] = 0.02;
	//	posteriors[4] = 0.02;
	//
	//
	//	vector<size_t> res = sort_indexes(posteriors);
	//
	//	printf("Testing \n");
	//	for (int i=0;i<5;i++){
	//
	//		printf("%d : %f \n", res[i], posteriors[res[i]]);
	//	};


	//	return 0;

	///////////////////////////////////////////////////////////////

	//	ChromosomePair a(1, 2);
	//	ChromosomePair b(3, 4);


	// rows x cols
	//	MatrixXf L(6,5);
	//    VectorXf u(6);
	//    VectorXf coefs_p0(5);
	//    VectorXf coefs_p1(5);
	//
	//    float e = 0.01;
	//	int p=3;
	//	int q=1;
	//
	//	L << 0,pow((1-e),2),pow((1-e),2)-1,0,pow((1-e),2),
	//		   0,2*(1-e)*e-1, 2*(1-e)*e, 0,2*(1-e)*e,
	//		   0,pow(e,2),pow(e,2),0,(pow(e,2))-1,
	//		   (1-e)*e,0,0,(1-e)*e*2-1,0,
	//		   1-(pow(e,2))-pow(1-e,2),0,0,-2*(pow(1-e,2) + pow(e,2)),0,
	//		   2*p*q,pow(p,2) +pow(q,2),pow(p,2) +pow(q,2),4*p*q,pow(p,2) +pow(q,2);
	//
	//   u << 0, 0, 0 , 0 , 0, 1;
	//     cout << "Here is the matrix L:\n" << L << endl;
	//     cout << "Here is the vector u:\n" << u << endl;
	//
	//     VectorXf x(5);
	//     x = L.bdcSvd(ComputeThinU | ComputeThinV).solve(u);
	//     cout << "The solution is:\n" << x << endl;
	//
	// 	coefs_p0 << q, p, p , 2*q , p;
	// 	coefs_p1 << p, q, q , 2*p , q;
	//
	// 	cout << "P0: " << coefs_p0.dot(x) << endl;
	//
	//
	//	return 0;

	////////// Fixed parameters and settings ////////////////

	string data_id = "5";
	string filetype = "_s";

	//	coverage = "0.010000";


	string par;
	if(coverage.IsEmpty()) {
		par = "";
	}
	else{
		par = "_" + string(coverage.c_str()) + "_" + to_string(seq_error);
	}

	string splitstr;
	if(!split) {
		splitstr = "";
	}
	else{
		splitstr = ".split." + to_string(split);
	}

	HaplotypePhaserSym phaser;
	string distance_code = "sym_Ne"+ss.str();

	//	HaplotypePhaser phaser;
	//	string distance_code = "orig_Ne"+ss.str();

	//////////////////// Setup //////////////////////////////

	phaser.Ne = Ne;
	phaser.error = error;

	// corresponds to git branch
	//	string executable = "lsmlsmla_x"+to_string(niter);
	//	string executable = "ls_st_om_mls_x"+to_string(niter);
	string executable = "ls_st_omla_mls_x"+to_string(niter);


	//	string sample_file=string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(individual) + par+ string(filetype) +".vcf.gz";
	//	string true_file = string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(individual) + "_snps.vcf.gz";
	//	string ref_file =  string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(ref_set)  + "_snps" + splitstr + ".vcf.gz";

	//	string result_file =string(res_dir) + string(subset_id)+ "/"+ string(ref_set) + "/" + string(data_id) + "_" + string(subset_id) + "_" + string(individual) + par  +
	//			string(filetype)+ ".phased_"+ executable + distance_code + splitstr;



	string sample_file = string(s_file);
	string ref_file = string(r_file);
	string result_file = string(res_dir) + string(s_file) + ".imputed";


	printf("Phasing file : %s \n", s_file.c_str());
	printf("Reference file : %s \n", r_file.c_str());
	printf("Writing to : %s \n", result_file.c_str());
	//
	printf("With: \nNe %f \n", phaser.Ne);
	printf("error %f \n", phaser.error);

	////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	phaser.LoadReferenceData(r_file.c_str());

	VcfFileReader reader;
	VcfHeader header_read;

	reader.open((sample_file).c_str(), header_read);
	int num_samples = header_read.getNumSamples();
	reader.close();

	// n_samples x n_markers array of most likely genotypes
	vector<vector<int>> ml_genotypes;

//	ml_genotypes.push_back({});
//	ml_genotypes.push_back({});
//	ml_genotypes[0].push_back(0);
//	ml_genotypes[0].push_back(1);
//	ml_genotypes[0].push_back(2);
//	ml_genotypes[1].push_back(0);
//	ml_genotypes[1].push_back(1);
//	ml_genotypes[1].push_back(2);
//	phaser.PrintGenotypesToVCF(ml_genotypes, result_file.c_str(), sample_file.c_str());
//	return 0;

	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;

	begin1 = chrono::steady_clock::now();

	for(int sample = 0; sample< num_samples; sample++) {

		cout << "Starting sample " << sample << endl;

		// Ugly solution to reset the pedigree
		printf("before: %d \n ", phaser.ped.count);
		Pedigree ped;
		phaser.ped = ped;
		printf("after: %d \n ", phaser.ped.count);

		phaser.LoadReferenceData(r_file.c_str());
		phaser.LoadSampleData(s_file.c_str(), sample);

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

		//	phaser.CalcScaledForward();
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

		vector<vector<double>> stats = phaser.GetPosteriorStatsMarginalized((result_file+"_stats_h1").c_str());
		cout << "Done Stats \n";

		for (int i = 0; i < phaser.num_markers; i++) {
			phaser.ml_states_h1[i] = stats[i][39];
			if(stats[i][41] >= stats[i][42]) {
				phaser.ml_alleles_h1[i] = 0;
			}
			else {
				phaser.ml_alleles_h1[i] = 1;
			}
		}

		cout << "Done ml states h1 \n";

		//	FILE * mout = fopen((result_file+"_mlstates").c_str(), "w");
		//	for(int i = 0; i<phaser.num_markers; sample++) {
		//		fprintf(mout,"%d\n", ml_states[sample]);
		//	}

//		phaser.PrintReferenceHaplotypesMarginalized((result_file+"_ref_haps").c_str());


		vector<vector<double>> stats2;



		// lsmls_mla
		//	for (int n = 1; n < niter; n++) {
		//		cout << "Starting iteration " << n << "\n";
		//
		//		int hap = (n % 2) + 1;
		//
		//		phaser.curr_hap = hap;
		//		cout << "Starting Forward haplotype " << hap << "\n";
		//		phaser.CalcScaledForward();
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Backward haplotype " << hap << "\n";
		//		phaser.CalcScaledBackward();
		//
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Stats haplotype" << hap << "\n";
		//
		//		if(n==1) {
		//			stats2 = phaser.GetPosteriorStatsSpecial((result_file+"_stats_h2").c_str());
		//		}
		//		else{
		//			stats2 = phaser.GetPosteriorStats((result_file+"_stats_h2").c_str());
		//
		//		}
		//
		//		cout << "Done Stats \n";
		//		cout << "Done ml states h2 \n";
		//
		//	};


		//	// ls_st_mlsmla
		//	for (int n = 1; n < niter; n++) {
		//		cout << "Starting iteration " << n << "\n";
		//
		//		int hap = (n % 2) + 1;
		//
		//		phaser.curr_hap = hap;
		//		cout << "Starting Forward haplotype " << hap << "\n";
		//		phaser.CalcScaledForwardSeparate();
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Backward haplotype " << hap << "\n";
		//		phaser.CalcScaledBackwardSeparate();
		//
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Stats haplotype" << hap << "\n";
		//
		//		if(n==1) {
		//			stats2 = phaser.GetPosteriorStatsSpecial((result_file+"_stats_h2").c_str());
		//		}
		//		else{
		//			stats2 = phaser.GetPosteriorStats((result_file+"_stats_h2").c_str());
		//
		//		}
		//
		//		cout << "Done Stats \n";
		//		cout << "Done ml states h2 \n";
		//
		//	};

		//	// ls_st_mls
		//	for (int n = 1; n < niter; n++) {
		//		cout << "Starting iteration " << n << "\n";
		//
		//		int hap = (n % 2) + 1;
		//
		//		phaser.curr_hap = hap;
		//		cout << "Starting Forward haplotype " << hap << "\n";
		//		phaser.CalcScaledForwardSeparate();
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Backward haplotype " << hap << "\n";
		//		phaser.CalcScaledBackwardSeparate();
		//
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Stats haplotype" << hap << "\n";
		//
		//		stats2 = phaser.GetPosteriorStatsSpecial((result_file+"_stats_h2").c_str());
		//
		//		cout << "Done Stats \n";
		//		cout << "Done ml states h2 \n";
		//
		//	};



		// ls_st_mla
		//	for (int n = 1; n < niter; n++) {
		//		cout << "Starting iteration " << n << "\n";
		//
		//		int hap = (n % 2) + 1;
		//
		//		phaser.curr_hap = hap;
		//		cout << "Starting Forward haplotype " << hap << "\n";
		//		phaser.CalcScaledForwardSeparate();
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Backward haplotype " << hap << "\n";
		//		phaser.CalcScaledBackwardSeparate();
		//
		//		end= std::chrono::steady_clock::now();
		//		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//		begin = chrono::steady_clock::now();
		//		cout << "Starting Stats haplotype" << hap << "\n";
		//
		//		stats2 = phaser.GetPosteriorStats((result_file+"_stats_h2").c_str());
		//
		//		cout << "Done Stats \n";
		//		cout << "Done ml states h2 \n";
		//
		//	};



		//	// ls_stm1_mls
		//		for (int n = 1; n < niter; n++) {
		//			cout << "Starting iteration " << n << "\n";
		//
		//			int hap = (n % 2) + 1;
		//
		//			phaser.curr_hap = hap;
		//			cout << "Starting Forward haplotype " << hap << "\n";
		//			phaser.CalcScaledForwardIntegrated();
		//			end= std::chrono::steady_clock::now();
		//			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//			begin = chrono::steady_clock::now();
		//			cout << "Starting Backward haplotype " << hap << "\n";
		//			phaser.CalcScaledBackwardIntegrated();
		//
		//			end= std::chrono::steady_clock::now();
		//			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
		//			begin = chrono::steady_clock::now();
		//			cout << "Starting Stats haplotype" << hap << "\n";
		//
		//			stats2 = phaser.GetPosteriorStatsSpecial((result_file+"_stats_h2").c_str());
		//			cout << "Done Stats \n";
		//			cout << "Done ml states h2 \n";
		//
		//		};



		//	// ls_st_om_mls
		for (int n = 1; n < niter; n++) {
			cout << "Starting iteration " << n << "\n";

			int hap = (n % 2) + 1;

			phaser.curr_hap = hap;
			cout << "Starting Forward haplotype " << hap << "\n";
			phaser.CalcScaledForwardSeparate();
			end= std::chrono::steady_clock::now();
			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
			begin = chrono::steady_clock::now();
			cout << "Starting Backward haplotype " << hap << "\n";
			phaser.CalcScaledBackwardSeparate();

			end= std::chrono::steady_clock::now();
			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
			begin = chrono::steady_clock::now();
			cout << "Starting Stats haplotype" << hap << "\n";

			stats2 = phaser.GetPosteriorStatsSpecial((result_file+"_stats_h2").c_str());

			cout << "Done Stats \n";
			cout << "Done ml states h2 \n";

		};

		// push back a vector for this sample
		ml_genotypes.push_back({});

		for(int m = 0; m < phaser.num_markers; m++) {

			float max_geno_prob = 0.0;
			int max_geno_code;
			for(int i = 0; i < 3; i++) {
				if(stats2[m][41+i] > max_geno_prob) {
					max_geno_prob = stats[m][41+i];
					max_geno_code = i;
				};
			};
			//			printf("Max geno code: %d \n ", max_geno_code);
			ml_genotypes[sample].push_back(max_geno_code);

		};


		// This was when there was only 1 sample
		//		HaplotypePair mine_genos = phaser.PrintGenotypesToFile(stats2, (result_file+"_genos").c_str(), s_file.c_str());
		//		HaplotypePair mine_haps1 = phaser.PrintHaplotypesToFile(result_file.c_str(), s_file.c_str());

		end = std::chrono::steady_clock::now();
		cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;

	};


//	printf("ml genos size: %d %d \n ", ml_genotypes.size(), ml_genotypes[0].size());
//
//
//	for (int s = 0; s < ml_genotypes.size() ; s++) {
//		for (int m = 0; m < 4 ; m++) {
//			printf("sample: %d  marker: %d geno: %d \n" , s, m, ml_genotypes[s][m]);
//		}
//
//	}
	printf("wtf \n");
	phaser.PrintGenotypesToVCF(ml_genotypes, result_file.c_str(), sample_file.c_str());


	//
	//	phaser.curr_hap = 1;
	//	cout << "Starting Forward haplotype 1 again\n";
	//	phaser.CalcScaledForward();
	//	end= std::chrono::steady_clock::now();
	//	cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
	//	begin = chrono::steady_clock::now();
	//	cout << "Starting Backward haplotype 1 again \n";
	//	phaser.CalcScaledBackward();
	//
	//		end= std::chrono::steady_clock::now();
	//	cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
	//	begin = chrono::steady_clock::now();
	//	cout << "Starting Stats haplotype 1 again\n";
	//
	//	vector<vector<double>> stats3 = phaser.GetPosteriorStats((result_file+"_stats_h3").c_str());
	//	cout << "Done Stats \n";
	//	cout << "Done ml states h1 again \n";
	//
	//
	//	HaplotypePair mine_genos = phaser.PrintGenotypesToFile(stats3, (result_file+"_genos").c_str(), sample_file.c_str());
	//	HaplotypePair mine_haps1 = phaser.PrintHaplotypesToFile(result_file.c_str(), sample_file.c_str());

	end = std::chrono::steady_clock::now();
	cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

	return 0;
}

