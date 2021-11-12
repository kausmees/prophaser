#include <chrono>
//#include "HaplotypePhaserSym.h"
#include "HaplotypePhaser.h"
#include "Parameters.h"
#include "Pedigree.h"



int main(int argc, char ** argv){

	String dir;
	String res_file;
	String s_file;
	String r_file;
	String map_file;

	double Ne;
	double error;
	int niter;

	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Files");
	long_parameters.addString("directory", &dir);
	long_parameters.addString("outfile", &res_file);
	long_parameters.addString("sample_file", &s_file);
	long_parameters.addString("reference_file", &r_file);
	long_parameters.addString("map", &map_file);

	long_parameters.addGroup("Parameters");
	long_parameters.addDouble("Ne", &Ne);
	long_parameters.addDouble("Error", &error);
	long_parameters.addInt("iterations", &niter);

	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));
	parameter_list.Read(argc, argv);
	parameter_list.Status();

	// Defaults
	Ne = Ne ? Ne : 11418.0;
	error = error ? error : 1e-12;

	dir = !dir.IsEmpty() ? dir : "../../Data/1KGData/vcfs/chrom20/";
	res_file = !res_file.IsEmpty() ? res_file : "./Results/" + s_file;


// TODO implement base class and have algorithm choice in parameters

	HaplotypePhaser phaser;


	phaser.Ne = Ne;
	phaser.error = error;


	string sample_file = string(dir) + string(s_file);
	string sfilebase = string(s_file).substr(0, s_file.Find(".vcf.gz"));
	string ref_file = string(r_file);
	string result_file = string(res_file);

	string vcf_template = string(dir) + sfilebase + ".template.vcf";
    string vcf_template_gtgp = string(dir) + sfilebase + ".template_gtgp.vcf";

	printf("Phasing file : %s \n", sample_file.c_str());
	printf("----------------\n");

	printf("Reference file : %s \n", ref_file.c_str());
	printf("Map file : %s \n", map_file.c_str());
	//
	printf("With: \nNe %f \n", phaser.Ne);
	printf("error %.13e \n", phaser.error);

	printf("\nTemplate GT vcf: %s \n\n ", vcf_template.c_str());
    printf("\nTemplate GT:GP vcf: %s \n\n ", vcf_template_gtgp.c_str());

	printf("----------------\n");
	printf("Writing to : %s \n\n", (result_file).c_str());


	////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	phaser.LoadData(r_file.c_str(), sample_file.c_str(), map_file);

	VcfFileReader reader;
	VcfHeader header_read;


	// n_samples x n_markers array of most likely genotypes
	vector<vector<int>> ml_genotypes;

    // n_samples x n_markers x 3 array of most posterior genotypes probabilities
    vector<vector<vector<phaserreal>>> ml_postprobs;

	// n_samples x n_markers array of most likely states
	vector<vector<int>> ml_states;

	// n_markers array of most likely state for the current sample
	vector<int> final_ml_states_this_sample(phaser.num_markers, -1);

	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;

	begin1 = chrono::steady_clock::now();

	// first iteration: we ignore the sample haplotypes, they havent been initalized to anything good yet (just 0s)
	phaser.SetNumHaps(phaser.num_ref_inds * 2);

	for(int iter = 0; iter < niter; iter++){
		cout << "\n\n------------- Starting iteration " << iter << "-------------"<< endl;

		for(int sample = 0; sample < phaser.num_sample_inds; sample++) {

			phaser.current_sample = sample;
			phaser.current_sample_hap_index = (phaser.num_ref_inds + sample) * 2;


			cout << "\n\n------------- Starting sample " << sample << "-------------"<< endl;
			//(Index startRow, Index startCol, Index blockRows, Index blockCols)
			cout << "Haplotypes:\n" << phaser.haplotypes.block(0,0,4,10) << endl;
			cout << "Sample GLs:\n" << endl;

			for (int m = 0; m < 4 ; m++) {
				printf("--- marker %d : %f %f %f \n" , m, phaser.sample_gls(sample,m*3), phaser.sample_gls(sample,m*3 +1), phaser.sample_gls(sample,m*3+2));
			}
			for (int m = phaser.num_markers-5; m < phaser.num_markers-1 ; m++) {
				printf("--- marker %d : %f %f %f \n" , m, phaser.sample_gls(sample,m*3), phaser.sample_gls(sample,m*3 +1), phaser.sample_gls(sample,m*3+2));
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Old code
//			cout << "Starting Stats \n";
//			vector<vector<phaserreal>> stats = phaser.GetPosteriorStats((result_file+"_"+to_string(sample)+"_stats").c_str(), true);
//			cout << "Done Stats \n";
//
//			// push back a vector for this sample
//			ml_genotypes.push_back({});
//			ml_states.push_back({});
//			ml_postprobs.push_back({});
//
//			for(int m = 0; m < phaser.num_markers; m++) {
//				float max_geno_prob = 0.0;
//				int max_geno_code;
//				ml_postprobs[sample].push_back({});
//				ml_postprobs[sample][m].resize(3);
//				for(int i = 0; i < 3; i++) {
//					if(stats[m][41+i] > max_geno_prob) {
//						max_geno_prob = stats[m][41+i];
//						max_geno_code = i;
//					};
//				};
//				ml_genotypes[sample].push_back(max_geno_code);
//
//				for(int g = 0; g < 3; g++) {
//					ml_postprobs[sample][m][g] = stats[m][41+g];
//				}
//	//            if (m % 1000 == 0) {
//	//                printf("\nSample %d", sample);
//	//                printf("\nposterior probabilities RR,RA,AA at marker %d = %f,%f,%f", m, ml_postprobs[sample][m][0], ml_postprobs[sample][m][1], ml_postprobs[sample][m][2]);
//	//                printf("\n");
//	//            }
//
//				ml_states[sample].push_back(stats[m][39]);
//			};
//
//			end = std::chrono::steady_clock::now();
//			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////7

			// if its not the last iteration, we just get the ML states - no ml genos
			if (iter < niter -1) {

				cout << "Starting Get ML States \n";

				vector<int> ml_states_this_sample = phaser.GetMLStates();

				end = std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec    size: " <<  ml_states_this_sample.size() << endl<<endl;
				///////////////////////////////////////////////////////////////////////////////
				cout << "ml states: " << endl << endl;
				for (int m = 0; m < 10 ; m++) {
					cout << ml_states_this_sample[m] << endl;
				}
				begin = chrono::steady_clock::now();
				cout << "Starting updating haplotypes" << endl << endl;

				// Update the current estimate of this samples haplotypes according to the most likley states
				phaser.UpdateSampleHaplotypes(ml_states_this_sample);
				end = std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl<<endl;
			}

			// final iteration: we also get the ML genos
			// and store ml_states and ml_genos for all samples for printing to file
			else {
				///////////////////////////////////////////////////////////////////////////////////////////
				cout << "Starting Stats \n";
				vector<vector<phaserreal>> stats = phaser.GetPosteriorStats((result_file+"_"+to_string(sample)+"_stats").c_str(), iter == niter -1);

				end = std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<  endl<<endl;


				begin = chrono::steady_clock::now();
				cout << "Starting updating states \n";

				// push back a vector for this sample
				ml_genotypes.push_back({});
				ml_postprobs.push_back({});


				for(int m = 0; m < phaser.num_markers; m++) {
					float max_geno_prob = 0.0;
					int max_geno_code;
		            ml_postprobs[sample].push_back({});
		            ml_postprobs[sample][m].resize(3);
					for(int i = 0; i < 3; i++) {
						if(stats[m][41+i] > max_geno_prob) {
							max_geno_prob = stats[m][41+i];
							max_geno_code = i;
						};
					};
					ml_genotypes[sample].push_back(max_geno_code);

					for(int g = 0; g < 3; g++) {
						ml_postprobs[sample][m][g] = stats[m][41+g];
					}
					final_ml_states_this_sample[m] = stats[m][39];
				};
				ml_states.push_back(final_ml_states_this_sample);

				end = std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl<<endl;

				begin = chrono::steady_clock::now();
				cout << "Starting updating haplotypes \n";

				// Update the current estimate of this samples haplotypes according to the most likley states
				phaser.UpdateSampleHaplotypes(final_ml_states_this_sample);
				end = std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl<<endl;

				///////////////////////////////////////////////////////////////////////////////////////////


				end = std::chrono::steady_clock::now();
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl<<endl;

			}



		}// end sample loop


		if (iter == 0) {

			begin = chrono::steady_clock::now();
			cout << "Starting set num haplotypes \n";

			// Number of reference haplotypes that will be considered when handling one sample
			// After first iteration, sample haps have been filled and will be used for other samples.
			phaser.SetNumHaps((phaser.num_inds - 1) * 2);
			end = std::chrono::steady_clock::now();
			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl<<endl;

		}
	};


	phaser.PrintGenotypesToVCF(ml_genotypes, (result_file + ".genos" ).c_str(), sample_file.c_str(), vcf_template.c_str());
    phaser.PrintPostGenotypesToVCF(ml_genotypes, ml_postprobs, (result_file + ".postgenos" ).c_str(), sample_file.c_str(), vcf_template_gtgp.c_str());

	end = std::chrono::steady_clock::now();
	cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

	return 0;
}

