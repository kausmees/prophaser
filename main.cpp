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
	String ref_set;
	String map_file;
	String coverage;
	double Ne;
	double error;
	double seq_error;
	int split;

	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Input files");
	long_parameters.addString("directory", &dir);
	long_parameters.addString("results_directory", &res_dir);
	long_parameters.addString("subset", &subset_id);
	long_parameters.addString("individual", &individual);
	long_parameters.addString("reference", &ref_set);
	long_parameters.addString("cov", &coverage);
	long_parameters.addDouble("seq_error", &seq_error);
	long_parameters.addString("map", &map_file);
	long_parameters.addGroup("Parameters");
	long_parameters.addDouble("Ne", &Ne);
	long_parameters.addDouble("error", &error);
	long_parameters.addInt("split", &split);


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
	string executable = "ls_";


	string sample_file=string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(individual) + par+ string(filetype) +".vcf.gz";
	string true_file = string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(individual) + "_snps.vcf.gz";
	string ref_file =  string(dir) + string(subset_id)+"/" + string(data_id) +"_" + string(subset_id) + "_" + string(ref_set)  + "_snps" + splitstr + ".vcf.gz";

	string result_file =string(res_dir) + string(subset_id)+ "/"+ string(ref_set) + "/" + string(data_id) + "_" + string(subset_id) + "_" + string(individual) + par  +
			string(filetype)+ ".phased_"+ executable +distance_code + splitstr;

	printf("Phasing file : %s \n", sample_file.c_str());
	printf("Reference file : %s \n", ref_file.c_str());
	printf("Writing to : %s \n", result_file.c_str());
//
	printf("With: \nNe %f \n", phaser.Ne);
	printf("error %f \n", phaser.error);

	////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

//	ChromosomePair cp (1,2);
//	int s = sizeof(cp);
//	printf("Size of chrom pair: % d \n" ,s);
//
//	return 0;

	phaser.LoadData(ref_file.c_str(), sample_file.c_str(), 0, map_file.c_str());

	cout << "Haplotypes:\n" << phaser.haplotypes.block(0,0,22,10) << endl;
	phaser.CalcAlleleCounts();
	cout << "Allele counts:\n" << phaser.allele_counts.segment(0, 10) << endl;

	return 0;

	chrono::steady_clock::time_point begin1;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;


	begin1 = chrono::steady_clock::now();
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

	vector<vector<double>> stats = phaser.GetPosteriorStats((result_file+"_stats").c_str());
	cout << "Done Stats \n";
	int * ml_states_h1 = new int[phaser.num_markers];
	for (int i = 0; i < phaser.num_markers; i++) {
		ml_states_h1[i] = stats[i][39];
	}
	cout << "Done ml states h1 \n";









//	FILE * mout = fopen((result_file+"_mlstates").c_str(), "w");
//	for(int i = 0; i<phaser.num_markers; i++) {
//		fprintf(mout,"%d\n", ml_states[i]);
//	}

	HaplotypePair mine_haps1 = phaser.PrintHaplotypesToFile(ml_states_h1, result_file.c_str(), sample_file.c_str());
	phaser.PrintReferenceHaplotypes(ml_states_h1,(result_file+"_ref_haps").c_str());
	HaplotypePair mine_genos = phaser.PrintGenotypesToFile(stats, (result_file+"_genos").c_str(), sample_file.c_str());

//	mine_haps1.print();

	cout << "Done Haps \n";

//	if(!mine_haps1.isEqual(mine_haps_byref)){
//		printf("Haps not equal \n");
//	}

	end= std::chrono::steady_clock::now();
	cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
	cout << "Total: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin1).count() << " millisec" <<endl<<endl<<endl;

}

