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


	String dir = "../../Data/1KGData/vcfs/chrom20/";
	String data_id = "4";


//	vector<double> coverages = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
//	vector<double> coverages = {0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

	//	vector<double> coverages = {0.1};
	vector<double> coverages = {};


	String subset_id = "B";
	double error = 0.001;
	String individual = "NA12890";
	String ref_set = "CEU_50";
	String distance_code = "4";



	String file_name_in;
	String file_name_out;

//	String sample_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_snps.vcf";
//	int sample_ind = ???;

	String sample_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_" + individual + "_snps.vcf";
	int sample_ind = 0;

	String ref_file = dir+subset_id+"/" + data_id +"_" + subset_id + "_" + ref_set + "_snps.vcf";

	HaplotypePhaser phaser;
	phaser.LoadReferenceData(ref_file, sample_file, sample_ind);

	HaplotypePair true_haps = loadHaplotypesFromVCF(sample_file, sample_ind);


	//////////////////////////////////////////// phasing start //////////////////////////////////////////////////////////////

	file_name_in = dir + subset_id + "/" + data_id + "_" + subset_id + "_" + individual + "_snps.vcf" ;
	file_name_out = "./Results/" + subset_id+ "/"+ ref_set + "/" + data_id + "_" + subset_id + "_" + individual + "_snps.phased_"+ distance_code;

	phaser.LoadSampleData(ref_file, file_name_in, 0);
	printf("Num states = %d Num markers = %d Num inds = %d \n", phaser.num_states, phaser.num_markers, phaser.num_inds);

//	for(int i=0; i<phaser.num_markers;i++){
//		printf("Marker %d : %f %f %f \n", i, phaser.sample_gls[3*i],phaser.sample_gls[3*i+1],phaser.sample_gls[3*i+2]);
//
//	}


	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
	cout << "Starting Forward \n";
	phaser.CalcScaledForward();

	cout << "Starting Backward \n";
	phaser.CalcScaledBackward();


	cout << "Starting Posterior \n";
	int * ml_states = new int[phaser.num_markers];
	phaser.GetMLHaplotypes(ml_states);
	cout << "DOMNE \n";


	HaplotypePair mine_haps = phaser.PrintHaplotypesToFile(ml_states, file_name_out);
	HaplotypePair mine_haps2 = phaser.PrintReferenceHaplotypes(ml_states,file_name_out+"_ref_haps");
	mine_haps.print();
	printf("Equal = %d \n", mine_haps.isEqual(mine_haps2));
	delete [] ml_states;

	chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

	cout << "Time original = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;

	for(auto cov : coverages) {
		file_name_in = dir + subset_id + "/" + data_id + "_" + subset_id + "_" + individual + "_snps_" + (to_string(cov)).c_str() + "_" + (to_string(error)).c_str() + ".vcf" ;
		file_name_out ="./Results/" + subset_id+ "/"+ ref_set + "/" + data_id + "_" + subset_id + "_" + individual + "_snps_" + (to_string(cov)).c_str() + "_" + (to_string(error)).c_str() + ".phased_"+ distance_code;

		phaser.LoadSampleData(ref_file, file_name_in , 0);

		begin = chrono::steady_clock::now();
		cout << "Starting Forward \n";
		phaser.CalcScaledForward();

		cout << "Starting Backward \n";
		phaser.CalcScaledBackward();

		cout << "Starting Posterior \n";
		int * ml_states = new int[phaser.num_markers];
		phaser.GetMLHaplotypes(ml_states);

		HaplotypePair mine_haps = phaser.PrintHaplotypesToFile(ml_states, file_name_out);
		HaplotypePair mine_haps2 = phaser.PrintReferenceHaplotypes(ml_states,file_name_out+"_ref_haps");
		mine_haps.print();
		printf("Equal = %d \n", mine_haps.isEqual(mine_haps2));

		delete [] ml_states;
		end= std::chrono::steady_clock::now();
		cout << "Time " << " coverage " << cov <<" = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;


//		true_haps.print();
	}

}
