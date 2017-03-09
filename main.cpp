#include <chrono>

#include "HaplotypePhaser.h"
#include "VcfUtils.h"
#include "Parameters.h"
#include "Pedigree.h"
//#include "VcfLoader.h"


using namespace std;


int main(int argc, char ** argv){
	String sample_file;
	String ref_file;

	ParameterList parameter_list;
	LongParamContainer long_parameters;
	long_parameters.addGroup("Input files");
	long_parameters.addString("sampleVCF", &sample_file);
	long_parameters.addString("referenceVCF", &ref_file);

	parameter_list.Add(new LongParameters("Options",long_parameters.getLongParameterList()));

	parameter_list.Read(argc, argv);
	parameter_list.Status();

	HaplotypePhaser phaser;
	//	NA18909
	phaser.LoadData(ref_file, sample_file, 626);
//	VcfUtils::CompareHaplotypes();

	printf("Num inds= %d num markers= %d num states = %d \n", phaser.num_inds, phaser.num_markers, phaser.num_states);
//	phaser.CalcForward();
//	phaser.CalcBackward();
//	phaser.CalcPosteriorOld();

//	Doing marker = 1 state = 34544

	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
	cout << "Starting Forward \n";
	phaser.CalcScaledForward();

	cout << "Starting Backward \n";
	phaser.CalcScaledBackward();

	cout << "Starting Posterior \n";
	phaser.CalcPosterior();
	chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

	cout << "Time difference = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " millisec" <<endl;
//
	int * ml_states = new int[phaser.num_markers];
	phaser.SampleHaplotypesNew(ml_states);
	phaser.PrintHaplotypes(ml_states);
	delete [] ml_states;
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
