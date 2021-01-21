// SoftwareModel.cpp : Defines the entry point for the application.
//

#include "Launcher.hpp"
#include "SoftwareModel.hpp"

using namespace std;

int is_valid_CodeParameters() {
	if (CodeParameters.K > CodeParameters.N) {
		std::cout << "K cannot be greater than N" << "\n\n";
		return -1;
	}

	int size = 1;
	int num_bin_kernels = 0;
	for (int i = 0; i < CodeParameters.generator_kernel_order.size(); i++) {
		if (CodeParameters.generator_kernel_order[i] != 2 && CodeParameters.generator_kernel_order[i] != 3) {
			std::cout << "Kernels sizes other than 2 and 3 not supported";
			return -1;
		}
		if (CodeParameters.generator_kernel_order[i] == 2) {
			num_bin_kernels++;
		}
		size *= CodeParameters.generator_kernel_order[i];
	}
	if (num_bin_kernels == 0) {
		std::cout << "Needs at least 1 binary kernel (i.e.  [1 0; 1 1])";
		return -1;
	}
	if (size != CodeParameters.N) {
		std::cout << "Generator Kernel Order and N don't match" << "\n\n";
		return -1;
	}

	return 1;
}

int main()
{
	CodeParameters.N = 8192;
	CodeParameters.K = (int)(0.5 * CodeParameters.N);
	CodeParameters.designSNRdb = 0;
	CodeParameters.generator_kernel_order = std::vector<unsigned int>{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
	CodeParameters.startSNRdB = 0;
	CodeParameters.endSNRdB = 5;
	CodeParameters.SNR_step_size = 1;

	if (is_valid_CodeParameters() < 0)
		return -1;

	Launcher<int, float> l(CodeParameters.N, 
		CodeParameters.K, 
		CodeParameters.designSNRdb, 
		CodeParameters.generator_kernel_order, 
		CODE_CONSTRUCTION::MINIMUM_DISTANCE, 
		CodeParameters.startSNRdB,
		CodeParameters.endSNRdB,
		CodeParameters.SNR_step_size);

	l.run_fast_ssc(1);

	return 0;
}
