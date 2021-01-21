// SoftwareModel.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#pragma warning( push )
#pragma warning( disable : 26451 )

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

struct Code_Params {
	int N = 0;
	int K = 0;
	double designSNRdb = 0;
	std::vector<unsigned int> generator_kernel_order;	// generator_kernel_order[0] -> rightmost stage in encoder xor graph
	std::vector<unsigned int> FrozenSet;
	double startSNRdB = 0;
	double endSNRdB = 0;
	double SNR_step_size = 0;
} CodeParameters;

// TODO: Reference additional headers your program requires here.
