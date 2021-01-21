#include "Encoder.hpp"
#include <cmath>
#include <iostream>
#include <sstream>

Polar_Encoder::Polar_Encoder(int N, int K, std::vector<unsigned int> kernel_order, std::vector<unsigned int> FrozenSet)
	: N(N), K(K), Rate(K/(double)N), kernel_order(kernel_order)
{
	frozen_indices = std::vector<int>(N);
	for (int i = 0; i < FrozenSet.size(); i++) {
		frozen_indices[FrozenSet[i]] = 1;
	}
	
	right = MatrixXd::Zero(N, N).cast<int>();
	mult = MatrixXd::Identity(N, N).cast<int>();

#define REAL_INDICES

	if (kernel_order.size() > 1) {
		generate_Q_matrices();
		generate_P_matrices();
	}

	std::cout << "\n";
	std::cout << "Generated Tanner Graph for code with parameters:" << "\n";
	std::cout << "    " << "N: " << N << "\n";
	std::cout << "    " << "K: " << K << "\n";
	std::cout << "    " << "Rate: " << Rate << "\n";

	std::cout << "    " << "G:";
	for (int i = 1; i <= kernel_order.size() - 1; i++) {
		std::cout << " T" << kernel_order[kernel_order.size() - i] << " *";
	}
	std::cout << " T" << kernel_order[0] << "\n";

	std::cout << "    " << "Frozen Set: { ";
	for (int i = 0; i < FrozenSet.size() - 1; i++) {
		std::cout << FrozenSet[i] << ", ";
	}
	std::cout << FrozenSet[FrozenSet.size() - 1] << " }" << "\n";

	/*std::cout << "\t" << "Interconnect Matrices: " << "\n";
	for (int i = 0; i < P_arr.size(); i++) {
		std::cout << "P_" << i + 1 << ": \n";
		std::cout << P_arr[i] << "\t";
		std::cout << "\n\n";
	}*/
}

void Polar_Encoder::generate_P_matrices() 
{
	std::cout << "Generating P matrices..." << "\n";
	P_arr = std::vector<Matrix<int, 2, Dynamic>>(kernel_order.size());

	// Top row of P is output index of stage on the left
	// Bottom row of P is intput index of stage on the right
	Matrix<int, 2, Dynamic> P = MatrixXd::Zero(2, N).cast<int>();;

	// Generate P_2 ... P_s
	for (int k = 1; k < kernel_order.size(); k++) {
		std::cout << "    " << " - Generating " << k << "/" << kernel_order.size() << "...\n";
		int idx = k + 1;
		int Ni_next = get_partial_product(idx + 1);

		Matrix<int, 2, Dynamic> Q = Q_arr[idx - 2];

		for (int i = 0; i < (N / Q.cols()); i++) {
			// All this because sizes need to match up when initializing P
			P << P.block(0, 0, 2, i * Q.cols()),				// Keep already updated section
				(Q.array() + i * (Ni_next)),					// New section
				P.block(0, 0, 2, N - ((i + 1) * Q.cols()));		// Keep section of old block

		}
		P_arr[k] = P;
	}

	// Generate P_1
	std::cout << "    " << " - Generating " << kernel_order.size() << "/" << kernel_order.size();
	for (int i = 1; i < P_arr.size(); i++) {
		std::cout << ".";
		Matrix<int, 2, Dynamic> P = P_arr[i];
		right = MatrixXd::Zero(N, N).cast<int>();
		for (int row_num = 0; row_num < N; row_num++) {
#ifdef REAL_INDICES
			right(row_num, P(1, row_num)) = 1;
#else
			right(row_num, P(1, row_num) - 1) = 1;
#endif
		}
		mult *= right;

	}

	mult = mult.cast<double>().inverse().cast<int>();
	
	for (int i = 0; i < mult.rows(); i++) {
		for (int j = 0; j < mult.cols(); j++) {
			if (mult(i, j) == 1) {
#ifdef REAL_INDICES
				P(1, i) = j;	
#else
				P(1, i) = j + 1;
#endif
			}
		}
	}

	std::cout << "DONE!" << "\n";
	P_arr[0] = P;
}

void Polar_Encoder::generate_Q_matrices()
{
	Q_arr = std::vector< Matrix<int, 2, Dynamic> >(kernel_order.size() - 1);
	for (int j = 0; j < kernel_order.size() - 1; j++) {
		std::cout << "Generating Q Matrix: " << j << "/" << kernel_order.size() - 2 << "\n";
		int idx = j + 2; 
		int col_idx_ctr = 0;
		int Ni = get_partial_product(idx);
		int ni = kernel_order[idx - 1];

		Matrix<int, 2, Dynamic> Q = MatrixXd::Zero(2, get_partial_product(idx + 1)).cast<int>();
		
		for (int i = 1; i <= Ni; i++) {
#ifdef REAL_INDICES
			Q(0, col_idx_ctr) = i - 1;
			Q(1, col_idx_ctr) = ((i - 1) * ni);
#else
			Q(0, col_idx_ctr) = i;
			Q(1, col_idx_ctr) = ((i - 1) * ni) + 1;
#endif
			col_idx_ctr++;
		}

		int iteration_ctr = 0;
		for (int i = Ni + 1; i < (ni - 1)*Ni + 1; i++) {
#ifdef REAL_INDICES
			Q(0, col_idx_ctr) = i - 1;
			Q(1, col_idx_ctr) = 1 + iteration_ctr * ni;
			iteration_ctr++;
			col_idx_ctr++;
#else
			Q(0, col_idx_ctr) = i;
			Q(1, col_idx_ctr) = 2 + iteration_ctr *ni;
			iteration_ctr++;
			col_idx_ctr++;
#endif
		}

		iteration_ctr = 0;
		for (int i = (ni - 1) * Ni + 1; i <= get_partial_product(idx+1); i++) {
#ifdef REAL_INDICES
			Q(0, col_idx_ctr) = i - 1;
			Q(1, col_idx_ctr) = (iteration_ctr + 1) * ni - 1;
			iteration_ctr++;
			col_idx_ctr++;
#else
			Q(0, col_idx_ctr) = i;
			Q(1, col_idx_ctr) = (iteration_ctr +1)*ni;
			iteration_ctr++;
			col_idx_ctr++;
#endif
		}
		Q_arr[j] = Q;

		/*
		for (int i = 2; i <= ni; i++) {
#ifdef REAL_INDICES
			Q(0, ctr) = ((i - 1) * Ni);
			Q(1, ctr) = (i - 1);
#else
			Q(0, ctr) = Ni+1 + (i-2);
			Q(1, ctr) = i;
			//Q(0, ctr) = ((i - 1) * Ni) + 1;
			//Q(1, ctr) = (i - 1) + 1;
#endif
			ctr++;
		}


		std::cout << "Q: \n";
		std::cout << Q << "\t";
		std::cout << "\n\n";


		for (int i = (ni - 1)*Ni + 1; i <= get_partial_product(idx + 1); i++) {
#ifdef REAL_INDICES
			Q(0, ctr) = ((ni - 1) * Ni + i);
			Q(1, ctr) = ((i + 1) * ni - 1);
#else
			Q(0, ctr) = i;
			Q(1, ctr) = 1;
			// Q(1, ctr) = ((i + 1) * ni - 1) + 1;
#endif
			ctr++;
		}

		*/
	}
}

int Polar_Encoder::get_partial_product(int idx)
{
	// Partial Product of kernel sizes in stage i
	int N_i = 1;

	for (int i = 0; i < idx - 1; i++) {
		N_i *= kernel_order[i];
	}

	return N_i;
	
}

std::vector<int> Polar_Encoder::uint_to_binary_vector(uint64_t data) {
	std::vector<int> data_vec( K, 0);
	
	for (int i = 0; i < data_vec.size(); i++) {
		data_vec[i] = (data >> i) & 0x01;
	}

	return data_vec;
}

std::vector<unsigned int> Polar_Encoder::encode(uint64_t data) {
	return encode(uint_to_binary_vector(data));
}

std::vector<unsigned int> Polar_Encoder::encode(std::vector<int> data) 
{
	// Place data vector into its appropriate position
	place_data_into_codeword(data);

	std::vector<unsigned int> cw_temp = codeword;
	int num_stages = kernel_order.size();
	Matrix<int, 2, Dynamic> P;

	// Encode
	for (int stage = 0; stage < num_stages; stage++) {
		if (kernel_order[kernel_order.size() - 1 - stage] == 2) {
			for (int j = 0; j < codeword.size() / 2; j++) {
				// x0 = u0 XOR u1
				codeword[2 * j] = codeword[2 * j] ^ codeword[2 * j + 1];
			}
		}
		else if (kernel_order[num_stages - stage - 1] == 3) {
			for (int j = 0; j < codeword.size() / 3; j++) {
				// x0 = u0 (XOR) u1
				// x1 = u0 (XOR) u2
				// x2 = u0 (XOR) u1 (XOR) u2

				int x0 = codeword[3 * j] ^ codeword[3 * j + 1];
				int x1 = codeword[3 * j] ^ codeword[3 * j + 2];
				int x2 = codeword[3 * j] ^ codeword[3 * j + 1] ^ codeword[3 * j + 2];

				codeword[3 * j] = x0;
				codeword[3 * j + 1] = x1;
				codeword[3 * j + 2] = x2;
			}
		}

		// Reshuffle indices
		if (kernel_order.size() > 1) {
			for (int i = 0; i < codeword.size(); i++) {
				P = P_arr[num_stages - stage - 1];
				cw_temp[P(0, i)] = codeword[P(1, i)];
			}
			codeword = cw_temp;
		}
	}

	return codeword;
}

void Polar_Encoder::place_data_into_codeword(std::vector<int> data) 
{
	int data_idx = 0;
	codeword = std::vector<unsigned int>(N);

	for (int i = 0; i < codeword.size(); i++) {
		if (frozen_indices[i] == 0)
			codeword[i] = data[data_idx++];
	}
}
