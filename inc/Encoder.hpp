#ifndef MULTI_KERNEL_ENCODER_HPP
#define MULTI_KERNEL_ENCODER_HPP

#include<vector>
#include <Eigen/Dense>

using namespace Eigen;

class Polar_Encoder {
private:	
	int N;
	int K;
	double Rate;
	std::vector<unsigned int> kernel_order;
	std::vector<int> frozen_indices;
	std::vector<unsigned int> codeword;
	
	// These can be changed form Dynamic to fixed if we set N_max (maybe N_max = 4096)
	std::vector< Matrix<int, 2, Dynamic> > Q_arr;
	std::vector< Matrix<int, 2, Dynamic> > P_arr;

	Matrix<int, Dynamic, Dynamic> right;
	Matrix<int, Dynamic, Dynamic> mult;

private:
	std::vector<int> uint_to_binary_vector(uint64_t data);
	
	void place_data_into_codeword(std::vector<int> data);

	int get_partial_product(int idx);

	void generate_Q_matrices();

	void generate_P_matrices();

public:
public:
	Polar_Encoder(int N, int K, std::vector<unsigned int> kernel_order, std::vector<unsigned int> FrozenSet);
	std::vector<unsigned int> encode(uint64_t data);
	std::vector<unsigned int> encode(std::vector<int> data);
};

#endif