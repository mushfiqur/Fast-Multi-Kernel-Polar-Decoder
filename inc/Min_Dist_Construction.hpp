#ifndef MINIMUM_DISTANCE_FROZEN_SET_CONSTRUCTOR
#define MINIMUM_DISTANCE_FROZEN_SET_CONSTRUCTOR

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using namespace Eigen;

class Min_Dist_Frozen_Set_Const {
private:
	std::vector<unsigned int> kernel_order;
	std::vector<unsigned int> Tp_order;

	Matrix<int, Dynamic, Dynamic> Tp;

	RowVectorXi Spectrum_Tp;
	std::vector<std::vector<int>> I;

	RowVectorXi Spectrum_T2_n;

	RowVectorXi Spectrum_Gn;

	RowVectorXi rN;

	int num_binary_kernel;
	int dimension_Tp;

	Matrix<int, 2, 2> T2;
	Matrix<int, 3, 3> T3;

	int K;

private:
	Matrix<int, Dynamic, Dynamic> kron(Matrix<int, Dynamic, Dynamic> left, Matrix<int, Dynamic, Dynamic> right);
	int factorial(int n);
	RowVectorXi convert_to_binary_vector(int num, int len);

	void generate_combinations_recursive(std::vector<std::vector<int>>& arr, std::vector<int>& combinations_arr, std::vector<int>& indexes, int initial_idx, int remaining_indexes);

	RowVectorXi  calc_min_dist_spectrum_Tp();
	RowVectorXi  calc_min_dist_spectrum_T2_n();
	RowVectorXi  calc_min_dist_spectrum_Gn();
	RowVectorXi  calc_rN();

public:
	Min_Dist_Frozen_Set_Const(std::vector<unsigned int> g_kernel_order, int K);
	std::vector<unsigned int> generate_frozen_set();
};

#endif
