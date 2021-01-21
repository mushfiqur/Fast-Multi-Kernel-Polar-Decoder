#include "Min_Dist_Construction.hpp"

Min_Dist_Frozen_Set_Const::Min_Dist_Frozen_Set_Const(std::vector<unsigned int> g_kernel_order, int K)
	: kernel_order(g_kernel_order), K(K), num_binary_kernel(0), dimension_Tp(1)
{
	T2 << 1, 0, 1, 1;
	T3 << 1, 1, 1, 1, 0, 1, 0, 1, 1;

	int i;

	for (i = 0; i < kernel_order.size(); i++) {
		if (g_kernel_order[i] == 2) {
			num_binary_kernel++;
		}
		else {
			break;
		}
	}

	if (i != kernel_order.size()) {
		Tp_order = std::vector<unsigned int>(kernel_order.size() - i);

		for (int j = i; j < kernel_order.size(); j++) {
			dimension_Tp *= kernel_order[j];
			Tp_order[j - i] = kernel_order[j];
		}

		if (Tp_order[0] == 2)
			Tp = T2;
		else if (Tp_order[0] == 3)
			Tp = T3;

		if (Tp_order.size() > 1) {
			for (int k = 1; k < Tp_order.size(); k++) {
				if (Tp_order[k] == 2)
					Tp = kron(Tp, T2);
				else if (Tp_order[k] == 3)
					Tp = kron(Tp, T3); 
			}
		}
	}

	Spectrum_Tp = calc_min_dist_spectrum_Tp();
	Spectrum_T2_n = calc_min_dist_spectrum_T2_n();
	Spectrum_Gn = calc_min_dist_spectrum_Gn();
	rN = calc_rN();
}

std::vector<unsigned int> Min_Dist_Frozen_Set_Const::generate_frozen_set()
{
	std::vector<unsigned int> information_set;
	std::vector<unsigned int> frozen_set = std::vector<unsigned int>(rN.cols());

	int l, c, q;

	for (int i = 0; i < rN.cols(); i++) {
		frozen_set[i] = i;
	}

	for (int j = 0; j < K; j++) {
		rN.maxCoeff(&l);
		c = l % dimension_Tp;
		q = l - c;

		rN[l] = 0;

		if (c > 0) {
			if (information_set.size() > 0) {
				// Remove all elements from  I[c - 1] + q
				for (int i = 0; i < I[c - 1].size(); i++) {
					information_set.erase(std::remove(information_set.begin(), information_set.end(), I[c - 1][i] + q), information_set.end());
				}
			}
			// Add all elements from  I[c] + q
			for (int i = 0; i < I[c].size(); i++) {
				information_set.push_back(I[c][i] + q);
			}
		}
		else {
			// Add all elements from  I[c] + q
			for (int i = 0; i < I[c].size(); i++) {
				information_set.push_back(I[c][i] + q);
			}
		}
	}

	std::sort(information_set.begin(), information_set.end());

	for (int i = 0; i < information_set.size(); i++) {
		frozen_set.erase(std::remove(frozen_set.begin(), frozen_set.end(), information_set[i]), frozen_set.end());
	}

	return frozen_set;
}

Matrix<int, Dynamic, Dynamic> Min_Dist_Frozen_Set_Const::kron(Matrix<int, Dynamic, Dynamic> left, Matrix<int, Dynamic, Dynamic> right)
{
	Matrix<int, Dynamic, Dynamic> prod = Matrix<int, Dynamic, Dynamic>(left.rows() * right.rows(), left.cols() * right.cols());
	Matrix<int, Dynamic, Dynamic> blk = Matrix<int, Dynamic, Dynamic>(right.rows(), right.cols());

	for (int row_idx = 0; row_idx < left.rows(); row_idx++) {
		for (int col_idx = 0; col_idx < left.cols(); col_idx++) {
			blk = left(row_idx, col_idx) * right;

			for (int i = 0; i < blk.rows(); i++) {
				for (int j = 0; j < blk.cols(); j++) {
					prod((row_idx * right.rows()) + i, (col_idx * right.cols()) + j) = blk(i, j);
				}
			}
		}		
	}

	// std::cout << prod << "\t";
	return prod;
}

int Min_Dist_Frozen_Set_Const::factorial(int n)
{
	int ret = 1;
	for (int i = n; i > 0; i--) {
		ret *= i;
	}

	return ret;
}

void Min_Dist_Frozen_Set_Const::generate_combinations_recursive(std::vector<std::vector<int>>& arr, std::vector<int>& indexes, std::vector<int>& combination, int initial_idx, int remaining_indexes) {
	if (remaining_indexes == 0) {
		arr.push_back(combination);
		return;
	}
	for (int i = initial_idx; i <= indexes.size() - remaining_indexes; ++i) {
		combination.push_back(indexes[i]);
		generate_combinations_recursive(arr, indexes, combination, i + 1, remaining_indexes - 1);
		combination.pop_back();
	}
}

RowVectorXi Min_Dist_Frozen_Set_Const::calc_min_dist_spectrum_Tp()
{
	//std::vector<int> spectrum = std::vector<int>(dimension_Tp);
	RowVectorXi spect = RowVectorXi(dimension_Tp);

	std::vector<std::vector<int>> rows_used_list = std::vector<std::vector<int>>(dimension_Tp);

	std::vector<int> indexes;
	std::vector<int> combination;
	std::vector<std::vector<int>> index_combinations;

	std::vector<RowVectorXi> dec_to_bin_LUT = std::vector<RowVectorXi>(pow(2, dimension_Tp));

	int max_min_dist;
	std::vector<int> max_min_dist_order = std::vector<int>();
	std::vector<int> min_dist_order = std::vector<int>();;

	for (int i = 0; i < dec_to_bin_LUT.size(); i++) {
		dec_to_bin_LUT[i] = convert_to_binary_vector(i, dimension_Tp);
	}

	for (int i = 0; i < dimension_Tp; i++) {
		indexes.push_back(i);
	}

	// For each possible code dimension
	for (int i = 0; i < dimension_Tp; i++) {

		// Choose (i + 1) rows from (dimension_Tp) rows;
		generate_combinations_recursive(index_combinations, indexes, combination, 0, i + 1);

		Matrix<int, Dynamic, Dynamic> Tp_sub_matrix = Matrix<int, Dynamic, Dynamic>(i + 1, Tp.cols()).setZero();

		max_min_dist = 0;

		// For each possible Generator
		for (int j = 0; j < index_combinations.size(); j++) {
			// Copy
			for (int k = 0; k < index_combinations[j].size(); k++) {
				for (int m = 0; m < Tp.cols(); m++) {
					Tp_sub_matrix(k, m) = Tp(index_combinations[j][k], m);
				}
			}

			int min_dist = Tp.cols() + 1;

			for (int k = 1; k < pow(2, i + 1); k++) {
				RowVectorXi data_word = dec_to_bin_LUT[k].tail(i + 1);
				RowVectorXi code_word = data_word * Tp_sub_matrix;

				int weight = 0;

				for (int m = 0; m < code_word.cols(); m++) {
					if (code_word[m] > 1)
						code_word[m] %= 2;

					if (code_word[m])
						weight++;
				}
				if (weight < min_dist) {
					min_dist = weight;
					min_dist_order = index_combinations[j];
				}
			}

			if (min_dist > max_min_dist) {
				max_min_dist = min_dist;
				max_min_dist_order = min_dist_order;
			}
		}

		spect[i] = max_min_dist;
		rows_used_list[i] = max_min_dist_order;

		max_min_dist_order.clear();
		max_min_dist_order.clear();
		index_combinations.clear();
	}

	I = rows_used_list;

	return spect;
}

RowVectorXi Min_Dist_Frozen_Set_Const::calc_min_dist_spectrum_T2_n()
{
	RowVector2i base;
	base << 2, 1;

	RowVectorXi spect = base;

	for (int i = 0; i < num_binary_kernel - 1; i++) {
		spect = kron(spect, base);
	}
	
	std::vector<int> spect_vec = std::vector<int>(spect.cols());

	for (int i = 0; i < spect.cols(); i++) {
		spect_vec[i] = spect[i];
	}

	std::sort(spect_vec.begin(), spect_vec.end(), std::greater<int>());

	for (int i = 0; i < spect_vec.size(); i++) {
		spect[i] = spect_vec[i];
	}

	return spect;
}

RowVectorXi Min_Dist_Frozen_Set_Const::calc_min_dist_spectrum_Gn()
{
	RowVectorXi spect;

	spect = kron(Spectrum_T2_n, Spectrum_Tp);
	std::vector<int> spect_Gn = std::vector<int>(Spectrum_Gn.cols());
	for (int i = 0; i < Spectrum_Gn.cols(); i++) {
		spect_Gn[i] = spect[i];
	}

	std::sort(spect_Gn.begin(), spect_Gn.end(), std::greater<int>());

	for (int i = 0; i < spect_Gn.size(); i++) {
		spect[i] = spect_Gn[i];
	}

	return spect;
}

RowVectorXi Min_Dist_Frozen_Set_Const::calc_rN()
{
	RowVector2i base;
	base << 1, 2;
	// base << 2, 1;

	RowVectorXi spect = base;
	RowVectorXi rN;

	for (int i = 0; i < num_binary_kernel - 1; i++) {
		spect = kron(spect, base);
	}

	rN = kron(spect, Spectrum_Tp);

	return rN;
}

RowVectorXi Min_Dist_Frozen_Set_Const::convert_to_binary_vector(int num, int len)
{
	RowVectorXi bin_vec = RowVectorXi(len);

	int i = 0;
	if (num == 0) {
		bin_vec = RowVectorXi::Zero(bin_vec.cols()).cast<int>();;
	}
	else {
		while (i < len) {
			bin_vec[len - i - 1] = num % 2;
			num /= 2;
			i++;
		}
	}

	return bin_vec;
}
