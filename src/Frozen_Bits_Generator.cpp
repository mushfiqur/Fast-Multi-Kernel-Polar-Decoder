#include "Frozen_Bits_Generator.hpp"

#include <cmath>
#include <queue>
#include <algorithm>
#include <iostream>

Frozen_Bits_GA_Generator::Frozen_Bits_GA_Generator(int N, int K, std::vector<unsigned int> kernel_order, double designSNRdb)
	: N(N), K(K), Rate(K/(double)N), kernel_order(kernel_order), designSNRdb(designSNRdb), designSNR(pow(10, designSNRdb/10))
{
	std::reverse(kernel_order.begin(), kernel_order.end());
}

std::vector<unsigned int> Frozen_Bits_GA_Generator::generate3()
{
	z = std::vector<double>(N, 4 * Rate * designSNR);

	int num_stages = kernel_order.size();
	int nodes_in_stage = 1;
	double T;

	for (int j = 0; j < num_stages; j++) {
		nodes_in_stage *= kernel_order[j];
		if (kernel_order[j] == 2) {
			for (int i = 1; i < nodes_in_stage; i++) {
				T = z[i];
				z[2*i - 1] = phi_inv(1 - pow(1 - phi(T), 2));
				z[2 * i] = 2.0 * T;
			}
		}
		else if (kernel_order[j] == 3) {
			for (int i = 0; i < nodes_in_stage; i++) {

			}
		}

	}

	return find_least_indices();
}

std::vector<unsigned int> Frozen_Bits_GA_Generator::generate2()
{
	z = std::vector<double>(N, 0);
	z[0] = 4 * Rate * designSNR;
	
	int num_stages = kernel_order.size();
	int nodes_in_stage = 1;
	double T;

	for (int j = 0; j < num_stages; j++) {
		 nodes_in_stage *= kernel_order[j];

		 if (kernel_order[j] == 2) {
			 for (int t = 0; t < nodes_in_stage / 2; t++) {
				 T = z[t];

				 z[t] = phi_inv(1 - pow(1 - phi(T), 2));
				 z[nodes_in_stage / 2 + t] = 2.0 * T;

				 //std::cout << "Using z[" << t << "] to update z[" << t << "] and z[" << nodes_in_stage / 2 + t << "]" << "\n";
			 }
		 }
		 else if (kernel_order[j] == 3) {
			 for (int t = 0; t < nodes_in_stage / 3; t++) {
				 T = z[t];

				 z[t] = phi_inv(1 - (1 - phi(phi_inv(1 - pow((1 - phi(T)), 2)))) * (1 - phi(T)));
				 z[nodes_in_stage / 3 + t] = phi_inv(1 - pow((1 - phi(T)), 2)) + T;
				 z[nodes_in_stage / 3 + 2*t] = 2.0 * T;


				 //std::cout << "Using z[" << t << "] to update z[" << t << "] and z[" << nodes_in_stage / 3 + t << "]" << "and z[" << 2 * nodes_in_stage / 3 + t << "]\n";
			 }
		 }
	}

	return find_least_indices();
}

std::vector<unsigned int> Frozen_Bits_GA_Generator::generate_frozen_set()
{
	z = std::vector<double> (N, 4 * Rate * designSNR);
	std::vector<double> temp(N, 4 * Rate * designSNR);	// Can be optimized to not need temp

	for (int j = 0; j < kernel_order.size(); j++) {
		for (int i = 0; i < N / kernel_order[j]; i++) {
			if (kernel_order[j] == 2) {
				z[2 * i] = phi_inv(1 - pow((1 - phi(temp[i])), 2));
				z[2 * i + 1] = 2.0 * temp[i];
			}
			else if (kernel_order[j] == 3) {
				z[3 * i] = phi_inv(1 - (1 - phi(phi_inv(1 - pow((1 - phi(temp[i])), 2)))) * (1 - phi(temp[i])));
				z[3 * i + 1] = phi_inv(1 - pow((1 - phi(temp[i])), 2)) + temp[i];
				z[3 * i + 2] = 2.0 * temp[i];
			}
		}
		temp = z;
	}

	return find_least_indices();
}

std::vector<unsigned int> Frozen_Bits_GA_Generator::find_least_indices() {
	std::vector<unsigned int> F_set;
	std::priority_queue<std::pair<double, int>> q;

	for (int i = 0; i < N; i++) {
		F_set.push_back(i);
		q.push(std::pair<double, int>(z[i], i));
	}

	for (int i = 0; i < K; i++) {
		F_set.erase(std::find(F_set.begin(), F_set.end(), q.top().second));
		q.pop();
	}

	return F_set;	// return F_set
}

double Frozen_Bits_GA_Generator::phi(double x) 
{
	if (x < phi_pivot)
		return std::exp(0.0564 * x * x - 0.48560 * x);
	else
		return std::exp(alpha * std::pow(x, gamma) + beta);
}

double Frozen_Bits_GA_Generator::phi_inv(double x) 
{
	if (x > phi_inv_pivot)
		return 4.304964539 * (1 - sqrt(1 + 0.9567131408 * std::log(x)));
	else
		return std::pow(a * std::log(x) + b, c);
}