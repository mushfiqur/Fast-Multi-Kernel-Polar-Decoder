#ifndef FROZEN_BITS_GENERATOR_HPP
#define FROZEN_BITS_GENERATOR_HPP

#include<vector>

class Frozen_Bits_GA_Generator {
private:
	// Taken from aff3ct (https://github.com/aff3ct/)
	const double alpha = -0.4527;
	const double beta = 0.0218;
	const double gamma = 0.8600;

	const double a = 1.0 / alpha;
	const double b = -beta / alpha;
	const double c = 1.0 / gamma;

	const double phi_pivot = 0.867861;
	const double phi_inv_pivot = 0.6845772418;

	int N;
	int K;
	double Rate;
	std::vector<unsigned int> kernel_order;
	double designSNRdb;
	double designSNR;

private:
	double phi(double x);
	double phi_inv(double x);

	std::vector<unsigned int> find_least_indices();
public:
	std::vector<double> z;
public:
	Frozen_Bits_GA_Generator(int N, int K, std::vector<unsigned int> kernel_order, double designSNRdb);
	std::vector<unsigned int> generate_frozen_set();
	std::vector<unsigned int> generate2();
	std::vector<unsigned int> generate3();
};

#endif