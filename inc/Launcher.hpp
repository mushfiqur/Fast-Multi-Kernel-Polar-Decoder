#ifndef LAUNCHER_HPP
#define LAUNCHER_HPP

#include "Min_Dist_Construction.hpp"
#include "Frozen_Bits_Generator.hpp"
#include "Encoder.hpp"
#include "SC_Decoder.hpp"
#include "Fast_SSC_Decoder.hpp"

#include <chrono>
#include <random>
#include <cmath>
#include <iomanip>

enum CODE_CONSTRUCTION {
	MINIMUM_DISTANCE,
	GAUSSIAN_APPROX
};

template <typename B, typename R>
class Launcher {
private:
	Min_Dist_Frozen_Set_Const* min_dist_f_set_generator;
	Frozen_Bits_GA_Generator* ga_f_set_generator;
	Polar_Encoder* encoder;
	Polar_SC_Decoder<B, R>* sc_decoder;
	Polar_Fast_SSC_Decoder<B, R>* fast_ssc_decoder;

	int N;
	int K;
	double rate;
	double designSNRdB;
	std::vector<unsigned int> generator_kernel_order;
	std::vector<unsigned int> FrozenSet;

	double snr_range_start;
	double snr_range_end;
	double snr_range_step_size;
public:
	Launcher(int N, int K, double designSNRdB, std::vector<unsigned int> kernel_order, CODE_CONSTRUCTION constr_method, double start_SNR, double end_SNR, double step_size)
		: N(N), K(K), rate((N * 1.0)/K), designSNRdB(designSNRdB), generator_kernel_order(kernel_order), snr_range_start(start_SNR), snr_range_end(end_SNR), snr_range_step_size(step_size),
		min_dist_f_set_generator(nullptr), ga_f_set_generator(nullptr), encoder(nullptr), sc_decoder(nullptr), fast_ssc_decoder(nullptr)
	{
		min_dist_f_set_generator = new Min_Dist_Frozen_Set_Const(generator_kernel_order, K);
		ga_f_set_generator = new Frozen_Bits_GA_Generator(N, K, generator_kernel_order, designSNRdB);

		if (constr_method == MINIMUM_DISTANCE)
			this->FrozenSet = min_dist_f_set_generator->generate_frozen_set();
		else if (constr_method == GAUSSIAN_APPROX)
			this->FrozenSet = ga_f_set_generator->generate_frozen_set();

		encoder = new Polar_Encoder(N, K, generator_kernel_order, FrozenSet);

		sc_decoder = new Polar_SC_Decoder<B, R>(N, K, generator_kernel_order, FrozenSet);
		fast_ssc_decoder = new Polar_Fast_SSC_Decoder<B, R>(N, K, generator_kernel_order, FrozenSet);
	}
	~Launcher() {
		delete min_dist_f_set_generator;
		delete ga_f_set_generator;
		delete encoder;
		delete sc_decoder;
		delete fast_ssc_decoder;
	}

	void run_sc(int tgt_frame_errors)
	{
		std::vector<int> modulated(N);
		std::vector<R> received_vec(N);
		std::vector<unsigned int> codeword;

		std::vector<int> data_to_encode(K);
		std::vector<B> decoded_data_vec(K);

		// FRA: Number of frames simulated
		// BE:  Number of bit errors
		// FE:  Number of frame errors
		// BER: Bit Error Rate (calculated by: BE/(FRA x K) K is message size)
		// FER: Frame Error Rate (calculated by: FE/FRA)

		double FRA = 0;
		double BE = 0;
		double FE = 0;
		double BER = 0;
		double FER = 0;

		double channel_variance = 0;

		bool do_FE_increment = false;

		std::default_random_engine generator;
		std::normal_distribution<double> dist;
		srand(time(NULL));
		
		// int tgt_frame_errors = 100;

		HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
		CONSOLE_SCREEN_BUFFER_INFO bufferInfo;

		print_header();

		for (double snrdB = snr_range_start; snrdB <= snr_range_end; snrdB += snr_range_step_size) {
			GetConsoleScreenBufferInfo(h, &bufferInfo);

			channel_variance = sqrt(1.0 / (2.0 * rate * (pow(10, 0.1 * snrdB))));
			dist = std::normal_distribution<double>(0, channel_variance);

			while (FE < tgt_frame_errors) {
				FRA++;
				// Generate Data
				for (int i = 0; i < CodeParameters.K; i++)
					data_to_encode[i] = rand() % 2;

				// Encode
				codeword = encoder->encode(data_to_encode);

				// Modulate
				for (int i = 0; i < codeword.size(); i++) {
					if (codeword[i] == 0)
						modulated[i] = 1;
					else if (codeword[i] == 1)
						modulated[i] = -1;
				}

				// Add noise
				for (int i = 0; i < modulated.size(); i++)
					received_vec[i] = modulated[i] + dist(generator);

				// Decode
				sc_decoder->_decode_siho(received_vec.data(), decoded_data_vec.data());

				// Calculate Metrics
				for (int i = 0; i < CodeParameters.K; i++) {
					if (decoded_data_vec[i] != data_to_encode[i]) {
						BE++;
						do_FE_increment = true;
					}
				}

				if (do_FE_increment)
					FE++;

				BER = (BE * 1.0) / (FRA * K);
				FER = (FE * 1.0) / FRA;

				print_result(snrdB, FRA, BE, FE, BER, FER, h, bufferInfo);
				do_FE_increment = false;
			}

			BE = 0;
			BER = 0;
			FE = 0;
			FER = 0;
			FRA = 0;
		}
	}

	void run_fast_ssc(int tgt_frame_errors)
	{
		std::vector<int> modulated(N);
		std::vector<R> received_vec(N);
		std::vector<unsigned int> codeword = std::vector<unsigned int>(N);

		std::vector<int> data_to_encode(K);
		std::vector<B> decoded_data_vec(K);

		// FRA: Number of frames simulated
		// BE:  Number of bit errors
		// FE:  Number of frame errors
		// BER: Bit Error Rate (calculated by: BE/(FRA x K) K is message size)
		// FER: Frame Error Rate (calculated by: FE/FRA)

		double FRA = 0;
		double BE = 0;
		double FE = 0;
		double BER = 0;
		double FER = 0;

		double channel_variance = 0;

		bool do_FE_increment = false;

		std::default_random_engine generator;
		std::normal_distribution<double> dist;
		srand(time(NULL));

		//int tgt_frame_errors = 100;

		HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
		CONSOLE_SCREEN_BUFFER_INFO bufferInfo;

		print_header();

		for (double snrdB = snr_range_start; snrdB <= snr_range_end; snrdB += snr_range_step_size) {
			GetConsoleScreenBufferInfo(h, &bufferInfo);

			// channel_variance = sqrt(1.0 / (2.0 * rate * (pow(10, 0.1 * snrdB))));

			channel_variance = 1.0 / pow(10, 0.1 * snrdB);
			dist = std::normal_distribution<double>(0, channel_variance);

			while (FE < tgt_frame_errors) {
				FRA++;
				// Generate Data
				for (int i = 0; i < CodeParameters.K; i++)
					data_to_encode[i] = rand() % 2;

				// Encode
				codeword = encoder->encode(data_to_encode);

				// Modulate
				for (int i = 0; i < codeword.size(); i++) {
					if (codeword[i] == 0)
						modulated[i] = 1;
					else if (codeword[i] == 1)
						modulated[i] = -1;
				}

				// Add noise
				for (int i = 0; i < modulated.size(); i++)
					received_vec[i] = modulated[i] + dist(generator);

				// Decode
				fast_ssc_decoder->_decode_siho(received_vec.data(), decoded_data_vec.data(), true);

				// Calculate Metrics
				for (int i = 0; i < CodeParameters.K; i++) {
					if (decoded_data_vec[i] != data_to_encode[i]) {
						BE++;
						do_FE_increment = true;
					}
				}

				if (do_FE_increment)
					FE++;

				BER = (BE * 1.0) / (FRA * K);
				FER = (FE * 1.0) / FRA;

				print_result(snrdB, FRA, BE, FE, BER, FER, h, bufferInfo);
				do_FE_increment = false;
			}

			BE = 0;
			BER = 0;
			FE = 0;
			FER = 0;
			FRA = 0;
		}
	}

	void print_header() 
	{
		std::cout << "\n";
		std::cout << "||============================================================================================||" << "\n";

		std::cout << "||" << std::setw(10) << "Eb/N0" 
			<< std::setw(7) << "||" 
			<< std::setw(8) << "FRA"
			<< std::setw(7) << "||"
			<< std::setw(7) << "BE"
			<< std::setw(7) << "||"
			<< std::setw(7) << "FE"
			<< std::setw(7) << "||"
			<< std::setw(9) << "BER"
			<< std::setw(8) << "||"
			<< std::setw(9) << "FER"
			<< std::setw(8) << "||"
			<< "\n";

		std::cout << "||============================================================================================||" << "\n";
	}

	void print_result(double snrdB, double FRA, double BE, double FE, double BER, double FER, HANDLE h, CONSOLE_SCREEN_BUFFER_INFO bufferInfo)
	{
		SetConsoleCursorPosition(h, bufferInfo.dwCursorPosition);
		
		std::cout << "||" 
			<< std::setw(10) << std::fixed << std::setprecision(3) << snrdB
			<< std::setw(7) << "||"
			<< std::setw(9) << std::fixed << std::setprecision(0) << FRA
			<< std::setw(6) << "||"
			<< std::setw(7) << std::fixed << std::setprecision(0) << BE
			<< std::setw(7) << "||"
			<< std::setw(7) << std::fixed << std::setprecision(0) << FE
			<< std::setw(7) << "||"
			<< std::setw(14) << std::scientific << BER
			<< std::setw(3) << "||"
			<< std::setw(14) << std::scientific << FER
			<< std::setw(3) << "||"
			<< "\n";

	}

	void exhaustive_search_find_min_dist_codewords() 
	{
		std::vector<int> modulated(CodeParameters.N);
		std::vector<unsigned int> codeword;

		std::vector<int> data_to_encode(CodeParameters.K);
		std::vector<int> decoded_data_vec(CodeParameters.K);

		bool decoding_successful = true;
		int num = 0;
		int j = 0;
		int weight = 0;
		int min_dist = K;

		for (int i = 0; i < pow(2, K); i++) {
			num = i;
			while (j < K) {
				data_to_encode[K - j - 1] = num % 2;
				num /= 2;
				j++;
			}

			codeword = encoder->encode(data_to_encode);

			for (int m = 0; m < codeword.size(); m++) {
				if (codeword[m] == 1)
					weight++;
			}
			if (weight < min_dist && weight != 0)
				min_dist = weight;

			weight = 0;
			j = 0;
		}

		std::cout << "Minimum distance: " << min_dist << "\n";

	}

	void test_random(int iterations)
	{
		std::vector<R> modulated(CodeParameters.N);
		std::vector<unsigned int> codeword;

		std::vector<int> data_to_encode(CodeParameters.K);
		std::vector<B> decoded_data_vec(CodeParameters.K);

		bool decoding_successful = true;
		srand(time(NULL));

		for (int iter = 0; iter < iterations; iter++) {
			decoding_successful = true;

			for (int i = 0; i < CodeParameters.K; i++) {
				data_to_encode[i] = rand() % 2;
			}

			codeword = encoder->encode(data_to_encode);

			// Modulate
			for (int i = 0; i < codeword.size(); i++) {
				if (codeword[i] == 0) {
					modulated[i] = 1;
				}
				else if (codeword[i] == 1) {
					modulated[i] = -1;
				}
			}

			// Decode
			std::chrono::time_point<std::chrono::steady_clock> t1, t2;
			long long duration;

			t1 = std::chrono::high_resolution_clock::now();

			sc_decoder->_decode_siho(modulated.data(), decoded_data_vec.data());
			
			t2 = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
			SC_runtimes.push_back(duration);

			for (int i = 0; i < CodeParameters.K; i++) {
				if (decoded_data_vec[i] != data_to_encode[i])
					decoding_successful = false;
			}

			if (!decoding_successful) {
				std::cout << "SC DECODING UNSUCCESSFUL!" << "\n";
				system("pause");
			}

			t1 = std::chrono::high_resolution_clock::now();

			fast_ssc_decoder->_decode_siho(modulated.data(), decoded_data_vec.data(), false);

			t2 = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
			FAST_SSC_runtimes.push_back(duration);

			for (int i = 0; i < CodeParameters.K; i++) {
				if (decoded_data_vec[i] != data_to_encode[i])
					decoding_successful = false;
			}

			if (!decoding_successful) {
				std::cout << "FAST SSC UNSUCCESSFUL!" << "\n";

				system("pause");
			}
			else {
				//std::cout << "Iteration " << iter+1 << "/5000: " << "SUCCESS\n";
			}
		}
		int SC_accumulator = 0;
		int FAST_SSC_accumulator = 0;
		for (int i = 0; i < FAST_SSC_runtimes.size(); i++) {
			SC_accumulator += SC_runtimes[i];
			FAST_SSC_accumulator += FAST_SSC_runtimes[i];
		}

		std::cout << "SC Decoding Avg. Runtime: " << (float)SC_accumulator / SC_runtimes.size() << "\n";
		std::cout << "FAST SSC Decoding Avg. Runtime: " << (float)FAST_SSC_accumulator / FAST_SSC_runtimes.size() << "\n";

	}
};

#endif