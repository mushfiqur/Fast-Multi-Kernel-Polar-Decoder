#ifndef MULTI_KERNEL_FAST_SSC_DECODER_HPP
#define MULTI_KERNEL_FAST_SSC_DECODER_HPP

#include <array>
#include <vector>
#include <algorithm>
#include "Tree.hpp"
#include "Decoder_Functions.hpp"
#include "Contents_SC.hpp"

#include <windows.h>

template <typename B, typename R>
class Polar_Fast_SSC_Decoder {
private:
	int N;
	int K;
	double Rate;
	std::vector<unsigned int> generator_kernel_order;
	std::vector<unsigned int> frozen_set;

	std::array<std::vector<int>, 3> Pv;

	Tree<Contents_SC<B, R>> polar_tree;

	int total_nodes_deleted;

	std::vector<B> u_hat;

private:
	void recursive_allocate_nodes_contents(Tree_Node<Contents_SC<B, R>>* node_curr, const int vector_size)
	{
		if (node_curr != nullptr) {
			node_curr->set_contents(new Contents_SC<B, R>(vector_size));

			if (node_curr->center == nullptr) {
				this->recursive_allocate_nodes_contents(node_curr->left, vector_size / 2);
				this->recursive_allocate_nodes_contents(node_curr->right, vector_size / 2);
			}
			else {
				this->recursive_allocate_nodes_contents(node_curr->left, vector_size / 3);
				this->recursive_allocate_nodes_contents(node_curr->center, vector_size / 3);
				this->recursive_allocate_nodes_contents(node_curr->right, vector_size / 3);
			}
		}
	};

	void recursive_initialize_frozen_bits(const Tree_Node<Contents_SC<B, R>>* node_curr)
	{
		auto* contents = node_curr->contents;

		if (!node_curr->is_leaf()) // stop condition
		{
			if (node_curr->center == nullptr) {
				recursive_initialize_frozen_bits(node_curr->left);
				recursive_initialize_frozen_bits(node_curr->right);
			}
			else {
				recursive_initialize_frozen_bits(node_curr->left);
				recursive_initialize_frozen_bits(node_curr->center);
				recursive_initialize_frozen_bits(node_curr->right);
			}
		}
		else
			contents->is_frozen_bit = frozen_set[node_curr->lane_id];
	};

	void recursive_deallocate_nodes_contents(Tree_Node<Contents_SC<B, R>>* node_curr)
	{
		if (node_curr != nullptr)
		{
			if (node_curr->center == nullptr) {
				this->recursive_deallocate_nodes_contents(node_curr->left);
				this->recursive_deallocate_nodes_contents(node_curr->right);
			}
			else {
				this->recursive_deallocate_nodes_contents(node_curr->left);
				this->recursive_deallocate_nodes_contents(node_curr->center);
				this->recursive_deallocate_nodes_contents(node_curr->right);

			}

			auto* contents = node_curr->contents;
			delete contents;
			node_curr->set_contents(nullptr);
		}
	};

	void recursive_delete_children(Tree_Node<Contents_SC<B, R>>* node_curr) 
	{
		if (node_curr->left != nullptr) {
			recursive_deallocate_nodes_contents(node_curr->left);
			polar_tree.delete_nodes(node_curr->left);
			node_curr->left = nullptr;
		}
		if (node_curr->center != nullptr) {
			recursive_deallocate_nodes_contents(node_curr->center);
			polar_tree.delete_nodes(node_curr->center);
			node_curr->center = nullptr;
		}
		if (node_curr->right != nullptr) {
			recursive_deallocate_nodes_contents(node_curr->right);
			polar_tree.delete_nodes(node_curr->right);
			node_curr->right = nullptr;
		}
	}

	void prune_tree(Tree_Node<Contents_SC<B, R>>* node_curr)
	{
		return;

		int code_len = 1;
		int curr_depth = node_curr->depth;
		int curr_lane_id = node_curr->lane_id;
		int num_children = curr_depth >= generator_kernel_order.size() ? 1 : generator_kernel_order[curr_depth];
		
		for (int j = generator_kernel_order.size(); j > curr_depth; j--) {
			code_len *= generator_kernel_order[j - 1];
		}

		int start_idx = code_len * curr_lane_id;
		int end_idx = start_idx + code_len - 1;

		//std::cout << "Depth: " << curr_depth << " Frozen Set: [" << start_idx << ", " << start_idx + code_len - 1 << "]" << "\n";

		int num_info_bits = 0;
		int num_frozen_bits = 0;

		// Figure out type of node
		for (int i = start_idx; i <= end_idx; i++) {
			if (frozen_set[i] == 1)
				num_frozen_bits++;
			else
				num_info_bits++;
		}

		if (num_info_bits == code_len && !node_curr->is_leaf()) {
			std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a Rate-1 Node" << "\n";
			node_curr->type = NODE_TYPE::RATE_1;

			recursive_delete_children(node_curr);
			
			int stage_nodes = 1;
			int size = 0;
			for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
				for (int j = curr_depth; j <= i; j++) {
					stage_nodes *= generator_kernel_order[j];
				}
				size += stage_nodes;

				stage_nodes = 1;
			}
			total_nodes_deleted += size;
			//std::cout << "Deleted " << size << " Nodes (Rate-1)" << "\n";

			return;
		}
		else if (num_frozen_bits == code_len && !node_curr->is_leaf()) {
			std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a Rate-0 Node" << "\n";
			node_curr->type = NODE_TYPE::RATE_0;

			recursive_delete_children(node_curr);

			int stage_nodes = 1;
			int size = 0;
			for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
				for (int j = curr_depth; j <= i; j++) {
					stage_nodes *= generator_kernel_order[j];
				}
				size += stage_nodes;

				stage_nodes = 1;
			}
			total_nodes_deleted += size;
			//std::cout << "Deleted " << size << " Nodes (Rate-0)" << "\n";

			return;
		}
		else if ((num_info_bits == code_len - 1) && frozen_set[start_idx] == 1 && !node_curr->is_leaf()) {
			std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a SPC Node" << "\n";
			node_curr->type = NODE_TYPE::SPC;

			recursive_delete_children(node_curr);

			int stage_nodes = 1;
			int size = 0;
			for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
				for (int j = curr_depth; j <= i; j++) {
					stage_nodes *= generator_kernel_order[j];
				}
				size += stage_nodes;

				stage_nodes = 1;
			}
			total_nodes_deleted += size;
			//std::cout << "Deleted " << size << " Nodes (SPC)" << "\n";

			return;
		}
		else if (num_frozen_bits == code_len - 1 && (frozen_set[end_idx] == 0) && !node_curr->is_leaf()) {
			//std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a REP Node" << "\n";

			int num_binary_kernels = 0;
			int num_ternary_kernels = 0;

			for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
				if (generator_kernel_order[i] == 2)
					num_binary_kernels++;
				else if (generator_kernel_order[i] == 3)
					num_ternary_kernels++;
			}

			if (num_ternary_kernels == 0) {
				std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a REP2 Node" << "\n";
				node_curr->type = NODE_TYPE::REP2;
				
				recursive_delete_children(node_curr);

				int stage_nodes = 1;
				int size = 0;
				for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
					for (int j = curr_depth; j <= i; j++) {
						stage_nodes *= generator_kernel_order[j];
					}
					size += stage_nodes;

					stage_nodes = 1;
				}
				total_nodes_deleted += size;
				//std::cout << "Deleted " << size << " Nodes (REP2)" << "\n";

				return;
			}
			else if (num_binary_kernels == 0) {
				if (code_len <= 27) {
					std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a REP3A Node" << "\n";
					node_curr->type = NODE_TYPE::REP3A;
					
					recursive_delete_children(node_curr);

					int stage_nodes = 1;
					int size = 0;
					for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
						for (int j = curr_depth; j <= i; j++) {
							stage_nodes *= generator_kernel_order[j];
						}
						size += stage_nodes;

						stage_nodes = 1;
					}
					total_nodes_deleted += size;
					//std::cout << "Deleted " << size << " Nodes (REP3A)" << "\n";

					return;
				}
				else {
					std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a REP3A Node but code length is greater than 27" << "\n";
				}
			}
			else if (generator_kernel_order.back() == 2) {
				std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a REP3B Node" << "\n";
				node_curr->type = NODE_TYPE::REP3B;

				recursive_delete_children(node_curr);

				int stage_nodes = 1;
				int size = 0;
				for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
					for (int j = curr_depth; j <= i; j++) {
						stage_nodes *= generator_kernel_order[j];
					}
					size += stage_nodes;

					stage_nodes = 1;
				}
				total_nodes_deleted += size;
				//std::cout << "Deleted " << size << " Nodes (REP3B)" << "\n";

				return;
			}
			else if (generator_kernel_order.back() == 3) {
				std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a REP3C Node" << "\n";
				node_curr->type = NODE_TYPE::REP3C;

				recursive_delete_children(node_curr);

				int stage_nodes = 1;
				int size = 0;
				for (int i = curr_depth; i < generator_kernel_order.size(); i++) {
					for (int j = curr_depth; j <= i; j++) {
						stage_nodes *= generator_kernel_order[j];
					}
					size += stage_nodes;

					stage_nodes = 1;
				}
				total_nodes_deleted += size;
				//std::cout << "Deleted " << size << " Nodes (REP3C)" << "\n";

				return;
			}
		}
		else {
			//std::cout << "Node (" << curr_depth << ", " << curr_lane_id << ") is a Rate-N" << "\n";
		}
			
		if (!node_curr->is_leaf()) {
			if (generator_kernel_order[curr_depth] == 2) {
				prune_tree(node_curr->left);
				prune_tree(node_curr->right);
			}
			else if (generator_kernel_order[curr_depth] == 3) {
				prune_tree(node_curr->left);
				prune_tree(node_curr->center);
				prune_tree(node_curr->right);
			}
		}
	};

	void show_program(const Tree_Node<Contents_SC<B, R>>* node_curr, int code_len)
	{
		switch (node_curr->type)
		{
		case NODE_TYPE::RATE_0:
		{
			std::cout << "ENABLE RATE-0 Decoder Block" << "\n";
			// Update Partial Sums
			std::cout << "STORE " << node_curr->contents->s.size() << " elements into beta RAM" << "\n";
			std::fill(node_curr->contents->s.begin(), node_curr->contents->s.end(), 0);

			// Update Decoded Bits
			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":"
				<< (code_len * node_curr->lane_id) + code_len - 1 << "] "
				<< "in codeword register (Rate-0 Decoder)" << "\n";

			std::fill(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len - 1, 0);

			break;
		}

		case NODE_TYPE::RATE_1:
		{
			std::cout << "ENABLE rate-1 Decoder Block" << "\n";
			std::cout << "STORE " << node_curr->contents->s.size() << " elements into beta RAM" << "\n";
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(node_curr->contents->lambda[i]);
			}

			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":" << (code_len * node_curr->lane_id) + code_len - 1 << "] in codeword register (Rate-1 Decoder)" << "\n";
			std::vector<B> temp = node_curr->contents->s;
			int ctr = code_len * node_curr->lane_id;

			for (int depth = node_curr->depth; depth < generator_kernel_order.size(); depth++) {
				if (generator_kernel_order[depth] == 2) {
					for (int i = 0; i < code_len / 2; i++) {
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 2)];

						u_hat[ctr++] = temp[i + (code_len / 2)];
					}
				}
				else if (generator_kernel_order[depth] == 3) {
					for (int i = 0; i < code_len / 3; i++) {
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						u_hat[ctr++] = temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						u_hat[ctr++] = temp[i] ^ temp[i + 2 * (code_len / 3)];
					}
				}
				ctr = code_len * node_curr->lane_id;
				temp = std::vector<B>(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len);
			}

			break;
		}

		case NODE_TYPE::SPC:
		{
			std::cout << "ENABLE SPC Decoder Block" << "\n";
			int parity = 0;

			// Update Partial Sums
			std::cout << "STORE " << node_curr->contents->s.size() << " elements into beta RAM" << "\n";

			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(node_curr->contents->lambda[i]);
				parity ^= node_curr->contents->s[i];
			}

			if (parity) {
				int least_reliable_idx = std::min_element(node_curr->contents->lambda.begin(), node_curr->contents->lambda.end()) - node_curr->contents->lambda.begin();
				node_curr->contents->s[least_reliable_idx] = node_curr->contents->s[least_reliable_idx] ^ parity;
			}

			std::vector<B> temp = node_curr->contents->s;
			int ctr = code_len * node_curr->lane_id;

			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":" << (code_len * node_curr->lane_id) + code_len - 1 << "] in codeword register (SPC Decoder)" << "\n";

			for (int depth = node_curr->depth; depth < generator_kernel_order.size(); depth++) {
				if (generator_kernel_order[depth] == 2) {
					for (int i = 0; i < code_len / 2; i++) {
						//std::cout << "temp[" << ctr << "] is " << s0 << " (XOR) " << s1 << "\n";
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 2)];

						//std::cout << "temp[" << ctr << "] is " << s1 << "\n";
						u_hat[ctr++] = temp[i + (code_len / 2)];
					}
				}
				else if (generator_kernel_order[depth] == 3) {
					for (int i = 0; i < code_len / 3; i++) {
						//std::cout << "temp[" << ctr << "] is " << s0 << " (XOR) " << s1 << " (XOR) " << s2 << "\n";
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						//std::cout << "temp[" << ctr << "] is " << s1 << " (XOR) " << s2 << "\n";
						u_hat[ctr++] = temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						//std::cout << "temp[" << ctr << "] is " << s0 << " (XOR) " << s2 << "\n";
						u_hat[ctr++] = temp[i] ^ temp[i + 2 * (code_len / 3)];
					}
				}
				ctr = code_len * node_curr->lane_id;
				temp = std::vector<B>(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len);
			}

			//do_rate_1_decode(node_curr, code_len);

			break;
		}

		case NODE_TYPE::REP2:
		{
			std::cout << "ENABLE REP2 Decoder Block" << "\n";

			int sum_of_LLRs = 0;
			for (int i = 0; i < node_curr->contents->lambda.size(); i++) {
				sum_of_LLRs += node_curr->contents->lambda[i];
			}

			// Update Partial Sums
			std::cout << "STORE " << node_curr->contents->s.size() << " elements into beta RAM" << "\n";
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(sum_of_LLRs);
			}

			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":" << (code_len * node_curr->lane_id) + code_len - 2 << "] = 0 in codeword register (REP2 Decoder)" << "\n";
			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			std::cout << "STORE bit[" << (code_len * node_curr->lane_id) + code_len - 1 << "] in codeword register (REP2 Decoder)" << "\n";
			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}
		case NODE_TYPE::REP3A:
		{
			std::cout << "ENABLE REP3A Decoder Block" << "\n";

			int sum_of_LLRs = 0;

			// Update Partial Sums

			std::cout << "STORE " << node_curr->contents->s.size() << " elements into beta RAM" << "\n";

			if (code_len == 3) {
				for (int i = 0; i < code_len; i++) {
					sum_of_LLRs += node_curr->contents->lambda[i] * Pv[0][i];
				}

				for (int i = 0; i < node_curr->contents->s.size(); i++) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs) * Pv[0][i];
				}
			}
			else if (code_len == 9) {
				for (int i = 0; i < code_len; i++) {
					sum_of_LLRs += node_curr->contents->lambda[i] * Pv[1][i];
				}

				for (int i = 0; i < node_curr->contents->s.size(); i++) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs) * Pv[1][i];
				}
			}
			else { // code_len == 27
				for (int i = 0; i < code_len; i++) {
					sum_of_LLRs += node_curr->contents->lambda[i] * Pv[2][i];
				}

				for (int i = 0; i < node_curr->contents->s.size(); i++) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs) * Pv[2][i];
				}
			}

			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":" << (code_len * node_curr->lane_id) + code_len - 2 << "] = 0 in codeword register (REP3A Decoder)" << "\n";

			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			std::cout << "STORE bit[" << (code_len * node_curr->lane_id) + code_len - 1 << "] in codeword register (REP3 Decoder)" << "\n";
			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}

		case NODE_TYPE::REP3B:
		{
			std::cout << "ENABLE REP3B Decoder Block" << "\n";

			int sum_of_LLRs = 0;

			for (int i = code_len / 3; i < code_len; i++) {
				sum_of_LLRs += node_curr->contents->lambda[i];
			}

			// Update Partial Sums
			std::cout << "STORE " << node_curr->contents->s.size() - (code_len / 3) << " elements into beta RAM" << "\n";

			for (int i = code_len / 3; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(sum_of_LLRs);
			}

			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":" << (code_len * node_curr->lane_id) + code_len - 2 << "] = 0 in codeword register (REP3B Decoder)" << "\n";
			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			std::cout << "STORE bit[" << (code_len * node_curr->lane_id) + code_len - 1 << "] in codeword register (REP3B Decoder)" << "\n";
			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}

		case NODE_TYPE::REP3C:
		{
			std::cout << "ENABLE REP3C Decoder Block" << "\n";

			int sum_of_LLRs = 0;

			for (int i = code_len / 3; i < code_len; i++) {
				sum_of_LLRs += node_curr->contents->lambda[i];
			}

			// Update Partial Sums
			std::cout << "STORE " << node_curr->contents->s.size() << " elements into beta RAM" << "\n";
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				if (i % 3 != 0) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs);
				}
			}

			std::cout << "STORE bits[" << (code_len * node_curr->lane_id) << ":" << (code_len * node_curr->lane_id) + code_len - 2 << "] = 0 in codeword register (REP3C Decoder)" << "\n";
			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			std::cout << "STORE bit[" << (code_len * node_curr->lane_id) + code_len - 1 << "] in codeword register (REP3C Decoder)" << "\n";
			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}

		case NODE_TYPE::RATE_N:
		{
			if (!node_curr->is_leaf()) // stop condition
			{
				const auto size = (int)node_curr->contents->lambda.size();

				if (node_curr->center == nullptr) {
					const auto size_2 = size / 2;

					const auto* node_left = node_curr->left; // get left node
					const auto* node_right = node_curr->right; // get right node

					// apply f()
					std::cout << "LOAD " << size_2 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_2 << " elements from alpha RAM" << "\n";
					std::cout << "Calculate F()" << "\n";
					std::cout << "STORE " << size_2 << " elements into alpha RAM (left node)" << "\n";

					for (auto i = 0; i < size_2; i++)
						node_left->contents->lambda[i] = F(node_curr->contents->lambda[i],
							node_curr->contents->lambda[size_2 + i]);

					this->show_program(node_left, size_2); // recursive call

					std::cout << "LOAD " << size_2 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_2 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_2 << " elements from beta RAM" << "\n";
					std::cout << "Calculate G()" << "\n";
					std::cout << "STORE " << size_2 << " elements into alpha RAM (right node)" << "\n";
					// apply g()
					for (auto i = 0; i < size_2; i++)
						node_right->contents->lambda[i] = G(node_curr->contents->lambda[i],
							node_curr->contents->lambda[size_2 + i],
							node_left->contents->s[i]);

					this->show_program(node_right, size_2); // recursive call

					// Combine
					if (!node_curr->is_root()) {
						std::cout << "STORE " << 2 * size_2 << " elements into beta RAM" << "\n";
						for (auto i = 0; i < size_2; i++)
							node_curr->contents->s[i] = node_left->contents->s[i] ^ node_right->contents->s[i]; // bit xor

						for (auto i = 0; i < size_2; i++)
							node_curr->contents->s[size_2 + i] = node_right->contents->s[i]; // bit eq
					}
				}
				else {
					const auto size_3 = size / 3;

					const auto* node_left = node_curr->left; // get left node
					const auto* node_center = node_curr->center; // get center node
					const auto* node_right = node_curr->right; // get right node

					// Apply new F()
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "Calculate f_Ternary()" << "\n";
					std::cout << "STORE " << size_3 << " elements into alpha RAM (left node)" << "\n";

					for (auto i = 0; i < size_3; i++)
						node_left->contents->lambda[i] = F_Ternary(node_curr->contents->lambda[i],  // apply f()
							node_curr->contents->lambda[size_3 + i],
							node_curr->contents->lambda[(2 * size_3) + i]);

					this->show_program(node_left, size_3); // recursive call

					// Apply G1()
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from beta RAM" << "\n";
					std::cout << "Calculate G1()" << "\n";
					std::cout << "STORE " << size_3 << " elements into alpha RAM (center, node)" << "\n";

					for (auto i = 0; i < size_3; i++)
						node_center->contents->lambda[i] = G_1(node_curr->contents->lambda[i], // apply g()
							node_curr->contents->lambda[size_3 + i],
							node_curr->contents->lambda[(2 * size_3) + i],
							node_left->contents->s[i]);

					this->show_program(node_center, size_3);

					// Apply G2()
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from alpha RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from beta RAM" << "\n";
					std::cout << "LOAD " << size_3 << " elements from beta RAM" << "\n";
					std::cout << "Calculate G2()" << "\n";
					std::cout << "STORE " << size_3 << " elements into alpha RAM (right node)" << "\n";

					for (auto i = 0; i < size_3; i++)
						node_right->contents->lambda[i] = G_2(node_curr->contents->lambda[i], // apply g()
							node_curr->contents->lambda[size_3 + i],
							node_curr->contents->lambda[(2 * size_3) + i],
							node_left->contents->s[i],
							node_center->contents->s[i]);

					this->show_program(node_right, size_3);

					// Combine
					if (!node_curr->is_root()) {
						std::cout << "STORE " << 3 * size_3 << " elements into beta RAM" << "\n";
						for (auto i = 0; i < size_3; i++)
							node_curr->contents->s[i] = node_left->contents->s[i] ^ node_center->contents->s[i]; // bit xor

						for (auto i = 0; i < size_3; i++)
							node_curr->contents->s[size_3 + i] = node_left->contents->s[i] ^ node_right->contents->s[i]; // bit xor

						for (auto i = 0; i < size_3; i++)
							node_curr->contents->s[(2 * size_3) + i] = node_left->contents->s[i] ^ node_center->contents->s[i] ^ node_right->contents->s[i]; // bit xor
					}
				}
			}
			else // specific leaf treatment
			{
				node_curr->contents->s[0] = (!node_curr->contents->is_frozen_bit && // if this is a frozen bit then s == 0
					H<B, R>(node_curr->contents->lambda[0])); // apply h()
				u_hat[node_curr->lane_id] = node_curr->contents->s[0];
			}

			//do_SC_decode(node_curr, code_len);
			break;
		}
		}
	};

	void recursive_decode(const Tree_Node<Contents_SC<B, R>>* node_curr, int code_len)
	{
		switch (node_curr->type)
		{
		case NODE_TYPE::RATE_0:
		{
			// Update Partial Sums
			std::fill(node_curr->contents->s.begin(), node_curr->contents->s.end(), 0);

			// Update Decoded Bits
			std::fill(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len - 1, 0);

			break;
		}

		case NODE_TYPE::RATE_1:
		{
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(node_curr->contents->lambda[i]);
			}

			std::vector<B> temp = node_curr->contents->s;
			int ctr = code_len * node_curr->lane_id;

			for (int depth = node_curr->depth; depth < generator_kernel_order.size(); depth++) {
				if (generator_kernel_order[depth] == 2) {
					for (int i = 0; i < code_len / 2; i++) {
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 2)];

						u_hat[ctr++] = temp[i + (code_len / 2)];
					}
				}
				else if (generator_kernel_order[depth] == 3) {
					for (int i = 0; i < code_len / 3; i++) {
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						u_hat[ctr++] = temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						u_hat[ctr++] = temp[i] ^ temp[i + 2 * (code_len / 3)];
					}
				}
				ctr = code_len * node_curr->lane_id;
				temp = std::vector<B>(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len);
			}

			break;
		}

		case NODE_TYPE::SPC:
		{
			int parity = 0;

			// Update Partial Sums
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(node_curr->contents->lambda[i]);
				parity ^= node_curr->contents->s[i];
			}

			if (parity) {
				int least_reliable_idx = std::min_element(node_curr->contents->lambda.begin(), node_curr->contents->lambda.end()) - node_curr->contents->lambda.begin();
				node_curr->contents->s[least_reliable_idx] = node_curr->contents->s[least_reliable_idx] ^ parity;
			}

			std::vector<B> temp = node_curr->contents->s;
			int ctr = code_len * node_curr->lane_id;


			for (int depth = node_curr->depth; depth < generator_kernel_order.size(); depth++) {
				if (generator_kernel_order[depth] == 2) {
					for (int i = 0; i < code_len / 2; i++) {
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 2)];

						u_hat[ctr++] = temp[i + (code_len / 2)];
					}
				}
				else if (generator_kernel_order[depth] == 3) {
					for (int i = 0; i < code_len / 3; i++) {
						u_hat[ctr++] = temp[i] ^ temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						u_hat[ctr++] = temp[i + (code_len / 3)] ^ temp[i + 2 * (code_len / 3)];

						u_hat[ctr++] = temp[i] ^ temp[i + 2 * (code_len / 3)];
					}
				}
				ctr = code_len * node_curr->lane_id;
				temp = std::vector<B>(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len);
			}

			//do_rate_1_decode(node_curr, code_len);

			break;
		}

		case NODE_TYPE::REP2:
		{

			int sum_of_LLRs = 0;
			for (int i = 0; i < node_curr->contents->lambda.size(); i++) {
				sum_of_LLRs += node_curr->contents->lambda[i];
			}

			// Update Partial Sums
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(sum_of_LLRs);
			}

			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}
		case NODE_TYPE::REP3A:
		{

			int sum_of_LLRs = 0;

			// Update Partial Sums


			if (code_len == 3) {
				for (int i = 0; i < code_len; i++) {
					sum_of_LLRs += node_curr->contents->lambda[i] * Pv[0][i];
				}

				for (int i = 0; i < node_curr->contents->s.size(); i++) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs) * Pv[0][i];
				}
			}
			else if (code_len == 9) {
				for (int i = 0; i < code_len; i++) {
					sum_of_LLRs += node_curr->contents->lambda[i] * Pv[1][i];
				}

				for (int i = 0; i < node_curr->contents->s.size(); i++) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs) * Pv[1][i];
				}
			}
			else { // code_len == 27
				for (int i = 0; i < code_len; i++) {
					sum_of_LLRs += node_curr->contents->lambda[i] * Pv[2][i];
				}

				for (int i = 0; i < node_curr->contents->s.size(); i++) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs) * Pv[2][i];
				}
			}


			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}

		case NODE_TYPE::REP3B:
		{

			int sum_of_LLRs = 0;

			for (int i = code_len / 3; i < code_len; i++) {
				sum_of_LLRs += node_curr->contents->lambda[i];
			}

			// Update Partial Sums

			for (int i = code_len / 3; i < node_curr->contents->s.size(); i++) {
				node_curr->contents->s[i] = H<B, R>(sum_of_LLRs);
			}

			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}

		case NODE_TYPE::REP3C:
		{

			int sum_of_LLRs = 0;

			for (int i = code_len / 3; i < code_len; i++) {
				sum_of_LLRs += node_curr->contents->lambda[i];
			}

			// Update Partial Sums
			for (int i = 0; i < node_curr->contents->s.size(); i++) {
				if (i % 3 != 0) {
					node_curr->contents->s[i] = H<B, R>(sum_of_LLRs);
				}
			}

			for (int i = (code_len * node_curr->lane_id); i <= (code_len * node_curr->lane_id) + code_len - 2; i++) {
				u_hat[i] = 0;
			}

			u_hat[(code_len * node_curr->lane_id) + code_len - 1] = H<B, R>(sum_of_LLRs);

			break;
		}

		case NODE_TYPE::RATE_N:
		{
			if (!node_curr->is_leaf()) // stop condition
			{
				const auto size = (int)node_curr->contents->lambda.size();

				if (node_curr->center == nullptr) {
					const auto size_2 = size / 2;

					const auto* node_left = node_curr->left; // get left node
					const auto* node_right = node_curr->right; // get right node

					// apply f()

					for (auto i = 0; i < size_2; i++)
						node_left->contents->lambda[i] = F(node_curr->contents->lambda[i],
							node_curr->contents->lambda[size_2 + i]);

					this->recursive_decode(node_left, size_2); // recursive call

					// apply g()
					for (auto i = 0; i < size_2; i++)
						node_right->contents->lambda[i] = G(node_curr->contents->lambda[i],
							node_curr->contents->lambda[size_2 + i],
							node_left->contents->s[i]);

					this->recursive_decode(node_right, size_2); // recursive call

					// Combine
					if (!node_curr->is_root()) {
						for (auto i = 0; i < size_2; i++)
							node_curr->contents->s[i] = node_left->contents->s[i] ^ node_right->contents->s[i]; // bit xor

						for (auto i = 0; i < size_2; i++)
							node_curr->contents->s[size_2 + i] = node_right->contents->s[i]; // bit eq
					}
				}
				else {
					const auto size_3 = size / 3;

					const auto* node_left = node_curr->left; // get left node
					const auto* node_center = node_curr->center; // get center node
					const auto* node_right = node_curr->right; // get right node

					// Apply new F()

					for (auto i = 0; i < size_3; i++)
						node_left->contents->lambda[i] = F_Ternary(node_curr->contents->lambda[i],  // apply f()
							node_curr->contents->lambda[size_3 + i],
							node_curr->contents->lambda[(2 * size_3) + i]);

					this->recursive_decode(node_left, size_3); // recursive call

					// Apply G1()

					for (auto i = 0; i < size_3; i++)
						node_center->contents->lambda[i] = G_1(node_curr->contents->lambda[i], // apply g()
							node_curr->contents->lambda[size_3 + i],
							node_curr->contents->lambda[(2 * size_3) + i],
							node_left->contents->s[i]);

					this->recursive_decode(node_center, size_3);

					// Apply G2()

					for (auto i = 0; i < size_3; i++)
						node_right->contents->lambda[i] = G_2(node_curr->contents->lambda[i], // apply g()
							node_curr->contents->lambda[size_3 + i],
							node_curr->contents->lambda[(2 * size_3) + i],
							node_left->contents->s[i],
							node_center->contents->s[i]);

					this->recursive_decode(node_right, size_3);

					// Combine
					if (!node_curr->is_root()) {
						for (auto i = 0; i < size_3; i++)
							node_curr->contents->s[i] = node_left->contents->s[i] ^ node_center->contents->s[i]; // bit xor

						for (auto i = 0; i < size_3; i++)
							node_curr->contents->s[size_3 + i] = node_left->contents->s[i] ^ node_right->contents->s[i]; // bit xor

						for (auto i = 0; i < size_3; i++)
							node_curr->contents->s[(2 * size_3) + i] = node_left->contents->s[i] ^ node_center->contents->s[i] ^ node_right->contents->s[i]; // bit xor
					}
				}
			}
			else // specific leaf treatment
			{
				node_curr->contents->s[0] = (!node_curr->contents->is_frozen_bit && // if this is a frozen bit then s == 0
					H<B, R>(node_curr->contents->lambda[0])); // apply h()
				u_hat[node_curr->lane_id] = node_curr->contents->s[0];
			}

			//do_SC_decode(node_curr, code_len);
			break;
		}
		}
	};

	void do_rate_1_decode(const Tree_Node<Contents_SC<B, R>>* node_curr, int code_len)
	{
		std::vector<B> temp = node_curr->contents->s;
		int ctr = code_len * node_curr->lane_id;

		for (int depth = node_curr->depth; depth < generator_kernel_order.size(); depth++) {
			if (generator_kernel_order[depth] == 2) {
				for (int i = 0; i < code_len / 2; i++) {
					int s0 = temp[i];
					int s1 = temp[i + (code_len / 2)];

					//std::cout << "temp[" << ctr << "] is " << s0 << " (XOR) " << s1 << "\n";
					u_hat[ctr++] = s0 ^ s1;

					//std::cout << "temp[" << ctr << "] is " << s1 << "\n";
					u_hat[ctr++] = s1;
				}
			}
			else if (generator_kernel_order[depth] == 3) {
				for (int i = 0; i < code_len / 3; i++) {
					int s0 = temp[i];
					int s1 = temp[i + (code_len / 3)];
					int s2 = temp[i + 2 * (code_len / 3)];

					//std::cout << "temp[" << ctr << "] is " << s0 << " (XOR) " << s1 << " (XOR) " << s2 << "\n";
					u_hat[ctr++] = s0 ^ s1 ^ s2;

					//std::cout << "temp[" << ctr << "] is " << s1 << " (XOR) " << s2 << "\n";
					u_hat[ctr++] = s1 ^ s2;

					//std::cout << "temp[" << ctr << "] is " << s0 << " (XOR) " << s2 << "\n";
					u_hat[ctr++] = s0 ^ s2;
				}
			}
			ctr = code_len * node_curr->lane_id;
			temp = std::vector<B>(u_hat.begin() + (code_len * node_curr->lane_id), u_hat.begin() + (code_len * node_curr->lane_id) + code_len);

			/*for (int i = 0; i < temp.size(); i++) {
				std::cout << temp[i] << " ";
			}
			std::cout << "\n";*/
		}
	}

	void do_SC_decode(const Tree_Node<Contents_SC<B, R>>* node_curr, int code_len)
	{
		if (!node_curr->is_leaf()) // stop condition
		{
			const auto size = (int)node_curr->contents->lambda.size();

			if (node_curr->center == nullptr) {
				const auto size_2 = size / 2;

				const auto* node_left = node_curr->left; // get left node
				const auto* node_right = node_curr->right; // get right node

				// apply f()
				for (auto i = 0; i < size_2; i++)
					node_left->contents->lambda[i] = F(node_curr->contents->lambda[i],
						node_curr->contents->lambda[size_2 + i]);

				this->recursive_decode(node_left, size_2); // recursive call

				// apply g()
				for (auto i = 0; i < size_2; i++)
					node_right->contents->lambda[i] = G(node_curr->contents->lambda[i],
						node_curr->contents->lambda[size_2 + i],
						node_left->contents->s[i]);

				this->recursive_decode(node_right, size_2); // recursive call

				// Combine
				for (auto i = 0; i < size_2; i++)
					node_curr->contents->s[i] = node_left->contents->s[i] ^ node_right->contents->s[i]; // bit xor

				for (auto i = 0; i < size_2; i++)
					node_curr->contents->s[size_2 + i] = node_right->contents->s[i]; // bit eq
			}
			else {
				const auto size_3 = size / 3;

				const auto* node_left = node_curr->left; // get left node
				const auto* node_center = node_curr->center; // get center node
				const auto* node_right = node_curr->right; // get right node

				// Apply new F()
				for (auto i = 0; i < size_3; i++)
					node_left->contents->lambda[i] = F_Ternary(node_curr->contents->lambda[i],  // apply f()
						node_curr->contents->lambda[size_3 + i],
						node_curr->contents->lambda[(2 * size_3) + i]);

				this->recursive_decode(node_left, size_3); // recursive call

				// Apply G1()
				for (auto i = 0; i < size_3; i++)
					node_center->contents->lambda[i] = G_1(node_curr->contents->lambda[i], // apply g()
						node_curr->contents->lambda[size_3 + i],
						node_curr->contents->lambda[(2 * size_3) + i],
						node_left->contents->s[i]);

				this->recursive_decode(node_center, size_3);

				// Apply G2()
				for (auto i = 0; i < size_3; i++)
					node_right->contents->lambda[i] = G_2(node_curr->contents->lambda[i], // apply g()
						node_curr->contents->lambda[size_3 + i],
						node_curr->contents->lambda[(2 * size_3) + i],
						node_left->contents->s[i],
						node_center->contents->s[i]);

				this->recursive_decode(node_right, size_3);

				// Combine
				for (auto i = 0; i < size_3; i++)
					node_curr->contents->s[i] = node_left->contents->s[i] ^ node_center->contents->s[i]; // bit xor

				for (auto i = 0; i < size_3; i++)
					node_curr->contents->s[size_3 + i] = node_left->contents->s[i] ^ node_right->contents->s[i]; // bit xor

				for (auto i = 0; i < size_3; i++)
					node_curr->contents->s[(2 * size_3) + i] = node_left->contents->s[i] ^ node_center->contents->s[i] ^ node_right->contents->s[i]; // bit xor

			}


		}
		else // specific leaf treatment
		{
			node_curr->contents->s[0] = (!node_curr->contents->is_frozen_bit && // if this is a frozen bit then s == 0
				H<B, R>(node_curr->contents->lambda[0])); // apply h()
			u_hat[node_curr->lane_id] = node_curr->contents->s[0];
		}
	}
public:
	Polar_Fast_SSC_Decoder(int N, int K, std::vector<unsigned int> generator_kernel_order, std::vector<unsigned int> FrozenSet)
		: N(N), K(K), Rate(K / (double)N), total_nodes_deleted(0),  generator_kernel_order(generator_kernel_order),
		polar_tree(generator_kernel_order)
	{
		frozen_set = std::vector<unsigned int>(N);
		for (int i = 0; i < FrozenSet.size(); i++) {
			frozen_set[FrozenSet[i]] = 1;
		}

		// Supported Pv patterns
		// (See: A. Cavatassi, T. Tonnellier and W. J. Gross, "Fast Decoding of Multi-Kernel Polar Codes,")
		Pv[0] = std::vector<int>({ 0, 1, 1 });
		Pv[1] = std::vector<int>({ 0, 0, 0, 0, 1, 1, 0, 1, 1 });
		Pv[2] = std::vector<int>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1 });
		
		u_hat = std::vector<B>(N, 9);

		recursive_allocate_nodes_contents(this->polar_tree.get_root(), this->N);
		recursive_initialize_frozen_bits(this->polar_tree.get_root());

		double num_nodes_pre_pruning = polar_tree.num_nodes - 1;
		std::cout << "\n";
		prune_tree(this->polar_tree.get_root());

		std::cout << "\n";
		std::cout << "Nodes Before Pruning: " << num_nodes_pre_pruning << "\n";
		std::cout << "Nodes After Pruning: " << polar_tree.num_nodes << "\n";

		std::cout << "\n";
		std::cout << "Reduction: " << (float)(num_nodes_pre_pruning - polar_tree.num_nodes)/ num_nodes_pre_pruning * 100 << " %\n";
		//std::cout << "% Reduction: " << (double)( (num_nodes_post_pruning * 100.0) / num_nodes_pre_pruning) << " %\n";
	}

	~Polar_Fast_SSC_Decoder()
	{
		this->recursive_deallocate_nodes_contents(this->polar_tree.get_root());
	};

	void _decode_siho(const R* Y_N, B* V_K, bool show_prog_output)
	{

		// Place received real signal into root node
		this->_load(Y_N);

		// Start the decoding
		if (show_prog_output) {
			HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
			SetConsoleTextAttribute(hConsole, 15);

			show_program(this->polar_tree.get_root(), this->N);

			std::cout << "\n";
			SetConsoleTextAttribute(hConsole, 2);
		}
		else
			this->recursive_decode(this->polar_tree.get_root(), this->N);

		for (int i = 0; i < u_hat.size(); i++) {
			if (u_hat[i] == 9 && frozen_set[i] == 0) {
				std::cout << "UNASSIGNED\n";
			}
		}

		this->_store(V_K);
	};


	void _load(const R* Y_N)
	{
		auto* contents = this->polar_tree.get_root()->contents;
		
		for (auto i = 0; i < this->N; i++)
			contents->lambda[i] = Y_N[i];
	};


	void _store(B* V)
	{
		std::vector<Tree_Node<Contents_SC<B, R>>*> l = polar_tree.get_leaves();
		int k = 0;
		for (int i = 0; i < l.size(); i++) {
			if (!frozen_set[i]) {
				V[k++] = u_hat[i];
			}
		}
	}

};

#endif