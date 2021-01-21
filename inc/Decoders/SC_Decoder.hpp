#ifndef MULTI_KERNEL_SC_DECODER_HPP
#define MULTI_KERNEL_SC_DECODER_HPP

#include <vector>
#include "Tree.hpp"
#include "Decoder_Functions.hpp"
#include "Contents_SC.hpp"

template <typename B, typename R>
class Polar_SC_Decoder {
private:
	int N;
	int K;
	double Rate;
	std::vector<unsigned int> generator_kernel_order;
	std::vector<unsigned int> frozen_set;

	Tree<Contents_SC<B, R>> polar_tree;
private:
	void recursive_allocate_nodes_contents(Tree_Node<Contents_SC<B, R>>* node_curr, const int vector_size) 
	{
		if (node_curr != nullptr) {
			node_curr->set_contents(new Contents_SC<B, R>(vector_size));

			if (node_curr->get_center() == nullptr) {
				this->recursive_allocate_nodes_contents(node_curr->get_left(), vector_size / 2);
				this->recursive_allocate_nodes_contents(node_curr->get_right(), vector_size / 2);
			}
			else {
				this->recursive_allocate_nodes_contents(node_curr->get_left(), vector_size / 3);
				this->recursive_allocate_nodes_contents(node_curr->get_center(), vector_size / 3);
				this->recursive_allocate_nodes_contents(node_curr->get_right(), vector_size / 3);
			}
		}
	};

	void recursive_initialize_frozen_bits(const Tree_Node<Contents_SC<B, R>>* node_curr) 
	{
		auto* contents = node_curr->get_contents();

		if (!node_curr->is_leaf()) // stop condition
		{
			if (node_curr->get_center() == nullptr) {
				recursive_initialize_frozen_bits(node_curr->get_left());
				recursive_initialize_frozen_bits(node_curr->get_right());
			}
			else {
				recursive_initialize_frozen_bits(node_curr->get_left());
				recursive_initialize_frozen_bits(node_curr->get_center());
				recursive_initialize_frozen_bits(node_curr->get_right());
			}
		}
		else
			contents->is_frozen_bit = frozen_set[node_curr->get_lane_id()];
	};

	void recursive_deallocate_nodes_contents(Tree_Node<Contents_SC<B, R>>* node_curr) 
	{
		if (node_curr != nullptr)
		{
			if (node_curr->get_center() == nullptr) {
				this->recursive_deallocate_nodes_contents(node_curr->get_left());
				this->recursive_deallocate_nodes_contents(node_curr->get_right());
			}
			else {
				this->recursive_deallocate_nodes_contents(node_curr->get_left());
				this->recursive_deallocate_nodes_contents(node_curr->get_center());
				this->recursive_deallocate_nodes_contents(node_curr->get_right());

			}

			auto* contents = node_curr->get_contents();
			delete contents;
			node_curr->set_contents(nullptr);
		}
	};

public:
	Polar_SC_Decoder(int N, int K, std::vector<unsigned int> generator_kernel_order, std::vector<unsigned int> FrozenSet)
		: N(N), K(K), Rate(K / (double)N), generator_kernel_order(generator_kernel_order),
		polar_tree(generator_kernel_order)
	{
		frozen_set = std::vector<unsigned int>(N);
		for (int i = 0; i < FrozenSet.size(); i++) {
			frozen_set[FrozenSet[i]] = 1;
		}

		recursive_allocate_nodes_contents(this->polar_tree.get_root(), this->N);
		recursive_initialize_frozen_bits(this->polar_tree.get_root());
	}

	~Polar_SC_Decoder()
	{
		this->recursive_deallocate_nodes_contents(this->polar_tree.get_root());
	};

	void _load(const R* Y_N) 
	{

		auto* contents = this->polar_tree.get_root()->get_contents();

		for (auto i = 0; i < this->N; i++)
			contents->lambda[i] = Y_N[i];
	};

	void recursive_decode(const Tree_Node<Contents_SC<B, R>>* node_curr) 
	{
		if (!node_curr->is_leaf()) // stop condition
		{
			const auto size = (int)node_curr->get_contents()->lambda.size();

			if (node_curr->get_center() == nullptr) {
				const auto size_2 = size / 2;

				const auto* node_left = node_curr->get_left(); // get left node
				const auto* node_right = node_curr->get_right(); // get right node

				// apply f()
				for (auto i = 0; i < size_2; i++)
					node_left->get_contents()->lambda[i] = F(node_curr->get_contents()->lambda[i],
						node_curr->get_contents()->lambda[size_2 + i]);

				this->recursive_decode(node_left); // recursive call

				// apply g()
				for (auto i = 0; i < size_2; i++)
					node_right->get_contents()->lambda[i] = G(node_curr->get_contents()->lambda[i],
						node_curr->get_contents()->lambda[size_2 + i],
						node_left->get_contents()->s[i]);

				this->recursive_decode(node_right); // recursive call

				// Combine
				for (auto i = 0; i < size_2; i++)
					node_curr->get_contents()->s[i] = node_left->get_contents()->s[i] ^ node_right->get_contents()->s[i]; // bit xor

				for (auto i = 0; i < size_2; i++)
					node_curr->get_contents()->s[size_2 + i] = node_right->get_contents()->s[i]; // bit eq
			}
			else {
				const auto size_3 = size / 3;

				const auto* node_left = node_curr->get_left(); // get left node
				const auto* node_center = node_curr->get_center(); // get center node
				const auto* node_right = node_curr->get_right(); // get right node

				// Apply new F()
				for (auto i = 0; i < size_3; i++)
					node_left->get_contents()->lambda[i] = F_Ternary(node_curr->get_contents()->lambda[i],  // apply f()
						node_curr->get_contents()->lambda[size_3 + i],
						node_curr->get_contents()->lambda[(2 * size_3) + i]);

				this->recursive_decode(node_left); // recursive call

				// Apply G1()
				for (auto i = 0; i < size_3; i++)
					node_center->get_contents()->lambda[i] = G_1(node_curr->get_contents()->lambda[i], // apply g()
						node_curr->get_contents()->lambda[size_3 + i],
						node_curr->get_contents()->lambda[(2*size_3) + i],
						node_left->get_contents()->s[i]);

				this->recursive_decode(node_center);

				// Apply G2()
				for (auto i = 0; i < size_3; i++)
					node_right->get_contents()->lambda[i] = G_2(node_curr->get_contents()->lambda[i], // apply g()
						node_curr->get_contents()->lambda[size_3 + i],
						node_curr->get_contents()->lambda[(2 * size_3) + i],
						node_left->get_contents()->s[i],
						node_center->get_contents()->s[i]);

				this->recursive_decode(node_right);

				// Combine
				for (auto i = 0; i < size_3; i++)
					node_curr->get_contents()->s[i] = node_left->get_contents()->s[i] ^ node_center->get_contents()->s[i]; // bit xor

				for (auto i = 0; i < size_3; i++)
					node_curr->get_contents()->s[size_3 + i] = node_left->get_contents()->s[i] ^ node_right->get_contents()->s[i]; // bit xor

				for (auto i = 0; i < size_3; i++)
					node_curr->get_contents()->s[(2 * size_3) + i] = node_left->get_contents()->s[i] ^ node_center->get_contents()->s[i] ^ node_right->get_contents()->s[i]; // bit xor

			}

			
		}
		else // specific leaf treatment
		{
			node_curr->get_contents()->s[0] = (!node_curr->get_contents()->is_frozen_bit && // if this is a frozen bit then s == 0
				H<B, R>(node_curr->get_contents()->lambda[0])); // apply h()
		}
	};

	void _store(B* V) 
	{
		std::vector<Tree_Node<Contents_SC<B, R>>*> l = polar_tree.get_leaves();
		int k = 0;
		for (int i = 0; i < l.size(); i++) {
			if (!l[i]->get_contents()->is_frozen_bit) {
				V[k++] = l[i]->get_contents()->s[0];
			}
		}
	}

	void _decode_siho(const R* Y_N, B* V_K) 
	{
		// Place received real signal into root node
		this->_load(Y_N);

		// Start the decoding
		this->recursive_decode(this->polar_tree.get_root());

		this->_store(V_K);
	};

	void recursive_store(const Tree_Node<Contents_SC<B, R>>* node_curr, B* V_K, int& k) const 
	{
		auto* contents = node_curr->get_contents();

		if (!node_curr->is_leaf()) // stop condition
		{
			this->recursive_store(node_curr->get_left(), V_K, k); // recursive call
			this->recursive_store(node_curr->get_right(), V_K, k); // recursive call
		}
		else
			if (!frozen_set[node_curr->get_lane_id()])
				V_K[k++] = contents->s[0];
	};
};

#endif