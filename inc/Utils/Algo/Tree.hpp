#ifndef ALGO_TREE_TREE_HPP
#define ALGO_TREE_TREE_HPP

/*
	* Derived from aff3ct (https://github.com/aff3ct/aff3ct/)
	* Combines Tree.hpp and Tree.hxx
	* Renamed class to be 'Tree' instead of 'Tree' since
	*	this is no longer a binary tree
 */

#include <vector>
#include "Tree_Node.hpp"

template<typename T>
class Tree {
private:
	int             depth;               /*!< Depth of the tree. */
	Tree_Node<T>*   root;                /*!< Pointer to the root node (first node at the top in the tree). */
	std::vector<Tree_Node<T>*> leaves;	 /*!< Vector of the tree leave pointers. */

	std::vector<unsigned int> generator_kernel_order;
private:
	void Tree<T>::create_nodes(Tree_Node<T>* cur_node, int cur_depth, std::vector<int>& lanes)
	{
		if (cur_depth < this->depth)
		{
			if (generator_kernel_order[cur_depth-1] == 2) {
				cur_node->left = new Tree_Node<T>(cur_node, nullptr, nullptr, nullptr, nullptr, cur_depth, lanes[cur_depth]++);
				cur_node->right = new Tree_Node<T>(cur_node, nullptr, nullptr, nullptr, nullptr, cur_depth, lanes[cur_depth]++);
				num_nodes += 2;

				this->create_nodes(cur_node->left, cur_depth + 1, lanes);
				this->create_nodes(cur_node->right, cur_depth + 1, lanes);
			}
			else if (generator_kernel_order[cur_depth-1] == 3) {
				cur_node->left = new Tree_Node<T>(cur_node, nullptr, nullptr, nullptr, nullptr, cur_depth, lanes[cur_depth]++);
				cur_node->center = new Tree_Node<T>(cur_node, nullptr, nullptr, nullptr, nullptr, cur_depth, lanes[cur_depth]++);
				cur_node->right = new Tree_Node<T>(cur_node, nullptr, nullptr, nullptr, nullptr, cur_depth, lanes[cur_depth]++);
				num_nodes += 3;

				this->create_nodes(cur_node->left, cur_depth + 1, lanes);
				this->create_nodes(cur_node->center, cur_depth + 1, lanes);
				this->create_nodes(cur_node->right, cur_depth + 1, lanes);
			}

		}
	}


	void Tree<T>::recursive_get_leaves(Tree_Node<T>* cur_node) {
		if (cur_node->is_leaf())
			leaves.push_back(cur_node);
		else
		{
			if(cur_node->left != nullptr)
				recursive_get_leaves(cur_node->left);
			if (cur_node->center != nullptr)
				recursive_get_leaves(cur_node->center);
			if (cur_node->right != nullptr)
				recursive_get_leaves(cur_node->right);
		}
	}

	void Tree<T>::print_tree(Tree_Node<T>* r)
	{
		std::cout << "(" << r->depth << ", " << r->lane_id << ")\n";
		if (!r->is_leaf()) {
			if (r->left != nullptr) {
				std::cout << "--left--";
				print_tree(r->left);
			}
			std::cout << "(" << r->depth << ", " << r->lane_id << ")\n";
			if (r->center != nullptr) {
				std::cout << "--center--";
				print_tree(r->center);
			}
			std::cout << "(" << r->depth << ", " << r->lane_id << ")\n";
			if (r->right != nullptr) {
				std::cout << "--right--";
				print_tree(r->right);
			}
		}
		else {
			std::cout << "(" << r->father->depth << ", " << r->father->lane_id << ")\n";
		}
	}
public:
	int num_nodes;
public:
	Tree<T>::Tree(std::vector<unsigned int> generator_kernel_order)
		: depth(generator_kernel_order.size() + 1), generator_kernel_order(generator_kernel_order), num_nodes(0)
	{
		std::vector<int> lanes(depth);

		for (unsigned i = 0; i < lanes.size(); i++)
			lanes[i] = 0;

		this->root = new Tree_Node<T>(nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0);
		num_nodes++;
		this->create_nodes(this->root, 1, lanes);

		recursive_get_leaves(this->root);
	}

	void Tree<T>::delete_nodes(Tree_Node<T>* cur_node)
	{
		if (cur_node != nullptr)
		{
			if (cur_node->left != nullptr)
				this->delete_nodes(cur_node->left);

			if (cur_node->center != nullptr)
				this->delete_nodes(cur_node->center);

			if (cur_node->right != nullptr)
				this->delete_nodes(cur_node->right);

			delete cur_node;
			num_nodes--;
		}
	}

	std::vector<Tree_Node<T>*> Tree<T>::get_leaves() const {
		return this->leaves;
	}

	Tree_Node<T>* Tree<T>::get_root() {
		return this->root;
	}

	Tree<T>::~Tree()
	{
		this->delete_nodes(this->root);
	}

	void Tree<T>::print() {
		this->print_tree(this->root);
	}

	
};

#endif