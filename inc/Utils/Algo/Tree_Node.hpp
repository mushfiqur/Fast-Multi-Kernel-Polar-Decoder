#ifndef ALGO_TREE_TREE_NODE_HPP
#define ALGO_TREE_TREE_NODE_HPP

/* 
	* Derived from aff3ct (https://github.com/aff3ct/aff3ct/)
	* Combines binary_node.hpp and binary_node.hxx
	* Added 'center' node for ternary stages
	* Renamed class to be 'Tree_Node' instead of 'Binary_Node' since
	*	this is no longer a binary tree
 */

enum class NODE_TYPE {
	RATE_N,
	RATE_0,
	RATE_1,
	SPC,
	REP2,
	REP3A,
	REP3B,
	REP3C
};

template<typename T>
class Tree;

template<typename T>
class Tree_Node {
	friend Tree<T>;
public:
	Tree_Node<T>* father; /*!< Pointer to the father node       (nullptr if this node is the root).						  */
	Tree_Node<T>* left;   /*!< Pointer to the left child node   (nullptr if this node is a leaf).						  */
	Tree_Node<T>* right;  /*!< Pointer to the right child node  (nullptr if this node is a leaf).						  */
	Tree_Node<T>* center; /*!< Pointer to the center child node (nullptr if this node is a leaf, or doesn't exist).    */

	T* contents;		  /*!< Pointer to the node contents, could be anything. */

	const int depth;	  /*!< Depth   of this node (vertical   indexing). */
	const int lane_id;    /*!< Lane id of this node (horizontal indexing). */

	NODE_TYPE type;
public:
	Tree_Node(Tree_Node<T>* father_node,
		Tree_Node<T>* left_node,
		Tree_Node<T>* right_node,
		Tree_Node<T>* center_node,
		T* contents,
		int depth,
		int lane_id)
		: father(father_node), left(left_node), right(right_node), center(center_node),
		contents(contents), depth(depth), lane_id(lane_id), type(NODE_TYPE::RATE_N)
	{
	}

	bool Tree_Node<T>::is_root() const
	{
		return (this->father == nullptr);
	}

	bool Tree_Node<T>::is_leaf() const
	{
		return (this->left == nullptr && this->right == nullptr);
	}

	bool Tree_Node<T>::is_empty() const
	{
		return (this->contents == nullptr);
	}

	bool Tree_Node<T>::is_right() const
	{
		return ((!this->is_root()) && (this->father->get_right() == this));
	}

	bool Tree_Node<T>::is_center() const
	{
		return ((!this->is_root()) && (this->father->get_center() == this));
	}

	bool Tree_Node<T>::is_left() const
	{
		return ((!this->is_root()) && (this->father->get_left() == this));
	}

	Tree_Node<T>* Tree_Node<T>::get_father() const
	{
		return this->father;
	}

	Tree_Node<T>* Tree_Node<T>::get_left() const
	{
		return this->left;
	}

	Tree_Node<T>* Tree_Node<T>::get_center() const
	{
		return this->center;
	}

	Tree_Node<T>* Tree_Node<T>::get_right() const
	{
		return this->right;
	}

	T* Tree_Node<T>::get_contents() const
	{
		return this->contents;
	}

	void Tree_Node<T>::set_contents(T* contents)
	{
		this->contents = contents;
	}

	int Tree_Node<T>::get_depth() const
	{
		return this->depth;
	}

	int Tree_Node<T>::get_lane_id() const
	{
		return this->lane_id;
	}

	void Tree_Node<T>::cut_left()
	{
		delete this->left;
		this->left = nullptr;
	}

	void Tree_Node<T>::cut_right()
	{
		delete this->right;
		this->right = nullptr;
	}
};

#endif