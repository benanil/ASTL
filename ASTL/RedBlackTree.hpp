// use this when you don't want to use hash functions 
// and if you want sorted map you can use this otherwise use HashMap or hashset
// this algorithm is from Introduction to algorithms book by Thomas H.Cormen and Ronald L.Rivest

#pragma once

#include "Memory.hpp"

template<typename ValueT>
class RedBlackTree
{
public:

	typedef void (*TraverseFunc)(ValueT& val);

	struct Node
	{
		Node() {}

		// protector constructor
		Node(Node* protect) 
		: parent(protect), left(protect), right(protect) 
		{ }

        // set as red by default		
		Node(ValueT val, Node* parent_, Node* protect)
		: parent((Node*)(0x1ull | (uint64)parent_)), left(protect), right(protect), value((ValueT&&)val)
		{  }
		
		Node*  parent;
		Node*  left;
		Node*  right;
		ValueT value;

		bool IsRed()  const { return ((uint64_t)parent) & 0x1ull; }
		bool IsBlack() const { return !IsRed(); }

		void SetRed()   { parent = (Node*)((uint64_t)parent | 0x1ull); }  
		void SetBlack() { parent = (Node*)((uint64_t)parent & ~0x1ull); }

		bool GetColor() { return IsRed(); }
		void SetColor(bool isRed)
		{
			parent = (Node*)(uint64_t(isRed) | ((uint64_t)parent & ~0x1ull));
		}

		void SetParent(Node* newParent) {
			parent = (Node*)((uint64_t)newParent | (((uint64_t)parent) & 0x1));
		}

		Node* GetParent() const {
			return (Node*)(((uint64_t)parent) & ~0x1);
		}
	};

    Node* m_root = nullptr;
	FixedSizeGrowableAllocator<Node> m_allocator{};
	static Node m_protect;

	RedBlackTree()
	{
		// m_protect.SetBlack(); // by default this is black
		m_protect.left = m_protect.right = &m_protect;
		m_protect.SetParent(&m_protect);
		m_root = &m_protect;
	}

    ~RedBlackTree()
	{
		Clear();
	}

	const Node* Nil() const { return &m_protect; }

	template<typename...Args>
	Node* Insert(Args&&... args)
	{
		Node* node   = m_root;
		Node* parent = &m_protect;
		ValueT value(Forward<Args>(args)...);

		while (node != &m_protect)
		{
			parent = node;
			if (node->value < value) node = node->right;
			else if (value < node->value) node = node->left;
			else return node;
		}

		Node* added = AllocateNode(Forward<ValueT>(value), parent);
		
		if (parent != &m_protect)
			if (value < parent->value)
				parent->left = added;
			else
				parent->right = added;
		else 
			m_root = added;

		InsertFixup(added);
		Validate();
		return added;
	}

	Node* Insert(const ValueT& value)
	{
		Node* node   = m_root;
		Node* parent = &m_protect;

		while (node != &m_protect)
		{
			parent = node;
			if (node->value < value) node = node->right;
			else if (value < node->value) node = node->left;
			else return node;
		}

		Node* added = AllocateNode(value, parent);
		
        if (parent != &m_protect)
			if (value < parent->value)
				parent->left = added;
			else
				parent->right = added;
		else 
			m_root = added;

		InsertFixup(added);
		Validate();
		return added;
	}

	int Erase(const ValueT& key)
	{
		if (Node* toErase = FindNode(key))
		{
			Erase(toErase);
			return true;
		}
		return false;
	}

	void Erase(Node* node)
	{
		ASSERT(m_root); // tree is already empty
		Node* replacement;

		if (node->left == &m_protect || node->right == &m_protect)
			replacement = node;
		else
		{
			// get successor
			replacement = node->right;
			while (replacement->left != &m_protect)
				replacement = replacement->left;
		}

		Node* eraseChild = replacement->left != &m_protect ? replacement->left : replacement->right;
		eraseChild->SetParent(replacement->GetParent());
		
		if (replacement->GetParent() != &m_protect)
		{
			if (replacement == replacement->GetParent()->left)
				replacement->GetParent()->left = eraseChild;
			else
				replacement->GetParent()->right = eraseChild;
		}
		else
			m_root = eraseChild;

		node->value = replacement->value;

		if (replacement->IsBlack())
			FixupAfterErase(eraseChild);

		FreeNode(replacement, false);

		Validate();
	}

    Node* FindNode(const ValueT& key)
	{
		Node* node = m_root;
		while (node != &m_protect)
		{
			if (node->value < key)
				node = node->right;
			else if (key < node->value)
				node = node->left;
			else 
				return node;
		}
		return nullptr;
	}

	bool Contains(const ValueT& key)
	{
		return FindNode(key) != nullptr;
	}

	void Clear()
	{
		if (m_root)
		{
			FreeNode(m_root, true);
			m_root = &m_protect;
		}
	}

	void TraverseNode(Node* n, TraverseFunc func)
	{
		if (n->left != &m_protect)
			TraverseNode(n->left, func);
		func(n->value);
		if (n->right != &m_protect)
			TraverseNode(n->right, func);
	}

	void Traverse(TraverseFunc func)
	{
		if (m_root)
			TraverseNode(m_root, func);
	}

	Node* GetBeginNode() const
	{
		Node* iter = &m_protect;
		if (m_root != &m_protect)
		{
			iter = m_root;
			while (iter->left != &m_protect)
				iter = iter->left;
		}
		return iter;
	}

	Node* FindNextNode(Node* node) const
	{
		if (node == &m_protect) return &m_protect;
		Node* next = &m_protect;
		
		if (node->right != &m_protect)
		{
			next = node->right;
			while (next->left != &m_protect)
				next = next->left;
		}
		else if (node->GetParent() != &m_protect)
		{
			if (node == node->GetParent()->left)
				return node->GetParent();
			else
			{
				next = node;
				while (next->GetParent() != &m_protect)
				{
					if (next == next->GetParent()->right)
						next = next->GetParent();
					else
						return next->GetParent();
				}
				next = nullptr;
			}
		}
		else
		{
			ASSERT(node == m_root);
		}
		return next;
	}

	int NumNodes(const Node* n) const
	{
        return n == &m_protect ? nullptr : 1 + NumNodes(n->left) + NumNodes(n->right);
	}

	void InsertFixup(Node* new_node)
	{
		ASSERT(new_node->IsRed());
		Node* iter = new_node;
		
		while (iter->GetParent()->IsRed())
		{
			Node* grandparent = iter->GetParent()->GetParent();
			if (iter->GetParent() == grandparent->left)
			{
				Node* uncle = grandparent->right;
				// Both parent and uncle are red.
				// Repaint both, make grandparent red.
				if (uncle->IsRed())
				{
					iter->GetParent()->SetBlack();
					uncle->SetBlack();
					grandparent->SetRed();
					iter = grandparent;
				}
				else
				{
					if (iter == iter->GetParent()->right)
					{
						iter = iter->GetParent();
						RotateLeft(iter);
					}
					grandparent = iter->GetParent()->GetParent();
					iter->GetParent()->SetBlack();
					grandparent->SetRed();
					RotateRight(grandparent);
				}
			}
			else
			{
				Node* uncle = grandparent->left;
				if (uncle->IsRed())
				{
					grandparent->SetRed();
					iter->GetParent()->SetBlack();
					uncle->SetBlack();
					iter = grandparent;
				}
				else
				{
					if (iter == iter->GetParent()->left)
					{
						iter = iter->GetParent();
						RotateRight(iter);
					}
					grandparent = iter->GetParent()->GetParent();
					iter->GetParent()->SetBlack();
					grandparent->SetRed();
					RotateLeft(grandparent);
				}
			}
		}
		m_root->SetBlack();
	}

	void FixupAfterErase(Node* n)
	{
		Node* iter = n;
		while (iter != m_root && iter->IsBlack())
		{
			if (iter == iter->GetParent()->left)
			{
				Node* sibling = iter->GetParent()->right;
				if (sibling->IsRed())
				{
					sibling->SetBlack();
					iter->GetParent()->SetRed();
					RotateLeft(iter->GetParent());
					sibling = iter->GetParent()->right;
				}
				if (sibling->left->IsBlack() && sibling->right->IsBlack())
				{
					sibling->SetRed();
					iter = iter->GetParent();
				}
				else
				{
					if (sibling->right->IsBlack())
					{
						sibling->left->SetBlack();
						sibling->SetRed();
						RotateRight(sibling);
						sibling = iter->GetParent()->right;
					}
					sibling->SetColor(iter->GetParent()->GetColor());
					iter->GetParent()->SetBlack();
					sibling->right->SetBlack();
					RotateLeft(iter->GetParent());
					iter = m_root;
				}
			}
			else // iter == right child
			{
				Node* sibling = iter->GetParent()->left;
				if (sibling->IsRed())
				{
					sibling->SetBlack();
					iter->GetParent()->SetRed();
					RotateRight(iter->GetParent());
					sibling = iter->GetParent()->left;
				}

				if (sibling->left->IsBlack() && sibling->right->IsBlack())
				{
					sibling->SetRed();
					iter = iter->GetParent();
				}
				else
				{
					if (sibling->left->IsBlack())
					{
						sibling->right->SetBlack();
						sibling->SetRed();
						RotateLeft(sibling);
						sibling = iter->GetParent()->left;
					}
					sibling->SetColor(iter->GetParent()->GetColor());
					iter->GetParent()->SetBlack();
					sibling->left->SetBlack();
					RotateRight(iter->GetParent());
					iter = m_root;
				}
			}
		}
		iter->SetBlack();
	}

#if 1 && defined(_DEBUG)
    void ValidateNodeRec(Node* n) const
	{
		if (m_root == nullptr || m_root == &m_protect) return;
		ASSERT(n->GetParent() == &m_protect ||
			n->GetParent()->left == n || n->GetParent()->right == n);
		
		if (n->IsRed())
		{
			ASSERT(n->left->IsBlack());
			ASSERT(n->right->IsBlack());
		}
		if (n->left != &m_protect)
			ValidateNodeRec(n->left);
		if (n->right != &m_protect)
			ValidateNodeRec(n->right);
	}

	void Validate() const
	{
		ASSERT(m_root->IsBlack());
		ValidateNodeRec(m_root);
	}
#else
    void Validate() const { }
#endif

	void RotateLeft(Node* node)
	{
		Node* rightChild = node->right;
		node->right = rightChild->left;
		
		if (node->right != &m_protect)
			node->right->SetParent(node);

		rightChild->SetParent(node->GetParent());
		
		if (node->GetParent() == &m_protect)
			m_root = rightChild;
		else
		{
			if (node == node->GetParent()->left)
				node->GetParent()->left = rightChild;
			else
				node->GetParent()->right = rightChild;
		}
		rightChild->left = node;
		node->SetParent(rightChild);
	}

	void RotateRight(Node* node)
	{
		Node* leftChild = node->left;
		node->left = leftChild->right;

		if (node->left != &m_protect)
			node->left->SetParent(node);

		leftChild->SetParent(node->GetParent());
		if (node->GetParent() == &m_protect)
			m_root = leftChild;
		else
			if (node == node->GetParent()->left)
				node->GetParent()->left = leftChild;
			else
				node->GetParent()->right = leftChild;
		
		leftChild->right = node;
		node->SetParent(leftChild);
	}

	Node* AllocateNode(const ValueT& value, Node* parent)
	{
		return m_allocator.Allocate(value, parent, &m_protect);
	}
	
	template<typename... Args>
	Node* AllocateNode(Node* parent, Args&&... args)
	{
		ValueT val(Forward<Args>(args)...);
		return m_allocator.Allocate(Forward<ValueT>(val), parent, &m_protect);
	}

	void FreeNode(Node* n, bool recursive)
	{
		if (recursive)
		{
			if (n->left != &m_protect)
				FreeNode(n->left, true);
			if (n->right != &m_protect)
				FreeNode(n->right, true);
		}
		if (n != &m_protect)
		{
			m_allocator.Deallocate(n, 1);
		}
	}
};

template<typename ValueT>
typename RedBlackTree<ValueT>::Node RedBlackTree<ValueT>::m_protect(&m_protect);

template<typename ValueT>
using Set = RedBlackTree<ValueT>;

template<typename KeyT, typename ValueT>
using Map = RedBlackTree<KeyValuePair<KeyT, ValueT>>;
