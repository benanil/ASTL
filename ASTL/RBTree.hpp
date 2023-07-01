// an red black tree, 
// this algorithm is from Introduction to algorithms book by Thomas H.Cormen and Ronald L.Rivest
// child[0] is left, child[1] is right
// todo key value map
// todo custom allocator (better allocations)

#pragma once

#include "Algorithms.hpp"

typedef int ValueT;

class RedBlackTree 
{
public:
    struct Node 
    {
        Node* child[2]{};
        Node*  parent = nullptr;
        ValueT data;
        bool   isRed  : 1;

        Node(ValueT val) : data((ValueT&&)val), isRed(false) {  }
    };

    Node* root = nullptr;

public:
    
    RedBlackTree() : root(nullptr)
    { }

    ~RedBlackTree()
    {
        DeleteRec(root);
    }

    // copy constructor
    RedBlackTree(const RedBlackTree& other) : root(other.root)
    { }
    
    RedBlackTree(RedBlackTree&& other) noexcept // move constructor 
    { root = other.root; other.root = nullptr; }
    
    void Insert(const ValueT& val)
    {
        Node* z = new Node(val);
        Node* x = root;
        Node* y = nullptr;
        
        while (x)
        {
            y = x;
            if (z->data == x->data) return; // already exist
            x = x->child[z->data > x->data];
        }

        z->parent = y;
        
        if (y == nullptr)
            root = z;
        else
            y->child[z->data > y->data] = z;
        
        z->isRed = true;

        InsertFixup(z);
    }

    void Remove(const ValueT& val)
    {
        Node* node = FindNode(val);
        if (!node) return;

        Node* replacement = nullptr;

        if (node->child[0] == nullptr || node->child[1] == nullptr)
            replacement = node;
        else
            replacement = GetSuccessor(node);

        Node* child = nullptr;
        child = replacement->child[replacement->child[0] == nullptr];

        if (child != nullptr)
            child->parent = replacement->parent;

        if (replacement->parent == nullptr)
            root = child;
        else
        {
            bool isLeftChild = replacement != replacement->parent->child[0];
            replacement->parent->child[isLeftChild] = child;
        }

        if (replacement != node)
            node->data = replacement->data;

        if (replacement->isRed == 0)
            DeleteFixup(child);

        delete replacement;
    }

    bool Contains(const ValueT& val)
    {
        return FindNode(val) != nullptr;
    }

private:

    Node* Minimum(Node* n)
    {
        Node* x = nullptr;
        while (n)
        {
            x = n;
            n = n->child[0];
        }

        if (x->child[1])
        {
            x = x->child[1];
        }
        return x;
    }

    
    Node* GetSuccessor(Node* node) {
        if (node->child[0] != nullptr) {
            return Minimum(node);
        }

        Node* parent = node->parent;
        while (parent != nullptr && node == parent->child[1]) {
            node = parent;
            parent = parent->parent;
        }
        return parent;
    }

    Node* FindNode(const ValueT& val)
    {
        Node* n = root;
        while (n)
        {
            if (n->data == val) return n;
            n = n->child[val > n->data];
        }
        return nullptr;
    }

    void DeleteFixup(Node* x)
    {
        while (x && x != root && !x->isRed)
        {
            bool leftChild = x == x->parent->child[0];
            Node* w = x->parent->child[leftChild];
            
            if (w->isRed) // case 1
            {
                w->isRed = false;
                x->parent->isRed = true;
                Rotate(x->parent, !leftChild);
                w = x->parent->child[leftChild];
            }
    
            bool leftBlack  = w->child[!leftChild] != nullptr && !w->child[!leftChild]->isRed;
            bool rightBlack = w->child[ leftChild] != nullptr && !w->child[ leftChild]->isRed;
    
            if (leftBlack && rightBlack) // case 2 
            {
                w->isRed = true;
                x = x->parent;
            }
            else if (rightBlack) // case 3
            {
                w->child[!leftChild]->isRed = false;
                w->isRed = true;
                Rotate(w, leftChild);
                w = x->parent->child[leftChild];
            }
            
            w->isRed = x->parent->isRed; // case 4
            x->parent->isRed = false;
            w->child[leftChild]->isRed = false;
            Rotate(x->parent, !leftChild);
            x = root;
        }
    
        if (x) x->isRed = false;
    }

    void InsertFixup(Node* z)
    {
        if (z->parent == nullptr || z->parent->parent == nullptr)
        {
            root->isRed = false;
            return;
        }

        while (z->parent && z->parent->isRed)
        {
            Node* grandParent = z->parent->parent;
            bool isLeft       = z->parent == grandParent->child[0];
            Node* y           = grandParent->child[isLeft];
            
            if (y && y->isRed) // case 1
            {
                z->parent->isRed   = false;
                y->isRed           = false;
                grandParent->isRed = true;
                z = grandParent;
            }
            else 
            {
                if (z == z->parent->child[isLeft]) // case 2
                    z = z->parent,
                    Rotate(z, !isLeft);

                z->parent->isRed   = false; // case 3
                grandParent->isRed = true;
                Rotate(grandParent, isLeft);
            }
        }
        root->isRed = false;
    }

    void DeleteRec(Node* node)
    {
        if (!node) return;
        DeleteRec(node->child[0]);
        DeleteRec(node->child[1]);
        delete node;
    }

    void Rotate(Node* x, int dir)
    {
        Node* y = x->child[1-dir];
        x->child[1-dir] = y->child[dir];
        
        if (y->child[dir] != nullptr)
            y->child[dir]->parent = x;

        y->parent = x->parent;
        
        if (x->parent == nullptr)
            root = y;
        else
            x->parent->child[x != x->parent->child[dir]] = y;
        
        y->child[dir] = x;
        x->parent     = y;
    }
};
