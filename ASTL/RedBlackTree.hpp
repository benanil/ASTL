


class RedBlackTreeAI {
public:
    enum class Color { RED, BLACK };

    struct Node {
        int key;
        Color color;
        Node* left;
        Node* right;
        Node* parent;

        Node(int key) : key(key), color(Color::RED), left(nullptr), right(nullptr), parent(nullptr) {}
    };

    Node* root;

    void DestroyRec(Node* node)
    {
        if (!node) return;
        DestroyRec(node->left);
        DestroyRec(node->right);
        delete node;
    }

    ~RedBlackTreeAI()
    {
        DestroyRec(root);
    }

    void rotateLeft(Node* node) {
        Node* rightChild = node->right;
        node->right = rightChild->left;

        if (rightChild->left != nullptr)
            rightChild->left->parent = node;

        rightChild->parent = node->parent;

        if (node->parent == nullptr)
            root = rightChild;
        else if (node == node->parent->left)
            node->parent->left = rightChild;
        else
            node->parent->right = rightChild;

        rightChild->left = node;
        node->parent = rightChild;
    }

    void rotateRight(Node* node) {
        Node* leftChild = node->left;
        node->left = leftChild->right;

        if (leftChild->right != nullptr)
            leftChild->right->parent = node;

        leftChild->parent = node->parent;

        if (node->parent == nullptr)
            root = leftChild;
        else if (node == node->parent->left)
            node->parent->left = leftChild;
        else
            node->parent->right = leftChild;

        leftChild->right = node;
        node->parent = leftChild;
    }

    void removeFixup(Node* node) {
        
        while (node && node != root && node->color == Color::BLACK) {
            if (node == node->parent->left) {
                Node* sibling = node->parent->right;

                if (sibling->color == Color::RED) {
                    sibling->color = Color::BLACK;
                    node->parent->color = Color::RED;
                    rotateLeft(node->parent);
                    sibling = node->parent->right;
                }

                if (sibling->left->color == Color::BLACK && sibling->right->color == Color::BLACK) {
                    sibling->color = Color::RED;
                    node = node->parent;
                }
                else {
                    if (sibling->right->color == Color::BLACK) {
                        sibling->left->color = Color::BLACK;
                        sibling->color = Color::RED;
                        rotateRight(sibling);
                        sibling = node->parent->right;
                    }

                    sibling->color = node->parent->color;
                    node->parent->color = Color::BLACK;
                    sibling->right->color = Color::BLACK;
                    rotateLeft(node->parent);
                    node = root;
                }
            }
            else {
                Node* sibling = node->parent->left;

                if (sibling->color == Color::RED) {
                    sibling->color = Color::BLACK;
                    node->parent->color = Color::RED;
                    rotateRight(node->parent);
                    sibling = node->parent->left;
                }

                if (sibling->right->color == Color::BLACK && sibling->left->color == Color::BLACK) {
                    sibling->color = Color::RED;
                    node = node->parent;
                }
                else {
                    if (sibling->left->color == Color::BLACK) {
                        sibling->right->color = Color::BLACK;
                        sibling->color = Color::RED;
                        rotateLeft(sibling);
                        sibling = node->parent->left;
                    }

                    sibling->color = node->parent->color;
                    node->parent->color = Color::BLACK;
                    sibling->left->color = Color::BLACK;
                    rotateRight(node->parent);
                    node = root;
                }
            }
        }

        if (node) node->color = Color::BLACK;
    }


    Node* getSuccessor(Node* node) {
        if (node->right != nullptr) {
            node = node->right;
            while (node->left != nullptr)
                node = node->left;
            return node;
        }

        Node* parent = node->parent;
        while (parent != nullptr && node == parent->right) {
            node = parent;
            parent = parent->parent;
        }
        return parent;
    }

    Node* searchNode(int key) {
        Node* current = root;

        while (current != nullptr) {
            if (key == current->key)
                return current;
            else if (key < current->key)
                current = current->left;
            else
                current = current->right;
        }

        return nullptr; // Node not found
    }

    void remove(int key) {
        Node* node = searchNode(key);

        if (node == nullptr)
            return;

        Node* replacement = nullptr;

        if (node->left == nullptr || node->right == nullptr)
            replacement = node;
        else
            replacement = getSuccessor(node);

        Node* child = nullptr;
        if (replacement->left != nullptr)
            child = replacement->left;
        else
            child = replacement->right;

        if (child != nullptr)
            child->parent = replacement->parent;

        if (replacement->parent == nullptr)
            root = child;
        else if (replacement == replacement->parent->left)
            replacement->parent->left = child;
        else
            replacement->parent->right = child;

        if (replacement != node)
            node->key = replacement->key;

        if (replacement->color == Color::BLACK)
            removeFixup(child);

        delete replacement;
    }

    void insertFixup(Node* node) {
        while (node->parent != nullptr && node->parent->color == Color::RED) {
            if (node->parent == node->parent->parent->left) {
                Node* uncle = node->parent->parent->right;

                if (uncle != nullptr && uncle->color == Color::RED) {
                    node->parent->color = Color::BLACK;
                    uncle->color = Color::BLACK;
                    node->parent->parent->color = Color::RED;
                    node = node->parent->parent;
                }
                else {
                    if (node == node->parent->right) {
                        node = node->parent;
                        rotateLeft(node);
                    }

                    node->parent->color = Color::BLACK;
                    node->parent->parent->color = Color::RED;
                    rotateRight(node->parent->parent);
                }
            }
            else {
                Node* uncle = node->parent->parent->left;

                if (uncle != nullptr && uncle->color == Color::RED) {
                    node->parent->color = Color::BLACK;
                    uncle->color = Color::BLACK;
                    node->parent->parent->color = Color::RED;
                    node = node->parent->parent;
                }
                else {
                    if (node == node->parent->left) {
                        node = node->parent;
                        rotateRight(node);
                    }

                    node->parent->color = Color::BLACK;
                    node->parent->parent->color = Color::RED;
                    rotateLeft(node->parent->parent);
                }
            }
        }

        root->color = Color::BLACK;
    }

    RedBlackTreeAI() : root(nullptr) {}

    void insert(int key) {
        Node* newNode = new Node(key);

        Node* parent = nullptr;
        Node* current = root;

        while (current != nullptr) {
            parent = current;

            if (key < current->key)
                current = current->left;
            else
                current = current->right;
        }

        newNode->parent = parent;

        if (parent == nullptr)
            root = newNode;
        else if (key < parent->key)
            parent->left = newNode;
        else
            parent->right = newNode;

        insertFixup(newNode);
    }
};