#ifndef AVLTREE_H_INCLUDED
#define AVLTREE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <inttypes.h>


typedef struct _define_AVLTreeNode {
    uint64_t key;   /** < \note unique global key value */
    uint64_t data;
    uint64_t height;
    struct _define_AVLTreeNode* lChild;
    struct _define_AVLTreeNode* rChild;
} AVLNode;


/**
 * A collection of test in this header file.
 */
void _AVLTreeTestSet();


/*
 * Working functions.
 */

/**
 * Create an AVL tree node.
 *
 * @param key key value of the node
 * @param data data of the node
 * @param lChild left child of the node
 * @param rChild right child of the node
 * @return AVL tree node created
 */
AVLNode* createAVLNode(uint64_t key, uint64_t data, AVLNode* lChild, AVLNode* rChild);

/**
 * Get the height of AVL tree rooting at a node.
 *
 * @param node a node of AVL tree
 * @return height of the AVL tree
 */
uint64_t getHeightOfAVLTree(AVLNode* tree);

/**
 * Insert an AVL tree node into AVL tree.
 *
 * @param root AVL tree root
 * @param node AVL tree node to be inserted
 * @return new tree root of the AVL tree after insertion
 */
AVLNode* insertAVLNode(AVLNode* root, AVLNode* node);

/**
 * Find the AVL tree node according to key value.
 *
 * @param root AVL tree root
 * @param key key value of the node
 * @return AVL tree node if is found matching key; NULL otherwise
 */
AVLNode* findAVLNode(AVLNode* root, uint64_t key);

/** < \note no need to implement the method deleteAVLNode */

/**
 * Clear an AVL tree.
 *
 * @param root AVL tree root
 */
void clearAVLTree(AVLNode* root);


/**
 * Traverse AVL tree and Print it.
 *
 * @param root AVL tree root
 */
void traversePrintAVLTree(AVLNode* root);

/**
 * Traverse AVL tree and find those nodes that have the maximum data in the tree.
 *
 * @param root AVL tree root
 * @param / @return size of the returned array
 * @return an array of AVL tree nodes that have the maximum data in the tree
 */
AVLNode** findNodeswithMaxData(AVLNode* root, uint64_t* nodesNum);







#endif // AVLTREE_H_INCLUDED
