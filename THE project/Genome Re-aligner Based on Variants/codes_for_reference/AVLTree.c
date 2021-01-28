#include "AVLTree.h"

#include "AuxiliaryFunction.h"
#include "Queue.h"

static AVLNode* rotateAVL_LL(AVLNode* originalRoot);
static AVLNode* rotateAVL_RR(AVLNode* originalRoot);
static AVLNode* rotateAVL_LR(AVLNode* originalRoot);
static AVLNode* rotateAVL_RL(AVLNode* originalRoot);






static void _AVLTreeTest();
static void _usageOfAVLinSNAPTest();

/**
 * A collection of test in this header file.
 */
void _AVLTreeTestSet() {
    _AVLTreeTest();
    _usageOfAVLinSNAPTest();
}

/*
 * Tests for working functions.
 */

/**
 * Test the basic functions of AVL tree.
 */
static void _AVLTreeTest() {
    printf("\n**************** _AVLTreeTest ****************\n");
    AVLNode* tree = NULL;
    AVLNode* node = NULL;
    uint64_t key = 4;
    uint64_t totalNodesNum = 7;

    for(uint64_t i = 0; i < totalNodesNum; i++) {
        tree = insertAVLNode(tree, createAVLNode(i, rand() % 10, NULL, NULL));
    }
    traversePrintAVLTree(tree);

    key = 5;
//    key = 11;
    node = findAVLNode(tree, key);
    if(node == NULL) {
        printf("key value of node being %"PRIu64" not found. \n", key);
    } else {
        printf("find node with key %"PRIu64" -> data %"PRIu64"\n", node->key, node->data);
    }

    clearAVLTree(tree);
}

/**
 * Test the actual usage in SNAP.
 */
static void _usageOfAVLinSNAPTest() {
    printf("\n**************** _usageOfAVLinSNAPTest ****************\n");
    uint64_t positionCount = 10;
    AVLNode* tree = NULL;
    AVLNode* node = NULL;
    uint64_t* positions = (uint64_t*)malloc(sizeof(uint64_t) * positionCount);
    uint64_t* seedCounts = (uint64_t*)malloc(sizeof(uint64_t) * positionCount);
    if(seedCounts == NULL || positions == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }

    printf("\n>> add hits of the first seed\n");
    for(uint64_t i = 0; i < positionCount; i++) {
        positions[i] = i;
        seedCounts[i] = 1 + rand() % 5;    /** < \note at least >= 1 */
        printf("position:%"PRIu64", seedCount:%"PRIu64"\n", positions[i], seedCounts[i]);
        tree = insertAVLNode(tree, createAVLNode(positions[i], 1, NULL, NULL));
        seedCounts[i]--;
    }
    printf("> traverse and print:\n");
    traversePrintAVLTree(tree);

    printf("\n>> add hits of the following seeds\n");
    for(uint64_t i = 0; i < positionCount; i++) {
        while(seedCounts[i] > 0) {  /** < \note seedCounts[] is an unsigned array! */
            node = findAVLNode(tree, positions[i]);
            node->data = node->data + 1;
            seedCounts[i] = seedCounts[i] - 1;  /** < \note stop when seedCount[i] == 0 */
        }
    }
    node = findAVLNode(tree, positionCount * 2);
    printf("Not found position (%"PRIu64") -> returned AVL node: ", positionCount * 2);
    if(node == NULL) { printf("NULL");}
    else {printf("(address): %p", node);}
    printf("\n");
    printf("> traverse and print:\n");
    traversePrintAVLTree(tree);


    uint64_t nodesNum = 0;
    AVLNode** nodes = findNodeswithMaxData(tree, &nodesNum);
    printf("\n>> positions that have the maximum seedCount - total number: %"PRIu64"\n", nodesNum);
    for(uint64_t i = 0; i < nodesNum; i++) {
        printf("position: %"PRIu64", seedCount: %"PRIu64"\n", nodes[i]->key, nodes[i]->data);
    }

    clearAVLTree(tree);
    free(positions);
    free(seedCounts);
}



/*
 * Working functions.
 */

AVLNode* createAVLNode(uint64_t key, uint64_t data, AVLNode* lChild, AVLNode* rChild) {
    AVLNode* node = (AVLNode*)malloc(sizeof(AVLNode));
    if(node == NULL) {
        printf("ERROR: memory not enough when creating an AVL node. \n");
        exit(EXIT_FAILURE);
    }
    node->key = key;
    node->data = data;
    node->height = 0;
    node->lChild = lChild;
    node->rChild = rChild;
    return node;
}

uint64_t getHeightOfAVLTree(AVLNode* tree) {
    if(tree == NULL) {
//        printf("ERROR: null pointer occurs when getting height of an AVL tree. \n");
//        exit(EXIT_FAILURE);
        return 0;
    }
    return tree->height;
}

AVLNode* insertAVLNode(AVLNode* root, AVLNode* node) {
    if(root == NULL) {
        root = createAVLNode(node->key, node->data, NULL, NULL);
    } else if(root->key > node->key) {
        root->lChild = insertAVLNode(root->lChild, node);
        if(getHeightOfAVLTree(root->lChild) - getHeightOfAVLTree(root->rChild) >= 2) {
            if(node->key < root->lChild->key) { // lChild - lChild
                root = rotateAVL_LL(root);
            } else {    // lChild - rChild
                root = rotateAVL_LR(root);
            }
        }
    } else if(root->key < node->key) {
        root->rChild = insertAVLNode(root->rChild, node);
        if(getHeightOfAVLTree(root->rChild) - getHeightOfAVLTree(root->lChild) >= 2) {
            if(node->key > root->rChild->key) { // rChild - rChild
                root = rotateAVL_RR(root);
            } else {    // rChild - lChild
                root = rotateAVL_RL(root);
            }
        }
    } else {    // root->key == node->key
        printf("ERROR: duplicating key value for AVL tree node not allowed. \n");
        exit(EXIT_FAILURE);
    }
    root->height = max_uint64_t(getHeightOfAVLTree(root->lChild),
                                getHeightOfAVLTree(root->rChild)) + 1;
    return root;
}

AVLNode* findAVLNode(AVLNode* root, uint64_t key) {
    AVLNode* node = root;
    while(node != NULL) {
        if(key > node->key) {
            node = node->rChild;
        } else if(key < node->key) {
            node = node->lChild;
        } else {
            return node;
        }
    }
    return node;
}

void clearAVLTree(AVLNode* root) {
    if(root == NULL) {
        return;
    }

    if(root->lChild != NULL) {
        clearAVLTree(root->lChild);
    }
    if(root->rChild != NULL) {
        clearAVLTree(root->rChild);
    }
    free(root);
}


void traversePrintAVLTree(AVLNode* root) {  // breadth first traverse
    uint64_t AVLTreeMaxSize = 2;
    for(uint64_t i = 0; i < root->height; i++) {
        AVLTreeMaxSize = AVLTreeMaxSize * 2;
    }

    uint64_t front = 0, rear = 0;
    AVLNode** queue = (AVLNode**)malloc(sizeof(AVLNode*) * AVLTreeMaxSize);
    if(queue == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    AVLNode* qCell = NULL;

    queue[rear++] = root;
    while(front != rear) {
        qCell = queue[front++];
        printf("key(%"PRIu64") -> data(%"PRIu64")\n", qCell->key, qCell->data);
        if(qCell->lChild != NULL) { queue[rear++] = qCell->lChild; }
        if(qCell->rChild != NULL) { queue[rear++] = qCell->rChild; }
    }

    free(queue);
}


AVLNode** findNodeswithMaxData(AVLNode* root, uint64_t* nodesNum) {
    if(root == NULL){
        *nodesNum = 0;
        return NULL;
    }
    uint64_t AVLTreeMaxSize = 2;
    for(uint64_t i = 0; i < root->height; i++) {
        AVLTreeMaxSize = AVLTreeMaxSize * 2;
    }
    uint64_t front = 0, rear = 0;
    AVLNode** queue = (AVLNode**)malloc(sizeof(AVLNode*) * AVLTreeMaxSize);
    if(queue == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    AVLNode* qCell = NULL;
    uint64_t maxDataValue = 0;

    /** < \note traverse to confirm value of the maximum data */
    front = 0;
    rear = 0;
    queue[rear++] = root;
    while(front != rear) {
        qCell = queue[front++];
        if(qCell->data > maxDataValue) {
            maxDataValue = qCell->data;
        }
        if(qCell->lChild != NULL) { queue[rear++] = qCell->lChild; }
        if(qCell->rChild != NULL) { queue[rear++] = qCell->rChild; }
    }

    /** < \note traverse to find the nodes that have the maximum data */
    AVLNode** nodes = (AVLNode**)malloc(sizeof(AVLNode*) * AVLTreeMaxSize);
    if(nodes == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    *nodesNum = 0;
    front = 0;
    rear = 0;
    queue[rear++] = root;
    while(front != rear) {
        qCell = queue[front++];
        if(qCell->data == maxDataValue) {
            nodes[*nodesNum] = qCell;
            *nodesNum = *nodesNum + 1;
        }
        if(qCell->lChild != NULL) { queue[rear++] = qCell->lChild; }
        if(qCell->rChild != NULL) { queue[rear++] = qCell->rChild; }
    }

    free(queue);
    return nodes;
}












/**
 * Rotate AVL tree node because of inserting lChild to lChild.
 *
 * @param originalRoot original AVL tree root
 * @return new AVL tree root
 */
static AVLNode* rotateAVL_LL(AVLNode* originalRoot) {
    AVLNode* originalLChild = originalRoot->lChild;
    originalRoot->lChild = originalLChild->rChild;
    originalLChild->rChild = originalRoot;

    originalLChild->height = max_uint64_t(getHeightOfAVLTree(originalLChild->lChild),
                                          getHeightOfAVLTree(originalLChild->rChild)) + 1;
    originalRoot->height = max_uint64_t(getHeightOfAVLTree(originalRoot->lChild),
                                        getHeightOfAVLTree(originalRoot->rChild)) + 1;

    return originalLChild;
}

/**
 * Rotate AVL tree node because of inserting rChild to rChild.
 *
 * @param originalRoot original AVL tree root
 * @return new AVL tree root
 */
static AVLNode* rotateAVL_RR(AVLNode* originalRoot) {
    AVLNode* originalRChild = originalRoot->rChild;
    originalRoot->rChild = originalRChild->lChild;
    originalRChild->lChild = originalRoot;

    originalRChild->height = max_uint64_t(getHeightOfAVLTree(originalRChild->lChild),
                                          getHeightOfAVLTree(originalRChild->rChild)) + 1;
    originalRoot->height = max_uint64_t(getHeightOfAVLTree(originalRoot->lChild),
                                        getHeightOfAVLTree(originalRoot->rChild)) + 1;

    return originalRChild;
}

/**
 * Rotate AVL tree node because of inserting rChild to lChild.
 *
 * @param originalRoot original AVL tree root
 * @return new AVL tree root
 */
static AVLNode* rotateAVL_LR(AVLNode* originalRoot) {
    originalRoot->lChild = rotateAVL_RR(originalRoot->lChild);
    return rotateAVL_LL(originalRoot);
}

/**
 * Rotate AVL tree node because of inserting lChild to rChild.
 *
 * @param originalRoot original AVL tree root
 * @return new AVL tree root
 */
static AVLNode* rotateAVL_RL(AVLNode* originalRoot) {
    originalRoot->rChild = rotateAVL_LL(originalRoot->rChild);
    return rotateAVL_RR(originalRoot);
}










