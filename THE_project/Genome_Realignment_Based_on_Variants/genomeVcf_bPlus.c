#include "genomeVcf_bPlus.h"

/************************************
 *       Definitions: structures
 ************************************/

struct RecVcf_bplus {
  bcf1_t *data;          // Key is the record's POS
  VcfBPlusNode *bpnode;  // Bpnode where this record is kept
  RecVcf_bplus *next;    // Pointer to the next record that have the same POS
};

struct VcfBPlusNode {
  bool isLeaf;
  // Pointer to the parent of this bpnode
  VcfBPlusNode *parent;
  // Indexes of positions in the parent
  // Maybe useful when implementing a "remove" method. But that's not in the
  // plan. So the correctness of these 2 fields is not confirmed.
  // int idx_key_parent;
  // int idx_pointer_parent;

  // Number of keys filled in this node
  int cnt_key;
  /**
   * @brief  A 1-dimension array of VcfBPlusKey. The size of it depends on
   * whether this node is a leaf node. If it's a leaf node, size(keys) =
   * rank_leaf_node. If it's an inner node, size(keys) = rank_inner_node - 1.
   */
  VcfBPlusKey *keys;

  // Pointers to the node's brothers
  VcfBPlusNode *left;
  VcfBPlusNode *right;

  /**
   * @brief  A 1-dimension array of pointers. The size of it depends on whether
   * this node is a leaf node or an inner node. If it's a leaf node,
   * size(pointers) = rank_leaf_node. If it's an inner node, size(pointers) =
   * rank_inner_node. These pointers can point to either VcfBPlusNode, or
   * RecVcf_bplus. As for the recognization of the pointers, it depends on the
   * principle of a b plus tree.
   */
  Pointer *pointers;
};

struct VcfBPlusTree {
  uint16_t height;
  VcfBPlusNode *root;
  // The first leaf node (from left end of data nodes)
  VcfBPlusNode *first;
  // The last leaf node (from right end of data nodes)
  VcfBPlusNode *last;
};

struct ChromVcf_bplus {
  char *name;
  // Number of records within this chrom (bplus tree)
  uint64_t cnt_rec;
  // The tree structure for keeping data
  VcfBPlusTree *tree;
  ChromVcf_bplus *next;
};

struct GenomeVcf_bplus {
  bcf_hdr_t *hdr;
  uint32_t cnt_chrom;
  // This is a linked-list without empty header.
  ChromVcf_bplus *chroms;
};

/************************************
 *       Declarations: functions
 ************************************/

VcfBPlusKey vcfbplus_tree_find_leftmost_key(VcfBPlusNode *bpnode);

VcfBPlusKey vcfbplus_tree_find_rightmost_key(VcfBPlusNode *bpnode);

/**
 * @brief  Split the inner node. Select the middle key of the inner node as the
 * key pushed up into the parent inner node or the root. The middle key will not
 * be kept in the 2 new nodes after splitting.
 * @note   Sub-method for vcfbplus_node_split.
 * @param  *bptree:  bptree.
 * @param  *bpnode:  bpnode to be split.
 * @param  inserted_key:  bplus key to be inserted into "*bpnode"
 * @param  inserted_pointer:  pointer to the new bpnode of lower level. That
 * bpnode can be a leaf node or an inner node.
 */
void vcfbplus_node_split_inner(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                               VcfBPlusKey inserted_key,
                               Pointer inserted_pointer);
/**
 * @brief  Split the leaf node. Select the first key of the right half node as
 * the key pushed up into the parent inner node. The right half node contains 1
 * more key than the left half node.
 * @note   Sub-method for vcfbplus_node_split.
 * @todo   The splitting process can be optimized ...
 */
void vcfbplus_node_split_leaf(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                              VcfBPlusKey inserted_key,
                              Pointer inserted_pointer);

/**
 * @brief  Split a bplus node and set connections and keys in the parent node.
 * @note   Because splitting only happens when a new pointer is inserted and the
 * temporary node is full, we need the temporary node and the inserted poiner
 * together with its inserted position.
 * @param  *bptree: bptree. When splitting nodes, we might need to update the
 * root and height of the bptree.
 * @param  *bpnode: bpnode to be split
 * @param  inserted_key_idx: index of the inserted position for the key
 * @param  inserted_key: key corresponding to the inserted pointer
 * @param  inserted_pointer_idx: index of the inserted position for the pointer
 * @param  inserted_pointer: pointer to the inserted object (record of node)
 */
void vcfbplus_node_split(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                         int inserted_key_idx, VcfBPlusKey inserted_key,
                         int inserted_pointer_idx, Pointer inserted_pointer);

/**
 * @brief  Insert a pointer and its corresponding key into a leaf node. That
 * pointer points to only a vcfbplus record object.
 * @note   This method was created for vcfbplus_node_insert.
 */
void vcfbplus_node_insert_leaf(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                               int inserted_idx, VcfBPlusKey inserted_key,
                               Pointer inserted_pointer);

/**
 * @brief  Insert a pointer and its corresponding key into an inner node. That
 * pointer can point to either a leaf node or another inner node.
 * @note   This method was created for vcfbplus_node_insert at first, but was
 * later found useful for vcfbplus_node_split.
 */
void vcfbplus_node_insert_inner(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                                int inserted_key_idx, VcfBPlusKey inserted_key,
                                int inserted_pointer_idx,
                                Pointer inserted_pointer);

/**
 * @brief  Insert a pointer and its corresponding key into the bpnode. That
 * pointer can point to either a vcf record or another vcf bplus node. If there
 * already exists a record at the inserted position, this method will move the
 * record afterwards together with those records to the right of it in the
 * bpnode, and if the node is full, this method will also split the node.
 * @note
 * @param  *bptree: vcf bplus tree
 * @param  *bpnode: vcf bplus node to be inserted to
 * @param  inserted_key_idx: index of the inserted position for the key
 * @param  inserted_key: key corresponding to the inserted pointer
 * @param  inserted_pointer_idx: index of the inserted position for the pointer
 * @param  inserted_pointer: pointer to the inserted object (record of node)
 */
void vcfbplus_node_insert(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                          int inserted_key_idx, VcfBPlusKey inserted_key,
                          int inserted_pointer_idx, Pointer inserted_pointer);

/**
 * @brief  Locate pointer (path) to the vcf record with requested POS in the
 * bplus node.
 * @param  requested_key: requested key
 * @param  *localRoot: root node. May not be the root of the tree.
 * @retval  *ret_pointerIdx: pointer to the index of pointer (node value) in the
 * node. If this is an inner node, the pointer points to the sub-node which lies
 * on the path to the requested record. If this is a leaf node, the pointer
 * points to the record object if it exists, otherwise the index points to the
 * position where the new record should be inserted.
 */
void vcfbplus_node_locate(VcfBPlusKey requested_key, VcfBPlusNode *localRoot,
                          int *ret_pointerIdx);

/**
 * @brief  Insert a record into the bplus tree. There may already exist records
 * with the same POS (key) as the inserted record's POS. And in such cases, we
 * link the new record with the original list which contains records with the
 * same POS.
 * @param  *rv: inserted record. Note that this method will not create a copy of
 * * the input record object.
 * @param  *bptree: bplus tree
 */
void vcfbplus_tree_insertRec(RecVcf_bplus *rv, VcfBPlusTree *bptree);

/**
 * @brief  Synchronize all records' fields in genomeVcf (update pointers and
 * validate their connections).
 */
void genomeVcf_bplus_synchronize(GenomeVcf_bplus *gv);

/************************************
 *         Basic Structures
 ************************************/

/**
 * @brief  Initialize a vcf record object.
 * @param  *data: vcf record. This method will create a copy of the record, and
 * thus the original input data in the calling method can be destroyed.
 * @param  *next: pointer to the next vcf record that have the same POS. "next"
 * is allowed to be NULL.
 * @retval A vcf record object. Must be freed using destroy_RecVcf_bplus method.
 */
RecVcf_bplus *init_RecVcf_bplus(bcf1_t *data, RecVcf_bplus *next) {
  RecVcf_bplus *rv = (RecVcf_bplus *)malloc(sizeof(RecVcf_bplus));
  rv->bpnode = NULL;
  rv->data = bcf_dup(data);
  bcf_unpack(rv->data, BCF_UN_ALL);
  rv->next = next;
  return rv;
}

RecVcf_bplus *next_RecVcf_bplus(RecVcf_bplus *rv) {
  if (rv->next != NULL) {
    // There exists records with the same POS
    return rv->next;
  } else {
    // Find next record on bpnodes
    VcfBPlusNode *bpnode = rv->bpnode;
    if (bpnode == NULL) {
      fprintf(stderr,
              "Error: need to synchronize genomeVcf before iteration.\n");
      exit(EXIT_FAILURE);
    } else {
      int ret_pointerIdx = 0;
      vcfbplus_node_locate(rv_pos(rv), bpnode, &ret_pointerIdx);
      ret_pointerIdx = ret_pointerIdx + 1;  // Go on to the next record
      if (ret_pointerIdx < bpnode->cnt_key) {
        // Find next record on this bpnode
        return (RecVcf_bplus *)bpnode->pointers[ret_pointerIdx];
      } else {
        // If the next record is not on this bpnode
        bpnode = bpnode->right;
        if (bpnode == NULL) {
          // If this bpnode is the last bpnode
          return NULL;
        } else {
          // According to principles of bplus tree, the next bpnode must be
          // unempty if it exists.
          assert(bpnode->pointers[0] != NULL);
          return (RecVcf_bplus *)bpnode->pointers[0];
        }
      }
    }
  }
}

/**
 * @brief  Destroy a vcf record object. This method will not handle the
 * connection to the next record with the same POS.
 */
void destroy_RecVcf_bplus(RecVcf_bplus *rv) {
  bcf_destroy1(rv_object(rv));
  free(rv);
}

/**
 * @brief  Initialize a bpnode object. This method only intializes the node's
 * memory, and does not initialize connections to
 * the node's brothers or parent. That task should be handled by those methods
 * that call this initialization method.
 * @param  isLeaf:  whether this node is a leaf node.
 * @retval A bpnode object. You should destroy the whole bplus tree or remove
 * one of the records, instead of destroying a single node.
 */
VcfBPlusNode *init_VcfBPlusNode(bool isLeaf) {
  VcfBPlusNode *bpnode = (VcfBPlusNode *)malloc(sizeof(VcfBPlusNode));
  bpnode->isLeaf = isLeaf;
  bpnode->parent = NULL;
  // bpnode->idx_key_parent = -1;
  // bpnode->idx_pointer_parent = -1;
  bpnode->cnt_key = 0;
  bpnode->left = NULL;
  bpnode->right = NULL;
  // Assign memories for pointers
  // Use calloc to set all values either NULL or unavailable_keyValue
  if (isLeaf) {
    bpnode->keys = (VcfBPlusKey *)calloc(RANK_LEAF_NODE, sizeof(VcfBPlusKey));
    bpnode->pointers = (Pointer *)calloc(RANK_LEAF_NODE, sizeof(Pointer));
  } else {
    bpnode->keys =
        (VcfBPlusKey *)calloc(RANK_INNER_NODE - 1, sizeof(VcfBPlusKey));
    bpnode->pointers = (Pointer *)calloc(RANK_INNER_NODE, sizeof(Pointer));
  }
  return bpnode;
}

void destroy_VcfBPlusNode(VcfBPlusNode *bpnode) {
  if (bpnode->isLeaf) {
    // Free every node element of the leaf node
    for (int i = 0; i < bpnode->cnt_key; i++) {
      // Free all elements in the linked-list of a node element
      RecVcf_bplus *rv = bpnode->pointers[i];
      while (rv != NULL) {
        RecVcf_bplus *tmp_rv = rv->next;
        destroy_RecVcf_bplus(rv);
        rv = tmp_rv;
      }
    }
    free(bpnode);
  } else {
    for (int i = 0; i < bpnode->cnt_key + 1; i++) {
      VcfBPlusNode *tmp_bpnode = (VcfBPlusNode *)bpnode->pointers[i];
      destroy_VcfBPlusNode(tmp_bpnode);
    }
    free(bpnode);
  }
}

VcfBPlusTree *init_VcfBPlusTree() {
  VcfBPlusTree *bptree = (VcfBPlusTree *)calloc(1, sizeof(VcfBPlusTree));
  // Root node is also a leaf node at the beginning.
  VcfBPlusNode *root = init_VcfBPlusNode(true);

  bptree->root = root;
  bptree->first = root;
  bptree->last = root;
  bptree->height = 1;
  return bptree;
}

void destroy_VcfBPlusTree(VcfBPlusTree *bptree) {
  destroy_VcfBPlusNode(bptree->root);
  free(bptree);
}

/**
 * @brief  Initialize a ChromVcf_bplus object.
 * @param  *name: name of the chromosome. This method will create a copy of it,
 * so feel free to destroy the original string object.
 * @retval A ChromVcf_bplus object. Must be freed later by calling
 * destroy_ChromVcf_bplus.
 */
ChromVcf_bplus *init_ChromVcf_bplus(const char *name) {
  ChromVcf_bplus *cv = (ChromVcf_bplus *)calloc(1, sizeof(ChromVcf_bplus));
  cv->tree = init_VcfBPlusTree();
  cv->name = strdup(name);
  return cv;
}

void destroy_ChromVcf_bplus(ChromVcf_bplus *cv) {
  free(cv->name);
  destroy_VcfBPlusTree(cv->tree);
  free(cv);
}

GenomeVcf_bplus *init_GenomeVcf_bplus(int rank_inner_node, int rank_leaf_node,
                                      bcf_hdr_t *hdr) {
  RANK_INNER_NODE = rank_inner_node;
  RANK_LEAF_NODE = rank_leaf_node;
  GenomeVcf_bplus *gv = (GenomeVcf_bplus *)calloc(1, sizeof(GenomeVcf_bplus));
  gv->hdr = bcf_hdr_dup(hdr);
  return gv;
}

void destroy_GenomeVcf_bplus(GenomeVcf_bplus *gv) {
  ChromVcf_bplus *cv = gv->chroms;
  while (cv != NULL) {
    ChromVcf_bplus *tmp_cv = cv->next;
    destroy_ChromVcf_bplus(cv);
    cv = tmp_cv;
  }
  free(gv);
}

void genomeVcf_bplus_traverse(GenomeVcf_bplus *gv) {
  // Print info fields of the Genome
  printf("cnt(chrom): %" PRIu32 "\n", gv->cnt_chrom);
  ChromVcf_bplus *cv = gv->chroms;
  while (cv != NULL) {
    // Print info fields of the chromosome
    printf("chrom: %s, cnt(records): %" PRIu64 "\n", cv->name, cv->cnt_rec);
    VcfBPlusNode *node = cv->tree->first;  // start from the first node
    // Print all nodes in the chromosome
    while (node != NULL) {
      assert(node->isLeaf == true);
      // Print all arrays in the node
      for (int i = 0; i < node->cnt_key; i++) {
        RecVcf_bplus *rv = (RecVcf_bplus *)node->pointers[i];
        // Print all records in the element (linked-list)
        while (rv != NULL) {
          genomeVcf_bplus_printRec(gv, rv);
          rv = next_RecVcf_bplus(rv);
        }
      }
      printf("node key cnt: %d, node root: %p\n", node->cnt_key, node->parent);
      node = node->right;
    }
    cv = cv->next;
  }
}

/************************************
 * Methods for manipulating Structures
 ************************************/

void vcfbplus_tree_print(VcfBPlusTree *bptree) {
  static VcfBPlusNode *queue_bpnode[2048];
  int idx_queue_in = 0;
  int idx_queue_out = 0;

  VcfBPlusNode *bpnode = bptree->root;
  queue_bpnode[idx_queue_in++] = bpnode;
  while (idx_queue_out < idx_queue_in) {
    bpnode = queue_bpnode[idx_queue_out++];
    if (bpnode->isLeaf) {
      printf("*");
      for (int i = 0; i < bpnode->cnt_key - 1; i++) {
        printf("%" PRId64 ",", bpnode->keys[i]);
      }
      printf("%" PRId64 "", bpnode->keys[bpnode->cnt_key - 1]);
      printf("*");
      printf("\n");
    } else {
      printf("|");
      for (int i = 0; i < bpnode->cnt_key - 1; i++) {
        printf("%" PRId64 ",", bpnode->keys[i]);
      }
      printf("%" PRId64 "", bpnode->keys[bpnode->cnt_key - 1]);
      printf("|");
      printf("\n");
      // Enqueue sub-bpnodes
      for (int i = 0; i < bpnode->cnt_key + 1; i++) {
        queue_bpnode[idx_queue_in++] = (VcfBPlusNode *)bpnode->pointers[i];
      }
    }
  }
}

VcfBPlusKey vcfbplus_tree_find_leftmost_key(VcfBPlusNode *bpnode) {
  if (bpnode->isLeaf) {
    return bpnode->keys[0];
  } else {
    while (bpnode->pointers[0] != NULL) {
      bpnode = (VcfBPlusNode *)bpnode->pointers[0];
      if (bpnode->isLeaf) {
        return bpnode->keys[0];
      }  // Else, go on to the next left most sub-node
    }
  }
}

VcfBPlusKey vcfbplus_tree_find_rightmost_key(VcfBPlusNode *bpnode) {
  if (bpnode->isLeaf) {
    return bpnode->keys[bpnode->cnt_key - 1];
  } else {
    while (bpnode->pointers[bpnode->cnt_key] != NULL) {
      bpnode = (VcfBPlusNode *)bpnode->pointers[bpnode->cnt_key];
      if (bpnode->isLeaf) {
        return bpnode->keys[bpnode->cnt_key - 1];
      }  // Else, go on to the next right most sub-node
    }
  }
}

void vcfbplus_node_split_inner(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                               VcfBPlusKey inserted_key,
                               Pointer inserted_pointer) {
  assert(bpnode->cnt_key == RANK_INNER_NODE - 1);
  // Integrate the original keys and inserted keys. Do the same to pointers.
  VcfBPlusNode *new_bpnode = init_VcfBPlusNode(false);
  VcfBPlusKey keys_all[RANK_INNER_NODE];
  Pointer pointers_all[RANK_INNER_NODE + 1];
  bool if_inserted_copied = false;
  int idx_all = 0;
  // Although #(pointer) == #(key) + 1, the first pointer will not be changed
  // during integration or splitting.
  pointers_all[0] = bpnode->pointers[0];
  for (int i = 0; i < RANK_INNER_NODE - 1; i++) {
    if (inserted_key > bpnode->keys[i]) {
      keys_all[idx_all] = bpnode->keys[i];
      // The first pointer is ignored
      pointers_all[idx_all + 1] = bpnode->pointers[i + 1];
      idx_all++;
    } else if (inserted_key < bpnode->keys[i]) {
      if (if_inserted_copied) {
        keys_all[idx_all] = bpnode->keys[i];
        pointers_all[idx_all + 1] = bpnode->pointers[i + 1];
        idx_all++;
      } else {
        if_inserted_copied = true;
        keys_all[idx_all] = inserted_key;
        pointers_all[idx_all + 1] = inserted_pointer;
        idx_all++;
        i--;
      }
    } else {
      assert(false);
    }
  }
  if (if_inserted_copied == false) {
    keys_all[RANK_INNER_NODE - 1] = inserted_key;
    pointers_all[RANK_INNER_NODE] = inserted_pointer;
  }
  // Reset the keys and pointers of bpnode
  bpnode->pointers[0] = NULL;
  for (int i = 0; i < RANK_INNER_NODE - 1; i++) {
    bpnode->keys[i] = unavailable_keyValue;
    bpnode->pointers[i + 1] = NULL;
  }
  // Split integrated keys and pointers.
  int cnt_key_left_node = (RANK_INNER_NODE + 1) / 2;
  // Copy to the left node
  bpnode->cnt_key = 0;
  bpnode->pointers[0] = pointers_all[0];
  ((VcfBPlusNode *)pointers_all[0])->parent = bpnode;
  // ((VcfBPlusNode *)pointers_all[0])->idx_key_parent = 0;
  // ((VcfBPlusNode *)pointers_all[0])->idx_pointer_parent = 0;
  for (int i = 0; i < cnt_key_left_node - 1; i++) {
    ((VcfBPlusNode *)(pointers_all[i + 1]))->parent = bpnode;
    bpnode->keys[i] = keys_all[i];
    // ((VcfBPlusNode *)(pointers_all[i + 1]))->idx_key_parent = i;
    // bpnode->pointers[0] is ignored
    bpnode->pointers[i + 1] = pointers_all[i + 1];
    // ((VcfBPlusNode *)(pointers_all[i + 1]))->idx_pointer_parent = i + 1;
    bpnode->cnt_key = bpnode->cnt_key + 1;
  }
  // Copy to the right node
  int idx_new_bpnode = 0;
  new_bpnode->cnt_key = 0;
  new_bpnode->pointers[0] = pointers_all[cnt_key_left_node];
  ((VcfBPlusNode *)pointers_all[cnt_key_left_node])->parent = new_bpnode;
  // ((VcfBPlusNode *)pointers_all[cnt_key_left_node])->idx_key_parent = 0;
  // ((VcfBPlusNode *)pointers_all[cnt_key_left_node])->idx_pointer_parent = 0;
  for (int i = cnt_key_left_node + 1; i < RANK_INNER_NODE + 1; i++) {
    ((VcfBPlusNode *)(pointers_all[i]))->parent = new_bpnode;
    new_bpnode->keys[idx_new_bpnode] = keys_all[i - 1];
    // ((VcfBPlusNode *)(pointers_all[i]))->idx_key_parent = idx_new_bpnode;
    new_bpnode->pointers[idx_new_bpnode + 1] = pointers_all[i];
    // ((VcfBPlusNode *)(pointers_all[i]))->idx_pointer_parent =
    idx_new_bpnode + 1;
    new_bpnode->cnt_key = new_bpnode->cnt_key + 1;
    idx_new_bpnode++;
  }
  // Push up the middle key
  VcfBPlusKey key_mid = keys_all[cnt_key_left_node - 1];
  if (bpnode->parent != NULL) {
    // If this inner node has a parent node (an inner node)
    VcfBPlusNode *parent = bpnode->parent;
    int idx_pointer_parent = 0;
    int idx_key_parent = 0;
    vcfbplus_node_locate(key_mid, parent, &idx_key_parent);
    idx_pointer_parent = idx_key_parent + 1;
    new_bpnode->parent = parent;
    // new_bpnode->idx_pointer_parent = idx_pointer_parent;
    // new_bpnode->idx_key_parent = idx_key_parent;
    vcfbplus_node_insert_inner(bptree, parent, idx_key_parent, key_mid,
                               idx_pointer_parent, (Pointer)new_bpnode);
  } else {
    // If this inner node is the root.
    VcfBPlusNode *parent = init_VcfBPlusNode(false);
    parent->keys[0] = key_mid;
    parent->pointers[0] = bpnode;
    parent->pointers[1] = new_bpnode;
    parent->cnt_key = 1;
    bpnode->parent = parent;
    // bpnode->idx_key_parent = 0;
    // bpnode->idx_pointer_parent = 0;
    new_bpnode->parent = parent;
    // new_bpnode->idx_key_parent = 0;
    // new_bpnode->idx_pointer_parent = 1;
    bptree->root = parent;
    bptree->height = bptree->height + 1;
  }
}

void vcfbplus_node_split_leaf(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                              VcfBPlusKey inserted_key,
                              Pointer inserted_pointer) {
  assert(bpnode->cnt_key == RANK_LEAF_NODE);
  // Integrate the original keys and inserted key. Do the same to pointers.
  VcfBPlusNode *new_bpnode = init_VcfBPlusNode(true);
  VcfBPlusKey keys_all[RANK_LEAF_NODE + 1];
  Pointer pointers_all[RANK_LEAF_NODE + 1];
  bool if_inserted_copied = false;
  int idx_all = 0;
  for (int i = 0; i < RANK_LEAF_NODE; i++) {
    if (inserted_key > bpnode->keys[i]) {
      keys_all[idx_all] = bpnode->keys[i];
      pointers_all[idx_all] = bpnode->pointers[i];
      idx_all++;
    } else if (inserted_key < bpnode->keys[i]) {
      if (if_inserted_copied) {
        keys_all[idx_all] = bpnode->keys[i];
        pointers_all[idx_all] = bpnode->pointers[i];
        idx_all++;
      } else {
        if_inserted_copied = true;
        keys_all[idx_all] = inserted_key;
        pointers_all[idx_all] = inserted_pointer;
        idx_all++;
        i--;
      }
    } else {
      assert(false);
    }
  }
  if (if_inserted_copied == false) {
    keys_all[RANK_LEAF_NODE] = inserted_key;
    pointers_all[RANK_LEAF_NODE] = inserted_pointer;
  }
  // Reset keys and pointers of bpnode
  for (int i = 0; i < RANK_LEAF_NODE; i++) {
    bpnode->keys[i] = unavailable_keyValue;
    bpnode->pointers[i] = NULL;
  }
  // Split integrated keys and pointers into 2 nodes.
  int cnt_key_left_node = (RANK_LEAF_NODE + 1) / 2;
  bpnode->cnt_key = 0;
  for (int i = 0; i < cnt_key_left_node; i++) {
    bpnode->keys[i] = keys_all[i];
    bpnode->pointers[i] = pointers_all[i];
    bpnode->cnt_key = bpnode->cnt_key + 1;
  }
  int idx_key_new_node = 0;
  for (int i = cnt_key_left_node; i < RANK_LEAF_NODE + 1; i++) {
    new_bpnode->keys[idx_key_new_node] = keys_all[i];
    new_bpnode->pointers[idx_key_new_node] = pointers_all[i];
    new_bpnode->cnt_key = new_bpnode->cnt_key + 1;
    idx_key_new_node++;
  }
  // Assign connections for the new bpnode.
  // This might cause further splitting of the parent node.
  new_bpnode->right = bpnode->right;
  new_bpnode->left = bpnode;
  bpnode->right = new_bpnode;
  // Update parent keys and pointers
  if (bpnode->parent != NULL) {
    // This leaf node has a parent node (an inner node)
    VcfBPlusNode *parent = bpnode->parent;
    VcfBPlusKey new_key_parent = new_bpnode->keys[0];
    // new_bpnode->parent = parent;
    // new_bpnode->idx_key_parent = bpnode->idx_key_parent + 1;
    // new_bpnode->idx_pointer_parent = bpnode->idx_pointer_parent + 1;
    int idx_key_parent = 0;
    int idx_pointer_parent = 0;
    vcfbplus_node_locate(new_key_parent, parent, &idx_key_parent);
    idx_pointer_parent = idx_key_parent + 1;
    new_bpnode->parent = parent;
    // When the parent is split, these linkers(idx_key_parent and so on) for
    // "new_bpnode" might be changed. But that's not the work of this function.
    vcfbplus_node_insert_inner(bptree, bpnode->parent, idx_key_parent,
                               new_key_parent, idx_pointer_parent,
                               (Pointer)new_bpnode);
    if (new_bpnode->right == NULL) {
      bptree->last = new_bpnode;
    }
  } else {
    // Codes within this brackets run only once during the creation of the tree.
    // This leaf node is the root node. Split the root node.
    VcfBPlusNode *parent = init_VcfBPlusNode(false);
    parent->keys[0] = new_bpnode->keys[0];
    parent->pointers[0] = bpnode;
    parent->pointers[1] = new_bpnode;
    parent->cnt_key = 1;
    bpnode->parent = parent;
    // bpnode->idx_key_parent = 0;
    // bpnode->idx_pointer_parent = 0;
    new_bpnode->parent = parent;
    // new_bpnode->idx_key_parent = 0;
    // new_bpnode->idx_pointer_parent = 1;
    bptree->root = parent;
    bptree->first = bpnode;
    bptree->last = new_bpnode;
    bptree->height = bptree->height + 1;
  }
}

void vcfbplus_node_split(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                         int inserted_key_idx, VcfBPlusKey inserted_key,
                         int inserted_pointer_idx, Pointer inserted_pointer) {
  // There is no identical key as the inserted key in the split node.
  if (bpnode->isLeaf) {
    // If the node is a leaf node
    vcfbplus_node_split_leaf(bptree, bpnode, inserted_key, inserted_pointer);
  } else {
    // If the node is an inner node
    vcfbplus_node_split_inner(bptree, bpnode, inserted_key, inserted_pointer);
  }
  return;
}

void vcfbplus_node_insert_leaf(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                               int inserted_idx, VcfBPlusKey inserted_key,
                               Pointer inserted_pointer) {
  if (inserted_idx < RANK_LEAF_NODE) {
    // Insert the new key into one of the positions in the node
    VcfBPlusKey original_key = bpnode->keys[inserted_idx];
    if (original_key == unavailable_keyValue) {
      // Empty position. Directly insert.
      assert(bpnode->pointers[inserted_idx] == NULL);
      bpnode->keys[inserted_idx] = inserted_key;
      bpnode->pointers[inserted_idx] = inserted_pointer;
      bpnode->cnt_key = bpnode->cnt_key + 1;
    } else if (original_key == inserted_key) {
      // Inserted position is occupied by records with the same key (POS)
      RecVcf_bplus *rv = (RecVcf_bplus *)bpnode->pointers[inserted_idx];
      while (rv->next != NULL) {
        rv = rv->next;
      }
      rv->next = (RecVcf_bplus *)inserted_pointer;
    } else {
      // Inserted position is occupied by another record.
      if (bpnode->cnt_key == RANK_LEAF_NODE) {
        // This leaf node is full. Need to split the node
        vcfbplus_node_split_leaf(bptree, bpnode, inserted_key,
                                 inserted_pointer);
      } else if (bpnode->cnt_key < RANK_LEAF_NODE) {
        // This leaf node is not full. Insert and move records afterwards.
        VcfBPlusKey tmp_inserted_key = inserted_key;
        Pointer tmp_inserted_pointer = inserted_pointer;
        for (int i = inserted_idx; i < bpnode->cnt_key; i++) {
          VcfBPlusKey tmp_key = bpnode->keys[i];
          Pointer tmp_pointer = bpnode->pointers[i];
          bpnode->keys[i] = tmp_inserted_key;
          bpnode->pointers[i] = tmp_inserted_pointer;
          tmp_inserted_key = tmp_key;
          tmp_inserted_pointer = tmp_pointer;
        }
        bpnode->keys[bpnode->cnt_key] = tmp_inserted_key;
        bpnode->pointers[bpnode->cnt_key] = tmp_inserted_pointer;
        bpnode->cnt_key = bpnode->cnt_key + 1;
      } else {
        assert(false);
      }
    }
  } else {
    // Need to split the node
    // The inserted position is at the end of the node, but out of boundary.
    vcfbplus_node_split_leaf(bptree, bpnode, inserted_key, inserted_pointer);
  }
}

void vcfbplus_node_insert_inner(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                                int inserted_key_idx, VcfBPlusKey inserted_key,
                                int inserted_pointer_idx,
                                Pointer inserted_pointer) {
  if (inserted_key_idx < RANK_INNER_NODE - 1) {
    assert(inserted_pointer_idx < RANK_INNER_NODE);
    // Insert the new key into one of the positions for keys in the node
    if (bpnode->keys[inserted_key_idx] == unavailable_keyValue) {
      // Empty position. Directly insert.
      bpnode->keys[inserted_key_idx] = inserted_key;
      bpnode->pointers[inserted_pointer_idx] = inserted_pointer;
      bpnode->cnt_key = bpnode->cnt_key + 1;
    } else {
      // Inserted position is occupied by another pointer.
      if (bpnode->cnt_key == RANK_INNER_NODE - 1) {
        // This inner node is full. Need to split the node.
        vcfbplus_node_split_inner(bptree, bpnode, inserted_key,
                                  inserted_pointer);
      } else {
        // This inner node is not full. Insert and move pointers afterwards.
        VcfBPlusKey tmp_inserted_key = inserted_key;
        for (int i = inserted_key_idx; i < bpnode->cnt_key; i++) {
          VcfBPlusKey tmp_key = bpnode->keys[i];
          bpnode->keys[i] = tmp_inserted_key;
          tmp_inserted_key = tmp_key;
        }
        bpnode->keys[bpnode->cnt_key] = tmp_inserted_key;
        // For inner nodes, it always satisfies cnt_pointer == cnt_key + 1.
        Pointer tmp_inserted_pointer = inserted_pointer;
        for (int i = inserted_pointer_idx; i < bpnode->cnt_key + 1; i++) {
          Pointer tmp_pointer = bpnode->pointers[i];
          bpnode->pointers[i] = tmp_inserted_pointer;
          tmp_inserted_pointer = tmp_pointer;
        }
        bpnode->pointers[bpnode->cnt_key + 1] = tmp_inserted_pointer;
        bpnode->cnt_key = bpnode->cnt_key + 1;
      }
    }
  } else {
    // Need to split the node
    // The inserted position is at the end of the node, but out of boundary.
    vcfbplus_node_split_inner(bptree, bpnode, inserted_key, inserted_pointer);
  }
}

void vcfbplus_node_insert(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                          int inserted_key_idx, VcfBPlusKey inserted_key,
                          int inserted_pointer_idx, Pointer inserted_pointer) {
  if (bpnode->isLeaf) {
    // If the node is a leaf node
    vcfbplus_node_insert_leaf(bptree, bpnode, inserted_key_idx, inserted_key,
                              inserted_pointer);
  } else {
    // If the node is an inner node
    vcfbplus_node_insert_inner(bptree, bpnode, inserted_key_idx, inserted_key,
                               inserted_pointer_idx, inserted_pointer);
  }
}

void vcfbplus_node_locate(VcfBPlusKey requested_key, VcfBPlusNode *localRoot,
                          int *ret_pointerIdx) {
  *ret_pointerIdx = -1;
  int l_idx = 0;
  int r_idx = localRoot->cnt_key - 1;
  while (1) {
    // Binary search the position of the requested key
    if (l_idx <= r_idx) {
      int mid_idx = (l_idx + r_idx) / 2;
      VcfBPlusKey mid_key = localRoot->keys[mid_idx];
      if (requested_key == mid_key) {
        *ret_pointerIdx = mid_idx;
        if (localRoot->isLeaf == false) *ret_pointerIdx = *ret_pointerIdx + 1;
        return;
      } else if (mid_key < requested_key) {
        l_idx = mid_idx + 1;
      } else {
        r_idx = mid_idx - 1;
      }
    } else {  // If there is no such record with requested key
      *ret_pointerIdx = l_idx - 1;
      *ret_pointerIdx = *ret_pointerIdx + 1;
      return;
    }
  }
}

void vcfbplus_tree_insertRec(RecVcf_bplus *rv, VcfBPlusTree *bptree) {
  // Locate the position where the new record should be inserted
  VcfBPlusNode *bpnode = bptree->root;
  VcfBPlusKey inserted_key = rv_pos(rv);

  // Locate the leaf node where the record should be inserted
  int ret_pointerIdx = -1;
  while (bpnode->isLeaf == false) {
    vcfbplus_node_locate(inserted_key, bpnode, &ret_pointerIdx);
    bpnode = (VcfBPlusNode *)(bpnode->pointers[ret_pointerIdx]);
  }

  // Locate the position in the leaf node where the record should be inserted
  assert(bpnode->isLeaf == true);
  vcfbplus_node_locate(inserted_key, bpnode, &ret_pointerIdx);

  // Insert
  vcfbplus_node_insert_leaf(bptree, bpnode, ret_pointerIdx, inserted_key,
                            (Pointer)rv);
}

void genomeVcf_bplus_insertRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv) {
  const char *chromName = rv_chromName(rv, gv);
  // The static variable is used for keeping the last used chromosome.
  // This will accelerate the program by omitting the process of locating
  // chromosome, if requested chrom is the same as the last used one.
  static ChromVcf_bplus *lastUsed_chromVcf;
  if (lastUsed_chromVcf == NULL ||
      strcmp(chromName, lastUsed_chromVcf->name) != 0) {
    // Locate the chromosome (check from the 1st chrom in the linked-list)
    lastUsed_chromVcf = gv->chroms;
    while (lastUsed_chromVcf != NULL) {
      if (strcmp(lastUsed_chromVcf->name, chromName) == 0) {
        break;
      } else {
        lastUsed_chromVcf = lastUsed_chromVcf->next;
      }
    }
    // Possibly there exists no such chromosome
    if (lastUsed_chromVcf == NULL) {
      lastUsed_chromVcf = init_ChromVcf_bplus(chromName);
      ChromVcf_bplus *tmp_chromVcf = gv->chroms;
      if (tmp_chromVcf == NULL) {
        gv->chroms = lastUsed_chromVcf;
      } else {
        while (tmp_chromVcf != NULL) {
          if (tmp_chromVcf->next == NULL) {
            tmp_chromVcf->next = lastUsed_chromVcf;
            break;
          } else {
            tmp_chromVcf = tmp_chromVcf->next;
          }
        }
      }
      gv->cnt_chrom++;
    }
  }
  // Insert the record into vcf_bplus tree
  vcfbplus_tree_insertRec(rv, lastUsed_chromVcf->tree);
  lastUsed_chromVcf->cnt_rec++;
}

void genomeVcf_bplus_removeRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv) {
  fprintf(stderr,
          "Error: there is no need to remove a record from the structure - "
          "unimplemented method.\n");
  exit(EXIT_FAILURE);
}

RecVcf_bplus *genomeVcf_bplus_getRecAfterPos(GenomeVcf_bplus *gv,
                                             const char *chromName,
                                             int64_t pos) {
  assert(pos >= 1);
  // Locate the chromosome of vcf records
  // The following codes used "static" to accelerate the process. But in a
  // multi-thread program, I'm afraid that this may result in bugs. Because I
  // didn't add a thread-lock here. So we used the simple impelmention.
  /*
  static ChromVcf_bplus *lastUsed_chromVcf;
  if (lastUsed_chromVcf == NULL ||
      strcmp(lastUsed_chromVcf->name, chromName) != 0) {
    lastUsed_chromVcf = gv->chroms;
    while (lastUsed_chromVcf != NULL) {
      if (strcmp(lastUsed_chromVcf->name, chromName) == 0) {
        break;
      } else {
        lastUsed_chromVcf = lastUsed_chromVcf->next;
      }
    }
  }
  if (lastUsed_chromVcf == NULL) {
    // Possibly there is no requested chrom in genomeVcf
    return NULL;
  }
  // Locate and find record(s)
  VcfBPlusTree *bptree = lastUsed_chromVcf->tree;
  */

  ChromVcf_bplus *cf = gv->chroms;
  while (cf != NULL) {
    if (strcmp(cf->name, chromName) == 0) {
      break;
    } else {
      cf = cf->next;
    }
  }
  if (cf == NULL) {
    return NULL;
  }
  VcfBPlusTree *bptree = cf->tree;

  VcfBPlusNode *bpnode = bptree->root;

  int ret_pointerIdx = -1;
  while (bpnode->isLeaf == false) {
    vcfbplus_node_locate(pos, bpnode, &ret_pointerIdx);
    bpnode = (VcfBPlusNode *)(bpnode->pointers[ret_pointerIdx]);
  }

  vcfbplus_node_locate(pos, bpnode, &ret_pointerIdx);

  if (ret_pointerIdx >= bpnode->cnt_key) {
    // Relocate in the next bpnode
    bpnode = bpnode->right;
    if (bpnode == NULL) {
      return NULL;
    } else {
      assert(bpnode->pointers[0] != NULL);
    }
    vcfbplus_node_locate(pos, bpnode, &ret_pointerIdx);
  }
  return (RecVcf_bplus *)bpnode->pointers[ret_pointerIdx];
}

void genomeVcf_bplus_synchronize(GenomeVcf_bplus *gv) {
  ChromVcf_bplus *cv = gv->chroms;
  while (cv != NULL) {
    VcfBPlusNode *node = cv->tree->first;  // start from the first node
    // Synchronize all nodes in the chromosome
    while (node != NULL) {
      assert(node->isLeaf == true);
      // Synchronize all arrays in the node
      for (int i = 0; i < node->cnt_key; i++) {
        RecVcf_bplus *rv = (RecVcf_bplus *)node->pointers[i];
        // Synchronize all records in the element (linked-list)
        while (rv != NULL) {
          rv->bpnode = node;
          rv = rv->next;
        }
      }
      node = node->right;
    }
    cv = cv->next;
  }
}

GenomeVcf_bplus *genomeVcf_bplus_loadFile(char *filePath, int rank_inner_node,
                                          int rank_leaf_node) {
  htsFile *fp = hts_open(filePath, "r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init1();

  if (hdr == NULL) {
    fprintf(stderr, "Error: failed creating bcf header struct.\n");
    exit(EXIT_FAILURE);
  }
  if (rec == NULL) {
    fprintf(stderr, "Error: memory not enough for new bcf1_t object.\n");
    exit(EXIT_FAILURE);
  }

  GenomeVcf_bplus *gv =
      init_GenomeVcf_bplus(rank_inner_node, rank_leaf_node, hdr);
  uint32_t loadedCnt = 0;
  while (bcf_read1(fp, hdr, rec) >= 0) {
    // The "bcf_unpack" method must be called for every new bcf1_t object
    bcf_unpack(rec, BCF_UN_ALL);
    bcf1_t *tmpRec = bcf_dup(rec);
    bcf_unpack(tmpRec, BCF_UN_ALL);

    RecVcf_bplus *rv = init_RecVcf_bplus(tmpRec, NULL);
    bcf_destroy1(tmpRec);

    // Ignore those variants with tags or with empty ALT fields
    if (rv_alleleCnt(rv) == 1) {
      continue;
    } else {
      const char *vcf_ALT = rv_allele(rv, 1);
      if (vcf_ALT[0] == '<') {
        continue;
      }
    }

    genomeVcf_bplus_insertRec(gv, rv);

    // Debug lines: used for checking the correctness of the bplus tree
    // vcfbplus_tree_print(gv->chroms->tree);

    // genomeVcf_bplus_printRec(gv, rv);
    loadedCnt++;
  }

  genomeVcf_bplus_synchronize(gv);

  bcf_hdr_destroy(hdr);
  return gv;
}

void genomeVcf_bplus_writeFile(GenomeVcf_bplus *gv, char *filePath) {
  fprintf(stderr, "Error: this method is unnecessary to implement.\n");
  exit(EXIT_FAILURE);
}

/**********************************
 * Accessing data within structures
 **********************************/

inline int gv_cnt_rec(GenomeVcf_bplus *gv) {
  int cnt_rec = 0;
  ChromVcf_bplus *cv = gv->chroms;
  while (cv != NULL) {
    cnt_rec = cnt_rec + cv->cnt_rec;
    cv = cv->next;
  }
  return cnt_rec;
}

inline bcf1_t *rv_object(RecVcf_bplus *rv) {
  assert(rv->data != NULL);
  return rv->data;
}

const char *rv_ID(RecVcf_bplus *rv) { return rv->data->d.id; }

inline int64_t rv_pos(RecVcf_bplus *rv) {
  assert(rv->data != NULL);
  return rv->data->pos + 1;
}

inline const char *rv_chromName(RecVcf_bplus *rv, GenomeVcf_bplus *gv) {
  assert(gv->hdr != NULL);
  assert(rv->data != NULL);
  return bcf_seqname_safe(gv->hdr, rv->data);
}

inline int rv_alleleCnt(RecVcf_bplus *rv) {
  assert(rv->data != NULL);
  return rv->data->n_allele;
}

inline int rv_alleleCoverLength(RecVcf_bplus *rv, int alleleIdx) {
  assert(rv->data != NULL);
  assert(
      (alleleIdx < rv->data->n_allele) ||
      (fprintf(stderr,
               "Error: array out of boundary for rv_alleleCoverLength(...)\n") <
       0));
  // Variants using TAGs ignored
  assert((rv_allele(rv, alleleIdx)[0] != '<') ||
         (fprintf(stderr,
                  "Error: variants with tags should be ignored when "
                  "loading.\n") < 0));
  int len_ref = strlen(rv_allele(rv, 0));
  int len_allele = strlen(rv_allele(rv, alleleIdx));
  if (len_ref == 1) {  // Insertion
    return 1;
  } else if (len_allele == 1) {  // Deletion
    return len_ref;
  } else {  // Other cases (MNP, etc.)
    return len_ref;
  }
}

inline const char *rv_allele(RecVcf_bplus *rv, int alleleIdx) {
  assert(rv->data != NULL);
  assert(
      (alleleIdx < rv->data->n_allele) ||
      (fprintf(
           stderr,
           "Error: array out of boundary for rvDataAlleleCoverLength(...)\n") <
       0));
  return rv->data->d.allele[alleleIdx];
}

void genomeVcf_bplus_printRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv) {
  assert(gv->hdr != NULL);
  assert(rv->data != NULL);
  bcf1_t *rec = rv_object(rv);
  // chrom
  fprintf(stderr, "%s\t", bcf_seqname_safe(gv->hdr, rec));
  // pos
  fprintf(stderr, "%" PRId64 "\t", rec->pos + 1);
  // id
  fprintf(stderr, "%s\t", rec->d.id);
  // ref alt
  fprintf(stderr, "%s\t", rec->d.allele[0]);
  if (rec->n_allele == 1) {
    fprintf(stderr, ".");
  } else {
    for (int i = 1; i < rec->n_allele; i++) {
      if (i >= 2) {
        fprintf(stderr, ",%s(%d)", rec->d.allele[i],
                bcf_get_variant_type(rec, i));
      } else {
        fprintf(stderr, "%s(%d)", rec->d.allele[i],
                bcf_get_variant_type(rec, i));
      }
    }
  }
  fprintf(stderr, "\t");
  // qual
  fprintf(stderr, "%f\t", rec->qual);
  // TODO filter, info, format and other fields are ignored
  fprintf(stderr, "\n");
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/************************* Debug Methods ************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

/**********************************
 * Debugging Methods for GenomeVcf
 **********************************/

static int _test_BasicFunctions() { return 1; }

void _testSet_genomeVcf_bplus() { assert(_test_BasicFunctions()); }