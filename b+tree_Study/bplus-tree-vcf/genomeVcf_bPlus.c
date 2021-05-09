#include "genomeVcf_bPlus.h"

/************************************
 *       Definitions: structures
 ************************************/

struct RecVcf_bplus {
  bcf1_t *data;        // Key is the record's POS
  RecVcf_bplus *next;  // Pointer to the next record that have the same POS
};

struct VcfBPlusNode {
  bool isLeaf;
  // Pointer to the parent of this bpnode
  VcfBPlusNode *parent;
  // Indexes of positions in the parent
  int idx_key_parent;
  int idx_pointer_parent;

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
  uint64_t cnt_snp;
  uint64_t cnt_sv;
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
 * @brief  Sub-method for vcfbplus_node_split.
 * @todo   Finish the relationship of methods first before continue implementing this method.
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * // TODO
 * 
 */
void vcfbplus_node_split_inner(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                               int inserted_key_idx, VcfBPlusKey inserted_key,
                               int inserted_pointer_idx,
                               Pointer inserted_pointer);
/**
 * @brief  Sub-method for vcfbplus_node_split.
 * @todo   The splitting process can be optimized ...
 */
void vcfbplus_node_split_leaf(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                              int inserted_idx, VcfBPlusKey inserted_key,
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
                         int inserted_poiner_idx, Pointer inserted_pointer);

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
  rv->data = bcf_dup(data);
  rv->next = next;
  return rv;
}

inline RecVcf_bplus *next_RecVcf_bplus(RecVcf_bplus *rv) { return rv->next; }

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
 * @retval A bpnode object. There is no such method as "destroy_VcfBPlusNode".
 * You should destroy the whole bplus tree or remove one of the records, instead
 * of destroying a single node.
 */
VcfBPlusNode *init_VcfBPlusNode(bool isLeaf) {
  VcfBPlusNode *bpnode = (VcfBPlusNode *)malloc(sizeof(VcfBPlusNode));
  bpnode->isLeaf = isLeaf;
  bpnode->parent = NULL;
  bpnode->idx_key_parent = -1;
  bpnode->idx_pointer_parent = -1;
  bpnode->cnt_key = 0;
  bpnode->left = NULL;
  bpnode->right = NULL;
  /*
   * The following memory allocation is not exactly the same as the official
   * bplus tree principles. Because we have 2 extra pointers to a node's
   * brothers, and thus in a leaf node, there is no need to allocate 1 more
   * memory unit for pointer to brother.
   */
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

VcfBPlusTree *init_VcfBPlusTree() {
  VcfBPlusTree *bptree = (VcfBPlusTree *)calloc(1, sizeof(VcfBPlusTree));
  // Root node is also a leaf node at the beginning.
  VcfBPlusNode *root = init_VcfBPlusNode(true);

  bptree->first = root;
  bptree->last = root;
  bptree->height = 1;
  return bptree;
}

void destroy_VcfBPlusTree(VcfBPlusTree *bptree) {
  // TODO
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
  // TODO initialize tree structures
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

void destroy_GenomeVcf_bplus(GenomeVcf_bplus *gv_bPlus) {
  // TODO Implement this method after finishing all other construction methods
}

void genomeVcf_bplus_traverse(GenomeVcf_bplus *gv) {
  // Print info fields of the Genome
  printf("cnt(snp): %" PRIu64 ", cnt(sv): %" PRIu64 "\n", gv->cnt_snp,
         gv->cnt_sv);
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
      for (int i = 0; i < node->cnt_key + 1; i++) {
        RecVcf_bplus *rv = (RecVcf_bplus *)node->pointers[i];
        // Print all records in the element (linked-list)
        while (rv != NULL) {
          genomeVcf_bplus_printRec(gv, rv);
          rv = next_RecVcf_bplus(rv);
        }
      }
      node = node->right;
    }
    cv = cv->next;
  }
}

/************************************
 * Methods for manipulating Structures
 ************************************/

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

void vcfbplus_node_split_leaf(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                              int inserted_idx, VcfBPlusKey inserted_key,
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
      pointers_all[idx_all] = bpnode->keys[i];
      idx_all++;
    } else if (inserted_key < bpnode->keys[i]) {
      if (if_inserted_copied) {
        keys_all[idx_all] = bpnode->keys[i];
        pointers_all[idx_all] = bpnode->keys[i];
        idx_all++;
      } else {
        if_inserted_copied = true;
        keys_all[idx_all] = inserted_key;
        pointers_all[idx_all] = inserted_pointer;
        idx_all++;
      }
    } else {
      assert(false);
    }
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
  for (int i = cnt_key_left_node + 1; i < RANK_LEAF_NODE + 1; i++) {
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
  VcfBPlusNode *parent = bpnode->parent;
  VcfBPlusKey new_key_parent = new_bpnode->keys[0];
  new_bpnode->parent = parent;
  new_bpnode->idx_key_parent = bpnode->idx_key_parent;
  new_bpnode->idx_pointer_parent = bpnode->idx_pointer_parent;
  // TODO check if the following comment is correct / valid / makes sense ...
  // When the parent is split, these linkers(idx_key_parent and so on)
  // might be changed. But that's not the work of this function.
  vcfbplus_node_insert_inner(bptree, bpnode->parent, bpnode->idx_key_parent + 1,
                             new_key_parent, bpnode->idx_pointer_parent + 1,
                             (Pointer)new_bpnode);
}

void vcfbplus_node_split(VcfBPlusTree *bptree, VcfBPlusNode *bpnode,
                         int inserted_key_idx, VcfBPlusKey inserted_key,
                         int inserted_poiner_idx, Pointer inserted_pointer) {
  // There is no identical key as the inserted key in the split node.
  if (bpnode->isLeaf) {
    // If the node is a leaf node
    vcfbplus_node_split_leaf(bptree, bpnode, inserted_key_idx, inserted_key,
                             inserted_pointer);
  } else {
    // If the node is an inner node
    // TODO
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
        // TODO
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
    // TODO
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
        // TODO
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
    // TODO
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
  // TODO modify after changing the leaf node's elements into linked-list.
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
        l_idx = mid_key + 1;
      } else {
        r_idx = mid_key - 1;
      }
    } else {  // If there is no such record with requested key
      *ret_pointerIdx = l_idx - 1;
      *ret_pointerIdx = *ret_pointerIdx + 1;
      return;
    }
  }
}

void vcfbplus_tree_insertRec(RecVcf_bplus *rv, VcfBPlusTree *bptree) {
  // TODO modify after changing the leaf node's elements into linked-list.
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
}

void genomeVcf_bplus_insertRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv) {
  const char *chromName = rv_chromName(rv, gv);
  // The static variable is used for keeping the last used chromosome.
  // This will accelerate the program by omitting the process of locating
  // chromosome, if requested chrom is the same as the last used one.
  static ChromVcf_bplus *lastUsed_chromVcf;
  if (strcmp(chromName, lastUsed_chromVcf->name) != 0) {
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
    }
  }
  // Insert the record into vcf_bplus tree
  // TODO
}

void genomeVcf_bplus_removeRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv) {
  fprintf(stderr,
          "Error: there is no need to remove a record from the structure - "
          "unimplemented method.\n");
  exit(EXIT_FAILURE);
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

    genomeVcf_bplus_insertRec(gv, init_RecVcf_bplus(tmpRec, NULL));

    bcf_destroy1(tmpRec);
    loadedCnt++;
  }

  bcf_hdr_destroy(hdr);
  return gv;
}

void genomeVcf_bplus_writeFile(GenomeVcf_bplus *gv, char *filePath) {
  // TODO
}

/**********************************
 * Accessing data within structures
 **********************************/

inline bcf1_t *rv_object(RecVcf_bplus *rv) {
  assert(rv->data != NULL);
  return rv->data;
}

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
  // TODO This implementation overlooked special SVs with tags like <INV>.
  assert(
      (alleleIdx < rv->data->n_allele) ||
      (fprintf(
           stderr,
           "Error: array out of boundary for rvDataAlleleCoverLength(...)\n") <
       0));
  int len_ref = strlen(rv->data->d.allele[0]);
  int len_allele = strlen(rv->data->d.allele[alleleIdx]);
  if (len_allele == len_ref)
    return 1;
  else
    return len_ref;
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
  printf("%s\t", bcf_seqname_safe(gv->hdr, rec));
  // pos
  printf("%" PRId64 "\t", rec->pos + 1);
  // id
  printf("%s\t", rec->d.id);
  // ref alt
  printf("%s\t", rec->d.allele[0]);
  if (rec->n_allele == 1) {
    printf(".");
  } else {
    for (int i = 1; i < rec->n_allele; i++) {
      if (i >= 2) {
        printf(",%s(%d)", rec->d.allele[i], bcf_get_variant_type(rec, i));
      } else {
        printf("%s(%d)", rec->d.allele[i], bcf_get_variant_type(rec, i));
      }
    }
  }
  printf("\t");
  // qual
  printf("%f\t", rec->qual);
  // TODO filter, info, format and other fields are ignored
  printf("\n");
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