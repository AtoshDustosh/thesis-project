#include "genomeVcf_bPlus.h"

/************************************
 *          Definitions
 ************************************/

struct RecVcf_bplus {
  bcf1_t *data;  // Key is the record's POS
};

struct VcfBPlusNode {
  bool isLeaf;
  VcfBPlusNode *parent;
  // Number of keys filled in this node
  int cnt_key;
  /**
   * @brief  A 1-dimension array of VcfBPlusKey. The size of it depends on
   * whether this node is a leaf node. If it's a leaf node, size(keys) =
   * rank_leaf_node - 1. If it's an inner node, size(keys) = rank_inner_node
   * - 1.
   */
  VcfBPlusKey *keys;

  // Pointers to the node's brothers
  VcfBPlusNode *left;
  VcfBPlusNode *right;

  /**
   * @brief  A 1-dimension array of pointers. The size of it depends on whether
   * this node is a leaf node. If it's a leaf node, size(pointers) =
   * rank_leaf_node - 1. If it's an inner node, size(pointers) = rank_inner_node
   * - 1.
   * These pointers can point to either VcfBPlusNode, or RecVcf_bplus. As
   * for the recognization of the pointers, it depends on the principle of a b
   * plus tree.
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

/************************************
 *         Basic Structures
 ************************************/

/**
 * @brief  Initialize a vcf record object. This method will only initialize
 * memory for keeping vcf record, but will not handle connection to the next
 * record with the same POS. That task should be handled by those methods that
 * call this method.
 * @param  *data: vcf record. This method will create a copy of the record, and
 * thus the original input data in the calling method can be destroyed.
 * @retval A vcf record object. Must be freed using destroy_RecVcf_bplus method.
 */
RecVcf_bplus *init_RecVcf_bplus(bcf1_t *data) {
  RecVcf_bplus *rv = (RecVcf_bplus *)malloc(sizeof(RecVcf_bplus));
  rv->data = bcf_dup(data);
  return rv;
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
 * memory and connection to its parent, but does not initialize connections to
 * the node's brothers. That task should be handled by those methods that call
 * this initialization method.
 * @retval A bpnode object. There is no such method as "destroy_VcfBPlusNode".
 * You should destroy the whole bplus tree or remove one of the records, instead
 * of destroying a single node.
 */
VcfBPlusNode *init_VcfBPlusNode(bool isLeaf, VcfBPlusNode *parent) {
  VcfBPlusNode *bpnode = (VcfBPlusNode *)malloc(sizeof(VcfBPlusNode));
  bpnode->isLeaf = isLeaf;
  bpnode->parent = parent;
  bpnode->cnt_key = 0;
  bpnode->left = NULL;
  bpnode->right =NULL;
  // Assign memories for pointers
  // Use calloc to set all values either NULL or unavailable_keyValue
  if (isLeaf) {
    bpnode->keys =
        (VcfBPlusKey *)calloc(RANK_LEAF_NODE - 1, sizeof(VcfBPlusKey));
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
  VcfBPlusNode *root = init_VcfBPlusNode(true, NULL);
  
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
  gv->hdr = strdup(hdr);
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
        RecVcf_bplus *rv = (RecVcf_bplus*)node->pointers[i];
        genomeVcf_bplus_printRec(gv,rv);
      }
      node = node->right;
    }
    cv = cv->next;
  }
}

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/

// TODO search record in the bptree.
/**
 * @brief  Search for a vcf record with the same POS in the bptree.
 * @param  *rv: vcf record to be searched for
 * @param  *root: root node. Maybe not the root of the tree. But search will
 * start from this node.
 * @retval  **ret_node: bpnode where the vcf record is stored if it exists in
 * the bptree. Otherwise return NULL, which indicates that there is no such
 * record.
 * @retval  *ret_arrayIdx: index of arrays in the bpnode if found the record.
 * Otherwise return -1.
 */
void vcfbplus_searchRec(RecVcf_bplus *rv, VcfBPlusNode *root,
                        VcfBPlusNode **ret_node, int *ret_arrayIdx) {
  if (root->isLeaf) {
    //TODO
  } else {  // If the root node is an inner node
  }
}

// TODO insert records into arrays within nodes in the tree ... and adjust the
// structures of the tree
void vcfbplus_insertRec(RecVcf_bplus *rv, VcfBPlusTree *bptree) {
  // Search the record in bptree
  VcfBPlusNode *root = bptree->root;
  // TODO
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
  // TODO It is probably unnecessary to implement this method.
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

    genomeVcf_bplus_insertRec(gv, init_RecVcf_bplus(tmpRec));

    bcf_destroy1(tmpRec);
    loadedCnt++;
  }

  bcf_hdr_destroy(hdr);
  return gv;
}

void genomeVcf_bplus_writeFile(GenomeVcf_bplus *gv, char *filePath) {
  // TODO
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