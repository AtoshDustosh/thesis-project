#include "genomeVcf_bPlus.h"

/******************
 *   Definitions
 ******************/

struct RecVcf_bplus {
  bcf1_t *rec;
};

struct VcfBPlusNode {
  bool isLeaf;
  struct _define_VcfBPlusNode **parent;
  // Number of children (inner node) or elements (leaf node)
  int pointerCnt;
  /*
   * If this node is an inner node.
   */
  // Dynamically allocated 1-dimension array
  struct _define_VcfBPlusNode **children;
  /*
   * If this node is a leaf node.
   */
  // Dynamically allocated 1-dimension array
  struct _define_VcfBPlusNode **elements;
  // brother node
  struct _define_VcfBPlusNode *next;
};

struct ChromVcf_bplus {
  char *name;
  // Number of records within this chrom (bplus tree)
  uint64_t recCnt;
  // Root of the bplus tree
  VcfBPlusNode *bplus_root;
  // The left most node (contains the first record)
  VcfBPlusNode *bplus_start;
  struct _define_ChromVcf_bplus *next;
};

struct GenomeVcf_bplus {
  bcf_hdr_t *hdr;
  uint64_t snpCnt;
  uint64_t svCnt;
  uint32_t chromCnt;
  // This is a linked-list without empty header.
  ChromVcf_bplus *chroms;
};

/**********************************
 * Accessing data within structures
 **********************************/

inline bcf1_t *rv_object(RecVcf_bplus *rv) { return rv->rec; }

inline int64_t *rv_pos(RecVcf_bplus *rv) { return rv->rec->pos + 1; }

inline char *rv_contigName(RecVcf_bplus *rv, GenomeVcf_bplus *gv) {
  return strdup(bcf_seqname_safe(gv->hdr, rv->rec));
}

inline int rv_alleleCnt(RecVcf_bplus *rv) { return rv->rec->n_allele; }

inline int rv_alleleCoverLength(RecVcf_bplus *rv, int alleleIdx) {
  // TODO This implementation overlooked special SVs with tags like <INV>.
  assert(
      (alleleIdx < rv->rec->n_allele) ||
      (fprintf(
           stderr,
           "Error: array out of boundary for rvDataAlleleCoverLength(...)\n") <
       0));
  int len_ref = strlen(rv->rec->d.allele[0]);
  int len_allele = strlen(rv->rec->d.allele[alleleIdx]);
  if (len_allele == len_ref)
    return 1;
  else
    return len_ref;
}

inline const char *rv_allele(RecVcf_bplus *rv, int alleleIdx) {
  assert(
      (alleleIdx < rv->rec->n_allele) ||
      (fprintf(
           stderr,
           "Error: array out of boundary for rvDataAlleleCoverLength(...)\n") <
       0));
  return rv->rec->d.allele[alleleIdx];
}

void genomeVcf_bplus_printRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv) {
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

/******************
 * Basic Structures
 ******************/

GenomeVcf_bplus *init_GenomeVcf_bplus(int rank_inner_node, int rank_leaf_node) {
  RANK_INNER_NODE = rank_inner_node;
  RANK_LEAF_NODE = rank_leaf_node;
  return (GenomeVcf_bplus *)calloc(1, sizeof(GenomeVcf_bplus));
}

void destroy_GenomeVcf_bplus(GenomeVcf_bplus *gv_bPlus) {
  // TODO Implement this method after finishing all other construction methods
}

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/
void genomeVcf_bplus_insertRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv){
  static ChromVcf_bplus *lastUsed_chromVcf;
}

void genomeVcf_bplus_removeRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv);

void genomeVcf_bplus_loadFile(GenomeVcf_bplus *gv, char *filePath);

void genomeVcf_bplus_writeFile(GenomeVcf_bplus *gv, char *filePath);

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