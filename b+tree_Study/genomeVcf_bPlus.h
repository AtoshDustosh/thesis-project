#ifndef GENOMEVCF_BPLUS_H_INCLUDED
#define GENOMEVCF_BPLUS_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"

/**
 * @brief  Rank of bplus tree
 * @note   A global variable with "static const" feels like "public static
 * final" in java. But on second thought, I may need to input the rank of leaf
 * node and inner node, so let's just make them "static".
 */
static int RANK_INNER_NODE = 6;
static int RANK_LEAF_NODE = 4;

// Type of position in vcf records
typedef int64_t KeyType;

/**
 * @brief  Vcf record.
 */
typedef struct RecVcf_bplus RecVcf_bplus;

typedef struct VcfBPlusNode VcfBPlusNode;

/**
 * @brief  Use bplus tree to index the vcf records under this chromosome.
 */
typedef struct ChromVcf_bplus ChromVcf_bplus;

/**
 * @brief  Use linked-list to index all chroms.
 */
typedef struct GenomeVcf_bplus GenomeVcf_bplus;

/**********************************
 * Accessing data within structures
 **********************************/

/**
 * @brief  Get the pointer to the original data object within rv.
 */
extern bcf1_t *rv_object(RecVcf_bplus *rv);

/**
 * @brief  Get the 1-based position of the variant.
 */
extern int64_t *rv_pos(RecVcf_bplus *rv);

/**
 * @brief  Return a copy of chrom name got from RecVcf object.
 * @retval chrom name. The returned string must be freed manually later.
 */
extern char *rv_contigName(RecVcf_bplus *rv, GenomeVcf_bplus *gv);

/**
 * @brief  Return the count of alleles the RecVcf object contains.
 * @retval   Ref is also included as an allele, thus returned value >= 1
 */
extern int rv_alleleCnt(RecVcf_bplus *rv);

/**
 * * @brief  Return the covered length of this allele.
 * There are several cases (REF ALT):
 * 1. (A  ACCC) - INS - return 1
 * 2. (A  C) - SNP - return 1
 * 3. (ACGTA  A) - DEL - return 5
 * 4. (CCC  C,CCCCCCCC) - DEL,INS - return 3,3
 */
extern int rv_alleleCoverLength(RecVcf_bplus *rv, int alleleIdx);

/**
 * @brief  Get the string of idx-th allele of the RecVcf object.
 * @retval pointer to string of allele; NULL if the idx is invalid. Do not free
 * the string. It's not a copy of the original data.
 */
extern const char *rv_allele(RecVcf_bplus *rv, int alleleIdx);

void genomeVcf_bplus_printRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv);

/************************************
 *           Basic Structures
 ************************************/

GenomeVcf_bplus *init_GenomeVcf_bplus(int rank_inner_node, int rank_leaf_node);

void destroy_GenomeVcf_bplus(GenomeVcf_bplus *gv);

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/
void genomeVcf_bplus_insertRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv);

void genomeVcf_bplus_removeRec(GenomeVcf_bplus *gv, RecVcf_bplus *rv);

void genomeVcf_bplus_loadFile(GenomeVcf_bplus *gv, char *filePath);

void genomeVcf_bplus_writeFile(GenomeVcf_bplus *gv, char *filePath);

/**********************************
 * Debugging Methods for GenomeVcf
 **********************************/

void _testSet_genomeVcf_bplus();

#endif