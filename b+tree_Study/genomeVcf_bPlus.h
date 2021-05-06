#ifndef GENOMEVCF_BPLUS_H_INCLUDED
#define GENOMEVCF_BPLUS_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"

/******************
 * Basic Structures
 ******************/

/**
 * @brief  Rank of bplus tree
 * @note   A global variable with "static const" feels like "public static
 * final" in java.
 */
static const int rank_inner_node = 5;
static const int rank_leaf_node = 10;

typedef struct _define_RecVcf_bPlus {
  // TODO
} RecVcf_bPlus;

/**
 * @brief  Use bplus tree to index the vcf records under this chromosome.
 */
typedef struct _define_ChromVcf_bPlus {
  char *name;
  uint32_t recCnt;
  // TODO bplus structure
  struct _define_ChromVcf_bPlus *next;
} ChromVcf_bPlus;

/**
 * @brief  Use linked-list to index all chroms.
 */
typedef struct _define_GenomeVcf_bPlus {
  bcf_hdr_t *hdr;
  uint32_t chromCnt;
  ChromVcf_bPlus *chroms;
} GenomeVcf_bPlus;

GenomeVcf_bPlus *init_GenomeVcf_bPlus();

void destroy_GenomeVcf_bPlus(GenomeVcf_bPlus *gv_bPlus);

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/
void genomeVcf_bPlus_insertRec(GenomeVcf_bPlus *gv_bPlus,
                               RecVcf_bPlus *rv_bPlus);

void genomeVcf_bPlus_removeRec(GenomeVcf_bPlus *gv_bPlus,
                               RecVcf_bPlus *rv_bPlus);

void genomeVcf_bPlus_loadFile(GenomeVcf_bPlus *gv_bPlus, char *filePath);

void genomeVcf_bPlus_writeFile(GenomeVcf_bPlus *gv_bPlus, char *filePath);

/**********************************
 * Debugging Methods for GenomeVcf
 **********************************/

void genomeVcf_bPlus_printRec(GenomeVcf_bPlus *gv_bPlus,
                              RecVcf_bPlus *rv_bPlus);

void _testSet_genomeVcf_bPlus();

#endif