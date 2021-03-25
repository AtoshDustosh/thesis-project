#ifndef GENOMEVCF_H_INCLUDED
#define GENOMEVCF_H_INCLUDED

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
/*
 * At first I wanted to implement this structure using a B+ tree or multi-index
 * structure, but that was a little complicated and time-taking. And after some
 * discussion with my teacher, I decided to use the linked-list to finish it
 * first. No unit list length! That would make it complicated again.
 */

typedef struct _define_RecVcf {
  bcf1_t *rec;
  struct _define_RecVcf *next;
} RecVcf;

typedef struct _define_ChromVcf {
  /*
   * The records are sorted by their "pos".
   */
  char *name;
  uint32_t recCnt;
  /**
   * @brief A linked list with an empty header. 
   */
  RecVcf *rvs;
  struct _define_ChromVcf *next;
} ChromVcf;

typedef struct _define_GenomeVcf {
  // TODO optimize the structure
  uint32_t chromCnt;
  bcf_hdr_t *hdr;
  /**
   * @brief A linked list (no header).  
   */
  ChromVcf *cvs;
} GenomeVcf;

/**
 * @brief An iterator for ChromVcf. Note that this only iterates vcf records of
 * this specific ChromVcf object.
 */
typedef struct _define_ChromVcfIterator {
  ChromVcf *cv;
  RecVcf *tmpRv;
} ChromVcfIterator;

/**
 * @brief An iterator for GenomeVcf. Note that this iterates all vcf records of
 * all ChromVcf objects, but cannot iterate the vcf records of a specific
 * ChromVcf object. And this iterator has nothing to do with the
 * ChromVcfIterator.
 */
typedef struct _define_GenomeVcfIterator {
  GenomeVcf *gv;
  ChromVcf *tmpCv;
  RecVcf *tmpRv;
} GenomeVcfIterator;

static inline bcf1_t *rvData(RecVcf *rv) { return rv->rec; }

/**
 * @brief  Get the 1-based position of the variant.
 */
static inline int64_t rvDataPos(RecVcf *rv) { return 1 + rv->rec->pos; }

/**
 * @brief Max length of all alleles in a vcf record. For example, for (REF  ALT)
 * = (A  ACGT,C), the result is 3. For (REF  ALT) = (ACGT  A), the result = 0.
 * For (REF  ALT) = (A  C), the result is 1.
 */
uint32_t rvDataMaxVarLength(RecVcf *rv);

/**
 * @brief  Return a copy of chrom name got from RecVcf object.
 * @retval chrom name. The returned string must be freed manually later.
 */
static inline char *rvDataChromName(RecVcf *rv, GenomeVcf *gv) {
  return strdup(bcf_seqname_safe(gv->hdr, rv->rec));
}

static inline char *cvName(ChromVcf *cv) { return cv->name; }

static inline uint32_t cvRecCnt(ChromVcf *cv) { return cv->recCnt; }

static inline uint32_t gvChromCnt(GenomeVcf *gv) { return gv->chromCnt; }

RecVcf *init_RecVcf();

void destroy_RecVcf(RecVcf *rv);

ChromVcf *init_ChromVcf();

void destroy_ChromVcf(ChromVcf *cv);

/**
 * @brief  Create and initialize a GenomeVcf object.
 *
 * @retval pointer to the GenomeVcf object. Note that it must be freed using
 * destory_GenomeVcf() manually later.
 */
GenomeVcf *init_GenomeVcf();

void destroy_GenomeVcf(GenomeVcf *gv);

ChromVcfIterator *init_ChromVcfIterator(ChromVcf *cv);

void destroy_ChromVcfIterator(ChromVcfIterator *cvIt);

GenomeVcfIterator *init_GenomeVcfIterator(GenomeVcf *gv);

void destroy_GenomeVcfIterator(GenomeVcfIterator *gvIt);

/**
 * @brief This iterator returns the next vcf record to be iterated. 
 * @return RecVcf* pointer to the next RecVcf object to be iterated. NULL if there is no records left in the chromosome or the iterator is not initialized with a non-NULL ChromVcf. 
 */
RecVcf *cvItNextRec(ChromVcfIterator *cvIt);

/**
 * @brief  This iterator returns the next chrom to be iterated.
 * @retval pointer to the next ChromVcf object to be iterated; NULL if there
 * is no chroms left or the iterator is not initialized with a non-NULL
 * GenomeVcf.
 */
ChromVcf *gvItNextChrom(GenomeVcfIterator *gvIt);

/**
 * @brief  This iterator return the next vcf record to be iterated.
 * @retval pointer to the next RecVcf object to be iterated. NULL if there is
 * no records left in the temporary chrom, chrom is not selected for iteration
 * (use @viItNextChrom before using this), or the iterator is not intialized
 * with a non-NULL GenomeVcf.
 */
RecVcf *gvItNextRec(GenomeVcfIterator *gvIt);

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/

void addChromToGenomeVcf(ChromVcf *cv, GenomeVcf *gv);

/**
 * @brief  Add a vcf record into the ChromVcf object.
 * // TODO Please do not add any duplicated vcf records. This method does
 * not provide duplication detection. I didn't found any effective method to
 * judege whether 2 vcf records are the same. And there is actually no necessity
 * to do the check if you are loading data from a vcf file.
 */
void addRecToChromVcf(RecVcf *rv, ChromVcf *cv);

ChromVcf *getChromFromGenomeVcfbyName(const char *chromName, GenomeVcf *gv);

/**
 * @brief  Get the vcf record with designated index.
 * @param  idx: 0-based index/id for the vcf record.
 */
RecVcf *getRecFromChromVcfbyIdx(uint32_t idx, ChromVcf *cv);

void loadGenomeVcfFromFile(GenomeVcf *gv, char *filePath);

void writeGenomeVcfIntoFile(GenomeVcf *gv, char *filePath);

/**********************************
 * Debugging Methods for GenomeVcf
 **********************************/

void _testSet_genomeVcf();

/**
 * Print information about the header of a vcf/bcf file.
 */
void printVcfHeader(bcf_hdr_t *hdr);

/**
 * @brief Print detailed information with annotations. Do not use this for
 * printing huge amount of vcf records, but use "printVcfRecord_brief" instead.
 */
void printVcfRecord(bcf1_t *rec);

/**
 * @brief  Print vcf record using standard format of *.vcf files.
 */
void printVcfRecord_brief(GenomeVcf *gv, bcf1_t *rec);

void printGenomeVcf(GenomeVcf *gv);

void printChromVcf(GenomeVcf *gv, ChromVcf *cv);

#endif