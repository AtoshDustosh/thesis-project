#ifndef GENOMESAM_H_INCLUDED
#define GENOMESAM_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"

/******************
 * Basic Structures
 ******************/

typedef struct _define_RecSam {
  bam1_t *rec;
  struct _define_RecSam *next;
} RecSam;

typedef struct _define_ChromSam {
  /*
   * The records are sorted by their "pos".
   */
  char *name;
  uint32_t recCnt;
  RecSam *rss;
  struct _define_ChromSam *next;
} ChromSam;

typedef struct _define_GenomeSam {
  // TODO optimize the structure
  uint32_t chromCnt;
  bam_hdr_t *hdr;
  ChromSam *css;
} GenomeSam;

typedef struct _define_GenomeSamIterator {
  GenomeSam *gs;
  ChromSam *tmpCs;
  RecSam *tmpRs;
} GenomeSamIterator;

static inline bam1_t *rsData(RecSam *rs) { return rs->rec; }

/**
 * @brief  Get the 1-based position of a mapped read.
 */
static inline uint64_t rsDataPos(RecSam *rs) { return 1 + rs->rec->core.pos; }

static inline const char *rsDataRname(GenomeSam *gs, RecSam *rs) {
  return sam_hdr_tid2name(gs->hdr, rs->rec->core.tid);
}

static inline uint32_t rsDataSeqLength(RecSam *rs) {
  return rs->rec->core.l_qseq;
}

/**
 * @brief  Get the base sequence of the sam record. Note that the successfully
 * returned value must be freed later manually.
 * @note   This method contains a switch structure. Thus "static inline" is not
 * applied.
 */
char *rsDataSeq(RecSam *rs);

static inline char *csName(ChromSam *cs) { return cs->name; }

static inline uint32_t csRecCnt(ChromSam *cs) { return cs->recCnt; }

static inline uint32_t gsChromCnt(GenomeSam *gs) { return gs->chromCnt; }

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/

static RecSam *init_RecSam();

static void destroy_RecSam(RecSam *rs);

static ChromSam *init_ChromSam();

static void destroy_ChromSam(ChromSam *cs);

/**
 * @brief  Create and initialize a GenomeSam object.
 *
 * @retval pointer to the GenomeSam object. Note that it must be freed using
 * destroy_GenomeSam() manually later.
 */
GenomeSam *init_GenomeSam();

void destroy_GenomeSam(GenomeSam *gs);

GenomeSamIterator *init_GenomeSamIterator(GenomeSam *gs);

void destroy_GenomeSamIterator(GenomeSamIterator *gsIt);

/**
 * @brief  This iterator return the next chrom to be iterated.
 * @retval pointer to the next ChromSam object to be iterated; NULL if there
 * is no chroms left or the iterator is not initialized with a non-NULL
 * GenomeSam.
 */
ChromSam *gsItNextChrom(GenomeSamIterator *gsIt);

/**
 * @brief  This iterator return the next sam record to be iterated.
 * @retval pointer to the next RecSam object to be iterated. NULL if there is
 * no records left in the temporary chrom, chrom is not selected for iteration
 * (use @vsItNextChrom before using this), or the iterator is not intialized
 * with a non-NULL GenomeSam.
 */
RecSam *gsItNextRec(GenomeSamIterator *gsIt);

void addChromToGenomeSam(ChromSam *cs, GenomeSam *gs);

void addRecToChromSam(RecSam *rs, ChromSam *cs);

ChromSam *getChromFromGenomeSam(char *chromName, GenomeSam *gs);

/**
 * @brief  Get the sam record with designated index.
 * @param  idx: 0-based index/id for the sam record.
 */
RecSam *getRecFromChromSam(uint32_t idx, ChromSam *cs);

// TODO get methods for fields of a sam record
/**
 * @brief  Get the position field of the sam record
 */
#define getRecSam_pos(rs) rs->rec->core.pos

/**
 * @brief  Return a copy of chrom name got from RecSam object.
 * @retval chrom name. The returned string must be freed manually later.
 */
char *getRecSam_chNam(RecSam *rs, GenomeSam *gs);

/*******************
 * Methods for users
 *******************/

void loadGenomeSamFromFile(GenomeSam *gs, char *filePath);

void writeGenomeSamIntoFile(GenomeSam *gs, char *filePath);

/**********************************
 * Debugging Methods for GenomeVcf
 **********************************/

void _testSet_genomeSam();

/**
 * Print information about the header of a sam/bam file.
 */
void printSamHeader(bam_hdr_t *header);

/**
 * @brief Print detailed information with annotations. Do not use this for
 * printing huge amount of sam records, but use "printSamRecord_brief" instead.
 */
void printSamRecord(bam1_t *record);

/**
 * @brief  Print sam record using standard format of *.sam files.
 */
void printSamRecord_brief(GenomeSam *gs, bam1_t *record);

void printGenomeSam(GenomeSam *gs);

void printChromSam(GenomeSam *gs, ChromSam *cs);

#endif