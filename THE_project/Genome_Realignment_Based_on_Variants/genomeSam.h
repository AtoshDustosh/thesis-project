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
  RecSam *rs;
  struct _define_ChromSam *next;
}ChromSam;

typedef struct _define_GenomeSam {
  // TODO optimize the structure
  uint32_t chromCnt;
  bam_hdr_t *hdr;
  ChromSam *cs;
}GenomeSam;


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
void printSamRecord_brief(bam_hdr_t *hdr,bam1_t *record);

void printGenomeSam(GenomeSam *gs);

void printChromSam(GenomeSam *gs, ChromSam *cs);


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


#endif