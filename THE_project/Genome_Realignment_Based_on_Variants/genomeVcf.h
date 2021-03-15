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

/**
 * @brief  This is actually a linked-list with header. 
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
  RecVcf *rv;
  struct _define_ChromVcf *next;
} ChromVcf;

typedef struct _define_GenomeVcf {
  // TODO optimize the structure
  uint32_t chromCnt;
  bcf_hdr_t *hdr;
  ChromVcf *cv;
} GenomeVcf;

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
void printVcfRecord_brief(bcf_hdr_t *hdr, bcf1_t *rec);

void printGenomeVcf(GenomeVcf *gv);

void printChromVcf(GenomeVcf *gv, ChromVcf *cv);

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/

static RecVcf *init_RecVcf();

static void destroy_RecVcf(RecVcf *rv);

static ChromVcf *init_ChromVcf();

static void destroy_ChromVcf(ChromVcf *cv);

/**
 * @brief  Create and initialize a GenomeVcf object.
 *
 * @retval pointer to the GenomeVcf object. Note that it must be freed using
 * destory_GenomeVcf() manually later.
 */
GenomeVcf *init_GenomeVcf();

void destroy_GenomeVcf(GenomeVcf *gv);

void addChromToGenomeVcf(ChromVcf *cv, GenomeVcf *gv);

/**
 * @brief  Add a vcf record into the ChromVcf object.
 * // TODO Please do not add any duplicated vcf records. This method does
 * not provide duplication detection. I didn't found any effective method to
 * judege whether 2 vcf records are the same. And there is actually no necessity
 * to do the check if you are loading data from a vcf file.
 */
void addRecToChromVcf(RecVcf *rv, ChromVcf *cv);

ChromVcf *getChromFromGenomeVcf(char *chromName, GenomeVcf *gv);

/**
 * @brief  Get the vcf record with designated index. 
 * @param  idx: 0-based index/id for the vcf record.
 */
RecVcf *getRecFromChromVcf(uint32_t idx, ChromVcf *cv);

/**
 * @brief  Get the position field of the vcf record
 */
#define getRecVcf_pos(rv) rv->rec->pos

/**
 * @brief  Return a copy of chrom name got from RecVcf object.
 * @retval chrom name. The returned string must be freed manually later. 
 */
char *getRecVcf_chNam(RecVcf *rv, GenomeVcf *gv);

/*******************
 * Methods for users
 *******************/

void loadGenomeVcfFromFile(GenomeVcf *gv, char *filePath);

void writeGenomeVcfIntoFile(GenomeVcf *gv, char *filePath);

#endif