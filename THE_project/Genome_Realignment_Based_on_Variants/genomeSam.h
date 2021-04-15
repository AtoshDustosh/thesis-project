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

static inline const char *rsDataQname(RecSam *rs) {
  return bam_get_qname(rsData(rs));
}

static inline const char *rsDataRname(GenomeSam *gs, RecSam *rs) {
  return sam_hdr_tid2name(gs->hdr, rs->rec->core.tid);
}

static inline uint32_t rsDataSeqLength(RecSam *rs) {
  return rs->rec->core.l_qseq;
}

static inline uint8_t rsDataMapQ(RecSam *rs) { return rs->rec->core.qual; }

/**
 * @brief  Get the base sequence of the sam record. Note that the successfully
 * returned value must be freed later manually.
 * @note   This method contains a switch structure. Thus "static inline" is not
 * applied.
 */
char *rsDataSeq(RecSam *rs);

/**
 * @brief  Change the pos, cigar and mapq of a bam1_t object but maitain other
 fields. Then return a copy of the modified bam1_t.
 * @note
 * @param read_begin  0-based position
 * @param read_end    0-based position
 * @param cigarStr    cigar. Length(cigarStr) must equal (read_end - read_begin)
 * @retval new bam1_t object (changed pos, cigar and mapq). The bam1_t struct
 returned by a successful call should be freed via bam_destroy1() when it is no
 longer needed.
 */
static inline bam1_t *bamSetPosCigarMapq(bam1_t *rec, int64_t newPos,
                                         int32_t read_begin, int32_t read_end,
                                         const char *cigarStr,
                                         uint8_t newMapQ) {
  bam1_t *newRec = bam_init1();
  uint32_t *cigarBuf = NULL;
  char *end;
  size_t m = 0;

  const char *qname = bam_get_qname(rec);
  // It seems the l_qname retained using rec->core.l_qname is not correct.
  size_t l_qname = strlen(qname);
  uint16_t flag = rec->core.flag;
  int32_t tid = rec->core.tid;
  hts_pos_t pos = newPos;
  uint8_t mapq = newMapQ;
  size_t n_cigar = sam_parse_cigar(cigarStr, &end, &cigarBuf, &m);
  const uint32_t *cigar = cigarBuf;
  int32_t mtid = rec->core.mtid;
  hts_pos_t mpos = rec->core.mpos;
  hts_pos_t isize = rec->core.isize;
  // This "+1" is for (array index [1]-[0] = 1, array length = 1 + 1)
  size_t l_seq = read_end - read_begin + 1;
  // const char *seq = bam_get_seq(rec);
  // This "+1" is for the '\0' at the end of a string
  char *seq = (char *)calloc(l_seq + 1, sizeof(char));
  for (int i = 0; i < l_seq; i++) {
    switch (bam_seqi(bam_get_seq(rec), read_begin + i)) {
      case 1: {
        seq[i] = 'A';
        break;
      }
      case 2: {
        seq[i] = 'C';
        break;
      }
      case 4: {
        seq[i] = 'G';
        break;
      }
      case 8: {
        seq[i] = 'T';
        break;
      }
      case 15: {
        seq[i] = 'T';
        break;
      }
      default: {
        fprintf(stderr, "Error: unexpected error when decoding bases.\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  const char *qual = bam_get_qual(rec);
  size_t l_aux = bam_get_l_aux(rec);
  bam_set1(newRec, l_qname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid,
           mpos, isize, l_seq, seq, qual, l_aux);
  free(seq);
  return newRec;
}

static inline char *csDataName(ChromSam *cs) { return cs->name; }

static inline uint32_t csDataRecCnt(ChromSam *cs) { return cs->recCnt; }

static inline uint32_t gsDataChromCnt(GenomeSam *gs) { return gs->chromCnt; }

static inline bam_hdr_t *gsDataHdr(GenomeSam *gs) { return gs->hdr; }

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