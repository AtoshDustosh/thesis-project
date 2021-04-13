#ifndef ALIGNMENT_H_INCLUDED
#define ALIGNMENT_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "debug.h"
#include "ksw2.h"
#include "ssw.h"

#define SCORE_DEFAULT_MATCH 2
#define SCORE_DEFAULT_MISMATCH -2
#define SCORE_DEFAULT_GAPOPEN 3
#define SCORE_DEFAULT_GAPEXTENSION 1

typedef struct _define_AlignResult {
  int64_t pos;
  char *cigar;
  uint8_t mapq;
  int32_t ref_begin;
  int32_t ref_end;
  int32_t read_begin;
  int32_t read_end;
} AlignResult;

static inline int64_t arDataPos(AlignResult *ar) { return ar->pos; }
static inline const char *arDataCigar(AlignResult *ar) { return ar->cigar; }
static inline const uint8_t arDataMapQ(AlignResult *ar) { return ar->mapq; }
static inline int32_t arDataRefBegin(AlignResult *ar) { return ar->ref_begin; }
static inline int32_t arDataRefEnd(AlignResult *ar) { return ar->ref_end; }
static inline int32_t arDataReadBegin(AlignResult *ar) {
  return ar->read_begin;
}
static inline int32_t arDataReadEnd(AlignResult *ar) { return ar->read_end; }
/**
 * @brief  Fix the POS of the alignment result.
 * @param  fixPos: positive or negative integer; 0 is allowed but it's
 * meaningless
 */
static inline void fixPos_AlignResult(AlignResult *ar, int fixPos) {
  ar->pos += fixPos;
}

/**
 * @brief  Initialize the scoring matrix for alignment. Use this method before
 * executing any alignments.
 */
void alignInitialize(int match, int mismatch, int gapOpen, int gapExtension);

/**
 * @brief  Initialize an AlignResult object. Please note that this object must
 * be freed later using destroy_AlignResult(...).
 */
static inline AlignResult *init_AlignResult() {
  return (AlignResult *)calloc(1, sizeof(AlignResult));
}

static inline void print_AlignResult(AlignResult *ar) {
  printf("Align result:\n\t");
  printf("pos: %" PRId64 ", cigar: %s, mapq: %" PRIu8 "\n", ar->pos, ar->cigar,
         ar->mapq);
  printf("ref begin: %" PRId32 ", ref end: %" PRId32 "\n", ar->ref_begin,
         ar->ref_end);
  printf("read begin: %" PRId32 ", ref end: %" PRId32 "\n", ar->read_begin,
         ar->read_end);
}

static inline void destroy_AlignResult(AlignResult *ar) { free(ar); }

/**
 * @brief  Do alignment using ksw2.
 * @note
 * @param  *tseq: target sequence
 * @param  tlen: length of target sequence
 * @param  *qseq: query sequence
 * @param  qlen: length of query seqeunce
 * @param  *ar: pointer to the result of alignment
 * @retval None
 */
void align_ksw2(const char *tseq, const int tlen, const char *qseq,
                const int qlen, AlignResult *ar);

/**
 * @brief  Do alignment using ssw.
 * @note
 * @param  *tseq: target sequence
 * @param  tlen: length of target sequence
 * @param  *qseq: query sequence
 * @param  qlen: length of query seqeunce
 * @param  *ar: pointer to the result of alignment
 * @retval None
 */
void align_ssw(const char *tseq, const int tlen, const char *qseq,
               const int qlen, AlignResult *ar);

void _testSet_alignment();

#endif