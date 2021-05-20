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

/*
 * Be careful with these scores. Different score might result in different
 * results for alignment, especially when you mistake the sign of these scores.
 */
#define SCORE_DEFAULT_MATCH 4
#define SCORE_DEFAULT_MISMATCH -4
#define SCORE_DEFAULT_GAPOPEN -6
#define SCORE_DEFAULT_GAPEXTENSION 1

/*
 * Default parameters for ksw2 alignment.
 */
#define KSW2_DEFAULT_BANDWIDTH -1
#define KSW2_DEFAULT_ZDROP -1
#define KSW2_DEFAULT_FLAG 0

typedef struct _define_AlignResult {
  int64_t pos;
  int cnt_cigar;
  uint32_t *cigarLen;
  char *cigarOp;
  uint8_t mapq;
  int32_t ref_begin;   // 0-based, included
  int32_t ref_end;     // 0-based, included
  int32_t read_begin;  // 0-based, included
  int32_t read_end;    // 0-based, included
} AlignResult;

static inline int64_t arDataPos(AlignResult *ar) { return ar->pos; }
static inline const uint8_t arDataMapQ(AlignResult *ar) { return ar->mapq; }
static inline int32_t arDataRefBegin(AlignResult *ar) { return ar->ref_begin; }
static inline int32_t arDataRefEnd(AlignResult *ar) { return ar->ref_end; }
static inline int32_t arDataReadBegin(AlignResult *ar) {
  return ar->read_begin;
}
static inline int32_t arDataReadEnd(AlignResult *ar) { return ar->read_end; }

static inline int ar_cigar_cnt(AlignResult *ar) { return ar->cnt_cigar; }

static inline char ar_cigarOp(AlignResult *ar, int idx_cigar) {
  if (idx_cigar < ar->cnt_cigar) {
    return ar->cigarOp[idx_cigar];
  } else {
    fprintf(stderr, "Error: array out of bound when extracting cigar op.\n");
    exit(EXIT_FAILURE);
  }
}

static inline uint32_t ar_cigarlen(AlignResult *ar, int idx_cigar) {
  if (idx_cigar < ar->cnt_cigar) {
    return ar->cigarLen[idx_cigar];
  } else {
    fprintf(stderr, "Error: array out of bound when extracting cigar len.\n");
    exit(EXIT_FAILURE);
  }
}

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
 * @brief  Initialize parameters for ksw2 alignment.
 */
void alignInitialize_ksw2(int bandWidth, int zdrop, int flag);

/**
 * @brief  Initialize an AlignResult object. Please note that this object must
 * be freed later using destroy_AlignResult(...).
 */
static inline AlignResult *init_AlignResult() {
  return (AlignResult *)calloc(1, sizeof(AlignResult));
}

static inline void print_AlignResult(AlignResult *ar) {
  printf("Align result:\n");
  if (ar == NULL) return;
  printf("\tpos: %" PRId64 ", mapq: %" PRIu8 "\n", ar->pos, ar->mapq);
  printf("\tref begin: %" PRId32 ", ref end: %" PRId32 "\n", ar->ref_begin,
         ar->ref_end);
  printf("\tread begin: %" PRId32 ", read end: %" PRId32 "\n", ar->read_begin,
         ar->read_end);
  printf("\tcigar: ");
  for (int i = 0; i < ar->cnt_cigar; i++) {
    printf("%" PRIu32 "%c", ar->cigarLen[i], ar->cigarOp[i]);
  }
  printf("\n");
}

static inline void destroy_AlignResult(AlignResult *ar) {
  free(ar->cigarOp);
  free(ar->cigarLen);
  free(ar);
}

/**
 * @brief  Do alignment using ksw2.
 * @note  The following paramters needs using "alignIntialize_ksw2" to set:
 *  bandWidth: band with (< 0 to disable)
 *  zdrop: off-diagonal drop-off to stop extension (positive; <0 to disable)
 *  flag: flag (see KSW_EZ_* macros)
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
