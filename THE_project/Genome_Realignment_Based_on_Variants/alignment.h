#ifndef ALIGNMENT_H_INCLUDED
#define ALIGNMENT_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "debug.h"
#include "ksw2.h"

typedef struct _define_AlignResult {
  char *cigar;
  uint8_t score;
} AlignResult;

/**
 * @brief  Initialize an AlignResult object. Please note that this object must
 * be freed later using destroy_AlignResult(...).
 */
static inline AlignResult *init_AlignResult() {
  return (AlignResult *)calloc(1, sizeof(AlignResult));
}

static inline void destroy_AlignResult(AlignResult *ar) { free(ar); }

void align(const char *tseq, const char *qseq, AlignResult *ar);

void _testSet_alignment();

#endif