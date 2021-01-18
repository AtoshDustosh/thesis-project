#ifndef MYMACROS_H_INCLUDED
#define MYMACROS_H_INCLUDED


#include <time.h>

#include "myMacros.h"


#define EXIT_SUCCESS 0

#define EXIT_FAILURE 1


/*
 * A structure used to keep input arguments.
 *
 * \note useful but not necessary for now
 */
typedef struct _define_Options {
  char *fn_ref;
  int flag;
  int clevel;
  int ignore_sam_err;
  int nreads; // num_reads: limit the output to the first num_reads reads
  int extra_hdr_nuls;
  int benchmark;
  int nthreads;
  int multi_reg;
  char *index;
  int min_shift;
} Options;


#endif // MYMACROS_H_INCLUDED
