/*	example.c
 *	This is a simple example to show you how to use the SSW C library.
 *	To run this example:
 *	1) gcc -Wall -lz ssw.c example.c
 *	2) ./a.out
 *	Created by Mengyao Zhao on 07/31/12.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ssw.h"

//	Print the BLAST like output.
static void ssw_write(const s_align* a, const char* ref_seq,
                      const char* read_seq, const int8_t* table) {
  fprintf(stdout,
          "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t",
          a->score1, a->score2);
  if (a->ref_begin1 + 1)
    fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
  fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
  if (a->read_begin1 + 1)
    fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
  fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
  if (a->cigar) {
    int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
    uint32_t i;
    while (e < a->cigarLen || left > 0) {
      int32_t count = 0;
      int32_t q = qb;
      int32_t p = pb;
      fprintf(stdout, "Target: %8d    ", q + 1);
      for (c = e; c < a->cigarLen; ++c) {
        char letter = cigar_int_to_op(a->cigar[c]);
        uint32_t length = cigar_int_to_len(a->cigar[c]);
        uint32_t l = (count == 0 && left > 0) ? left : length;
        for (i = 0; i < l; ++i) {
          if (letter == 'I')
            fprintf(stdout, "-");
          else {
            fprintf(stdout, "%c", *(ref_seq + q));
            ++q;
          }
          ++count;
          if (count == 60) goto step2;
        }
      }
    step2:
      fprintf(stdout, "    %d\n                    ", q);
      q = qb;
      count = 0;
      for (c = e; c < a->cigarLen; ++c) {
        char letter = cigar_int_to_op(a->cigar[c]);
        uint32_t length = cigar_int_to_len(a->cigar[c]);
        uint32_t l = (count == 0 && left > 0) ? left : length;
        for (i = 0; i < l; ++i) {
          if (letter == 'M') {
            if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)])
              fprintf(stdout, "|");
            else
              fprintf(stdout, "*");
            ++q;
            ++p;
          } else {
            fprintf(stdout, "*");
            if (letter == 'I')
              ++p;
            else
              ++q;
          }
          ++count;
          if (count == 60) {
            qb = q;
            goto step3;
          }
        }
      }
    step3:
      p = pb;
      fprintf(stdout, "\nQuery:  %8d    ", p + 1);
      count = 0;
      for (c = e; c < a->cigarLen; ++c) {
        char letter = cigar_int_to_op(a->cigar[c]);
        uint32_t length = cigar_int_to_len(a->cigar[c]);
        uint32_t l = (count == 0 && left > 0) ? left : length;
        for (i = 0; i < l; ++i) {
          if (letter == 'D')
            fprintf(stdout, "-");
          else {
            fprintf(stdout, "%c", *(read_seq + p));
            ++p;
          }
          ++count;
          if (count == 60) {
            pb = p;
            left = l - i - 1;
            e = (left == 0) ? (c + 1) : c;
            goto end;
          }
        }
      }
      e = c;
      left = 0;
    end:
      fprintf(stdout, "    %d\n\n", p);
    }
  }
}

/**
 * @brief  Do Smith-Waterman alignment.
 * @note
 * @param  readSeq: pointer to the query sequence
 * @param  readLen: length of the query sequence
 * @param  refSeq: pointer to the target sequence
 * @param  refLen: length ofthe target sequece
 * @param  match: match score
 * @param  mismatch: mismatch score (negative number for penalty)
 * @param  gapOpen: gap open penalty
 * @param  gapExtension: gap extension penalty
 * * @retval
 */
s_align* align(const char* readSeq, const int readLen, const char* refSeq,
               const int refLen, const int match, const int mismatch,
               const int gapOpen, const int gapExtension) {
  /* This table is used to transform nucleotide letters into numbers. */
  static const int8_t nt_table[128] = {
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
      4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 0, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  s_profile* profile;
  int8_t* numRead = (int8_t*)malloc(readLen + 1);
  int8_t* numRef = (int8_t*)malloc(refLen + 1);
  int8_t* scoreMat = (int8_t*)calloc(25, sizeof(int8_t));
  // initialize scoring matrix for genome sequences, for example:
  //  A  C  G  T	N (or other ambiguous code)
  //  2 -2 -2 -2 	0	A
  // -2  2 -2 -2 	0	C
  // -2 -2  2 -2 	0	G
  // -2 -2 -2  2 	0	T
  //	0  0  0  0  0	N (or other ambiguous code)
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        scoreMat[i * 5 + j] = match;
      } else {
        scoreMat[i * 5 + j] = mismatch;
      }
    }
  }

  for (int i = 0; i < readLen; i++) numRead[i] = nt_table[(int)readSeq[i]];
  for (int i = 0; i < refLen; i++) numRef[i] = nt_table[(int)refSeq[i]];

  profile = ssw_init(numRead, readLen, scoreMat, 5, 2);

  s_align* result = ssw_align(profile, numRef, refLen, gapOpen, gapExtension, 1,
                              0, 0, readLen / 2);

  printf("ref:\t%s\n", refSeq);
  printf("read:\t%s\n", readSeq);
  fprintf(stdout,
          "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t",
          result->score1, result->score2);
  if (result->ref_begin1 + 1)
    fprintf(stdout, "target_begin: %d\t", result->ref_begin1 + 1);
  fprintf(stdout, "target_end: %d\t", result->ref_end1 + 1);
  if (result->read_begin1 + 1)
    fprintf(stdout, "query_begin: %d\t", result->read_begin1 + 1);
  fprintf(stdout, "query_end: %d\n", result->read_end1 + 1);
  for (int i = 0; i < result->cigarLen; i++) {
    uint32_t cigarLen = cigar_int_to_len(result->cigar[i]);
    char cigarOp = cigar_int_to_op(result->cigar[i]);
    printf("%u%c", cigarLen, cigarOp);
  }
  printf("\n");

  init_destroy(profile);
  free(scoreMat);
  free(numRef);
  free(numRead);
  ;

  return result;
}

//	Align a pair of genome sequences.
int main(int argc, char* const argv[]) {
  /**
   * Compile: $gcc main.c ssw.c -o main
   * Run: $./main
   */
  static int32_t match = 2, mismatch = -2, gapOpen = 3,
                 gapExtension =
                     1;  // default parameters for genome sequence alignment
  // reference sequence
  static const char* ref_seq = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
  static const char* read_seq = "CTGAGCCGGTAAATC";

  s_align* result = align(read_seq, strlen(read_seq), ref_seq, strlen(ref_seq),
                          match, mismatch, gapOpen, gapExtension);
  align_destroy(result);

  return (0);
}
