#include "alignment.h"

/*
 * This table is used to transform nucleotide letters into numbers.
 * Only used for "align_ssw()"
 */
static int8_t nt_table[256];

// Scoring matrix
static int8_t scoreMat[25];

static int score_match = SCORE_DEFAULT_MATCH;
static int score_mismatch = SCORE_DEFAULT_MISMATCH;
static int score_gapOpen = SCORE_DEFAULT_GAPOPEN;
static int score_gapExtension = SCORE_DEFAULT_GAPEXTENSION;

void alignInitialize(int match, int mismatch, int gapOpen, int gapExtension) {
  // Inappropriate scores will result in odd cigars (especially when you get
  // confused on whether gapOpen and gapExtension should be positive or
  // negative). match, gapOpen, gapExtension > 0, mismatch < 0
  score_match = match > 0 ? match : -match;
  score_mismatch = mismatch < 0 ? mismatch : -mismatch;
  score_gapOpen = gapOpen > 0 ? gapOpen : -gapOpen;
  score_gapExtension = gapExtension > 0 ? gapExtension : -gapExtension;
  // initialize scoring matrix for genome sequences. For example,
  // when match = 2, and mismatch = -2, the matrix is:
  //  A  C  G  T	N (or other ambiguous code)
  //  2 -2 -2 -2 	0	A
  // -2  2 -2 -2 	0	C
  // -2 -2  2 -2 	0	G
  // -2 -2 -2  2 	0	T
  //	0  0  0  0  0	N (or other ambiguous code)
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        scoreMat[i * 5 + j] = score_match;
      } else {
        scoreMat[i * 5 + j] = score_mismatch;
      }
    }
  }
  for (int i = 0; i < 5; i++) {
    scoreMat[i * 5 - 1] = 0;  // the last column
    scoreMat[20 + i] = 0;     // the last row
  }
  // Initialize the encoding matrix for bases
  memset(nt_table, 4, 256);
  nt_table['A'] = nt_table['a'] = 0;
  nt_table['C'] = nt_table['c'] = 1;
  nt_table['G'] = nt_table['g'] = 2;
  nt_table['T'] = nt_table['t'] = 3;
}

void align_ksw2(const char *tseq, const int tlen, const char *qseq,
                const int qlen, int bandWidth, int zdrop, int flag,
                AlignResult *ar) {
  // printf("target seq(%d): %s\n", tlen, tseq);
  // printf("query seq(%d): %s\n", qlen, qseq);
  if (ar == NULL) {
    fprintf(stderr, "Error: null pointer for AlignResult. \n");
    exit(EXIT_FAILURE);
  }

  // Code original sequences (char*) into matrix (uint8_t*)
  uint8_t *numTseq = (uint8_t *)calloc(tlen, sizeof(uint8_t));
  uint8_t *numQseq = (uint8_t *)calloc(qlen, sizeof(uint8_t));

  for (int i = 0; i < tlen; i++) numTseq[i] = nt_table[(uint8_t)tseq[i]];
  for (int i = 0; i < qlen; i++) numQseq[i] = nt_table[(uint8_t)qseq[i]];

  // Initialize alignment structures
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  void *km = 0;

  // Align
  // ksw_extz(km, qlen, numQseq, tlen, numTseq, 5, scoreMat, score_gapOpen,
  //          score_gapExtension, bandWidth, zdrop, 0, &ez);
  ksw_extz2_sse(km, qlen, numQseq, tlen, numTseq, 5, scoreMat, score_gapOpen,
                score_gapExtension, bandWidth, zdrop, flag, 0, &ez);
  // ez->score = ksw_gg2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m,
  // mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
  // ez.score = ksw_gg2_sse(km, qlen, numQseq, tlen, numTseq, 5, scoreMat,
  //                        score_gapOpen, score_gapExtension, bandWidth,
  //                        &ez.m_cigar, &ez.n_cigar, &ez.cigar);

  // Create CIGAR string
  static char globalCigar[256];
  memset(globalCigar, 0, sizeof(globalCigar));
  int cigar_cnt = ksw_extz_cigarCnt(&ez);
  for (int i = 0; i < cigar_cnt; ++i) {
    char cigar_i[10];
    int cigarLen = ksw_extz_cigarLen(&ez, i);
    char cigarOp = ksw_extz_cigarOp(&ez, i);
    // Copy global cigar operations
    // There is no itoa() under Linux. Use sprintf instead.
    sprintf(cigar_i, "%d%c", cigarLen, cigarOp);
    strcat(globalCigar, cigar_i);
  }
  char *cigar = (char *)calloc(strlen(globalCigar) + 1, sizeof(char));
  strcpy(cigar, globalCigar);

  // Assign results (global alignment)
  ar->pos = 0;
  ar->cigar = cigar;
  ar->mapq = ez.score;
  ar->ref_begin = 0;
  ar->ref_end = tlen - 1;
  ar->read_begin = 0;
  ar->read_end = qlen - 1;

  kfree(km, ez.cigar);
  free(numTseq);
  free(numQseq);
}

void align_ssw(const char *tseq, const int tlen, const char *qseq,
               const int qlen, AlignResult *ar) {
  // printf("target seq(%d): %s\n", tlen, tseq);
  // printf("query seq(%d): %s\n", qlen, qseq);
  if (ar == NULL) {
    fprintf(stderr, "Error: null pointer for AlignResult. \n");
    exit(EXIT_FAILURE);
  }
  // Code original sequences (char*) into matrix (uint8_t*)
  int8_t *numRead = (int8_t *)malloc(qlen + 1);
  int8_t *numRef = (int8_t *)malloc(tlen + 1);

  for (int i = 0; i < qlen; i++) numRead[i] = nt_table[(int)qseq[i]];
  for (int i = 0; i < tlen; i++) numRef[i] = nt_table[(int)tseq[i]];

  // Initialize alignment structures
  s_profile *profile;
  profile = ssw_init(numRead, qlen, scoreMat, 5, 2);

  // Align
  // See instructions for ssw_align about the value of "maskLen"
  int maskLen = qlen / 2 >= 15 ? qlen / 2 : 15;
  s_align *result = ssw_align(profile, numRef, tlen, score_gapOpen,
                              score_gapExtension, 1, 0, 0, maskLen);

  // Create CIGAR string
  char localCigar[256];
  memset(localCigar, 0, sizeof(localCigar));
  if (result->read_begin1 != 0) {  // There exists insertions at the beginning
    sprintf(localCigar, "%" PRId32 "I", result->read_begin1);
    result->read_begin1 = 0;
  }
  for (int i = 0; i < result->cigarLen; i++) {
    char cigar_i[10];
    uint32_t cigarLen = cigar_int_to_len(result->cigar[i]);
    char cigarOp = cigar_int_to_op(result->cigar[i]);
    sprintf(cigar_i, "%" PRIu32 "%c", cigarLen, cigarOp);
    strcat(localCigar, cigar_i);
  }
  if (result->read_end1 + 1 != qlen) {  // There exists insertions at the end
    char insertion[10];
    sprintf(insertion, "%dI", qlen - result->read_end1 - 1);
    strcat(localCigar, insertion);
    result->read_end1 = qlen - 1;
  }

  char *cigar = (char *)calloc(strlen(localCigar) + 1, sizeof(char));
  strcpy(cigar, localCigar);

  // Assign results (local alignment)
  ar->pos = result->ref_begin1;
  ar->mapq = result->score1;
  ar->cigar = cigar;
  ar->ref_begin = result->ref_begin1;
  ar->ref_end = result->ref_end1;
  ar->read_begin = result->read_begin1;
  ar->read_end = result->read_end1;

  // printf("ref:\t%s\n", tseq);
  // printf("read:\t%s\n", qseq);
  // fprintf(stdout,
  //         "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t",
  //         result->score1, result->score2);
  // if (result->ref_begin1 + 1)
  //   fprintf(stdout, "target_begin: %d\t", result->ref_begin1 + 1);
  // fprintf(stdout, "target_end: %d\t", result->ref_end1 + 1);
  // if (result->read_begin1 + 1)
  //   fprintf(stdout, "query_begin: %d\t", result->read_begin1 + 1);
  // fprintf(stdout, "query_end: %d\n", result->read_end1 + 1);
  // for (int i = 0; i < result->cigarLen; i++) {
  //   uint32_t cigarLen = cigar_int_to_len(result->cigar[i]);
  //   char cigarOp = cigar_int_to_op(result->cigar[i]);
  //   printf("%u%c", cigarLen, cigarOp);
  // }
  // printf("\n");

  init_destroy(profile);
  align_destroy(result);
  free(numRef);
  free(numRead);

  return;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/************************* Debug Methods ************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

static int _test_ksw2Alignment(const char *tseq, const char *qseq) {
  printf(" - ksw2 align\n");
  printf("tseq: %s\nqseq: %s\n", tseq, qseq);
  AlignResult *ar = init_AlignResult();
  const int tlen = strlen(tseq);
  const int qlen = strlen(qseq);
  align_ksw2(tseq, tlen, qseq, qlen, KSW2_DEFAULT_BANDWIDTH, KSW2_DEFAULT_ZDROP,
             KSW2_DEFAULT_FLAG, ar);
  print_AlignResult(ar);
  destroy_AlignResult(ar);
  return 1;
}

static int _test_sswAlignment(const char *tseq, const char *qseq) {
  printf(" - ssw align\n");
  printf("tseq: %s\nqseq: %s\n", tseq, qseq);
  AlignResult *ar = init_AlignResult();
  const int tlen = strlen(tseq);
  const int qlen = strlen(qseq);
  align_ssw(tseq, tlen, qseq, qlen, ar);
  print_AlignResult(ar);
  destroy_AlignResult(ar);
  return 1;
}

void _testSet_alignment() {
  // default parameters for genome sequence alignment
  static int32_t match = 2, mismatch = -4;
  static int32_t gapOpen = 3, gapExtension = 1;
  alignInitialize(match, mismatch, gapOpen, gapExtension);

  // test cases
  static const char *tseq = "AAAAAAAAACGTACGTACGTAAAAACCCCCGTGTGA";
  static const char *qseq = "TTTTACGTACGTACCCCCGTAAA";

  assert(_test_ksw2Alignment(tseq, qseq));
  assert(_test_sswAlignment(tseq, qseq));
}