#include "alignment.h"

/*
 * This table is used to transform nucleotide letters into numbers.
 * Only used for "align_ssw()"
 */
static int8_t nt_table[256];

// Scoring matrix
static int8_t scoreMat[25];

static int ksw2_bandWidth = KSW2_DEFAULT_BANDWIDTH;
static int ksw2_zdrop = KSW2_DEFAULT_ZDROP;
static int ksw2_flag = KSW2_DEFAULT_FLAG;

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

void alignInitialize_ksw2(int bandWidth, int zdrop, int flag) {
  ksw2_bandWidth = bandWidth;
  ksw2_zdrop = zdrop;
  ksw2_flag = flag;
}

void align_ksw2(const char *tseq, const int tlen, const char *qseq,
                const int qlen, AlignResult *ar) {
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
  //          score_gapExtension, ksw2_bandWidth, ksw2_zdrop, 0, &ez);
  ksw_extz2_sse(km, qlen, numQseq, tlen, numTseq, 5, scoreMat, score_gapOpen,
                score_gapExtension, ksw2_bandWidth, ksw2_zdrop, ksw2_flag, 0,
                &ez);
  // This one gets different result with other methods. Don't use this
  // ez.score = ksw_gg2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, 5,
  //                        scoreMat, score_gapOpen, score_gapExtension,
  //                        ksw2_bandWidth,
  //                        &(ez.m_cigar), &(ez.n_cigar), &(ez.cigar));
  // ez.score = ksw_gg2_sse(km, qlen, numQseq, tlen, numTseq, 5, scoreMat,
  //                        score_gapOpen, score_gapExtension, ksw2_bandWidth,
  //                        &(ez.m_cigar), &(ez.n_cigar), &(ez.cigar));

  // Create CIGAR in alignment result
  ar->cnt_cigar = ksw_extz_cigarCnt(&ez);
  ar->cigarLen = (uint32_t *)calloc(ar->cnt_cigar, sizeof(uint32_t));
  ar->cigarOp = (char *)calloc(ar->cnt_cigar, sizeof(char));
  for (int i = 0; i < ar->cnt_cigar; i++) {
    // Copy global cigar operations
    ar->cigarLen[i] = ksw_extz_cigarLen(&ez, i);
    ar->cigarOp[i] = ksw_extz_cigarOp(&ez, i);
  }

  // Assign results (global alignment)
  ar->pos = 0;
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

  // Create CIGAR in alignment result
  if (result->read_begin1 != 0) {  // There exists insetions at the beginning
    ar->cnt_cigar++;
  }
  if (result->read_end1 + 1 != qlen) {  // There exists insertions at the end
    ar->cnt_cigar++;
  }
  int idx_ar_cigar = 0;
  ar->cnt_cigar += result->cigarLen;
  ar->cigarLen = (uint32_t *)calloc(ar->cnt_cigar, sizeof(uint32_t));
  ar->cigarOp = (char *)calloc(ar->cnt_cigar, sizeof(char));
  if (result->read_begin1 != 0) {
    ar->cigarLen[0] = result->read_begin1;
    ar->cigarOp[0] = 'I';
    result->read_begin1 = 0;
    idx_ar_cigar++;
  }
  if (result->read_end1 + 1 != qlen) {
    ar->cigarLen[ar->cnt_cigar - 1] = qlen - result->read_end1 - 1;
    ar->cigarOp[ar->cnt_cigar - 1] = 'I';
    result->read_end1 = qlen - 1;
  }
  for (int i = 0; i < result->cigarLen; i++) {
    ar->cigarLen[idx_ar_cigar] = cigar_int_to_len(result->cigar[i]);
    ar->cigarOp[idx_ar_cigar] = cigar_int_to_op(result->cigar[i]);
    idx_ar_cigar++;
  }

  // Assign results (local alignment)
  ar->pos = result->ref_begin1;
  ar->mapq = result->score1;
  ar->ref_begin = result->ref_begin1;
  ar->ref_end = result->ref_end1;
  ar->read_begin = result->read_begin1;
  ar->read_end = result->read_end1;

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
  align_ksw2(tseq, tlen, qseq, qlen, ar);
  print_AlignResult(ar);
  destroy_AlignResult(ar);
  return 1;
}

static int _test_sswAlignment(const char *tseq, const char *qseq) {
  // printf(" - ssw align\n");
  // printf("tseq: %s\nqseq: %s\n", tseq, qseq);
  AlignResult *ar = init_AlignResult();
  const int tlen = strlen(tseq);
  const int qlen = strlen(qseq);
  align_ssw(tseq, tlen, qseq, qlen, ar);
  // print_AlignResult(ar);
  destroy_AlignResult(ar);
  return 1;
}

void _testSet_alignment() {
  // default parameters for genome sequence alignment
  static int32_t match = 2, mismatch = -2;
  static int32_t gapOpen = -6, gapExtension = -1;
  alignInitialize(match, mismatch, gapOpen, gapExtension);
  alignInitialize_ksw2(KSW2_DEFAULT_BANDWIDTH, KSW2_DEFAULT_ZDROP,
                       KSW2_FLAG_EXTENSION);

  // test cases
  printf("testing alignment ...\n");
  static const char *tseq = "ACTCTAGACGTAATGATTATATAATAAAAAACAAGCTTA";
  static const char *qseq = "ACTCTACCCCGACGTAA";

  assert(_test_ksw2Alignment(tseq, qseq));
  assert(_test_sswAlignment(tseq, qseq));
}
