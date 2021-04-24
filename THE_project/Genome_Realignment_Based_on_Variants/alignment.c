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
  // match > 0, mismatch, gapOpen, gapExtension < 0
  score_match = match > 0 ? match : -match;
  score_mismatch = mismatch < 0 ? mismatch : -mismatch;
  score_gapOpen = gapOpen < 0 ? gapOpen : -gapOpen;
  score_gapExtension = gapExtension < 0 ? gapExtension : -gapExtension;
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
        scoreMat[i * 5 + j] = match;
      } else {
        scoreMat[i * 5 + j] = mismatch;
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
                const int qlen, AlignResult *ar) {
  if (ar == NULL) {
    fprintf(stderr, "Warning: null pointer for AlignResult. Autofixed. \n");
    ar = init_AlignResult();
  }
  // Initialize alignment structures
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));

  // Code original sequences (char*) into matrix (uint8_t*)
  uint8_t *numTseq = (uint8_t *)calloc(tlen, sizeof(uint8_t));
  uint8_t *numQseq = (uint8_t *)calloc(qlen, sizeof(uint8_t));

  for (int i = 0; i < tlen; i++) numTseq[i] = nt_table[(int)tseq[i]];
  for (int i = 0; i < qlen; i++) numQseq[i] = nt_table[(int)qseq[i]];

  // Align
  ksw_extz(0, qlen, numQseq, tlen, numTseq, 5, scoreMat, score_gapOpen,
           score_gapExtension, -1, -1, 0, &ez);

  // Create CIGAR string
  char buf[256];
  memset(buf, 0, sizeof(buf));
  for (int i = 0; i < ksw_extz_cigarCnt(&ez); ++i) {
    char tmpStr[10];
    // There is no itoa() under Linux. Use sprintf instead.
    sprintf(tmpStr, "%d", ksw_extz_cigarLen(&ez, i));
    strcat(buf, tmpStr);
    sprintf(tmpStr, "%c", ksw_extz_cigarOp(&ez, i));
    strcat(buf, tmpStr);
  }
  char *cigar = (char *)calloc(strlen(buf) + 1, sizeof(char));
  strcpy(cigar, buf);

  // Assign results (global alignment)
  ar->pos = 0;
  ar->cigar = cigar;
  ar->mapq = ez.score;
  ar->ref_begin = 0;
  ar->ref_end = tlen - 1;
  ar->read_begin = 0;
  ar->read_end = qlen - 1;

  free(ez.cigar);
  free(numTseq);
  free(numQseq);
}

void align_ssw(const char *tseq, const int tlen, const char *qseq,
               const int qlen, AlignResult *ar) {
  if (ar == NULL) {
    fprintf(stderr, "Warning: null pointer for AlignResult. Autofixed. \n");
    ar = init_AlignResult();
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
  s_align *result = ssw_align(profile, numRef, tlen, score_gapOpen,
                              score_gapExtension, 1, 0, 0, qlen / 2);

  // Create CIGAR string
  char buf[256];
  memset(buf, 0, sizeof(buf));
  for (int i = 0; i < result->cigarLen; i++) {
    char tmpStr[10];
    uint32_t cigarLen = cigar_int_to_len(result->cigar[i]);
    char cigarOp = cigar_int_to_op(result->cigar[i]);
    sprintf(tmpStr, "%" PRIu32 "%c", cigarLen, cigarOp);
    strcat(buf, tmpStr);
  }

  char *cigar = (char *)calloc(strlen(buf) + 1, sizeof(char));
  strcpy(cigar, buf);

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
  free(numRef);
  free(numRead);

  return;
}

static int _test_ksw2Alignment() {
  // align
  static const char *tseq = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
  static const char *qseq = "CTGAGCCGGTAAATC";
  AlignResult *ar = init_AlignResult();
  align_ksw2(tseq, strlen(tseq), qseq, strlen(qseq), ar);
  // printf("ksw2 align\t");
  // print_AlignResult(ar);
  destroy_AlignResult(ar);
  return 1;
}

static int _test_sswAlignment() {
  // align
  static const char *tseq = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
  static const char *qseq = "CTGAGCCGGTAAATC";
  AlignResult *ar = init_AlignResult();
  align_ssw(tseq, strlen(tseq), qseq, strlen(qseq), ar);
  // printf("ssw align\t");
  // print_AlignResult(ar);
  destroy_AlignResult(ar);
  return 1;
}

void _testSet_alignment() {
  // default parameters for genome sequence alignment
  static int32_t match = 4, mismatch = -6;
  static int32_t gapOpen = 3, gapExtension = 1;
  alignInitialize(match, mismatch, gapOpen, gapExtension);
  assert(_test_ksw2Alignment());
  assert(_test_sswAlignment());
}