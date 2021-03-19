#include "alignment.h"

static void _align_ksw2(const char *tseq, const char *qseq, int sc_mch,
                        int sc_mis, int gapo, int gape, AlignResult *ar) {
  int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis;  // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  int tl = strlen(tseq), ql = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;

  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0;
  c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2;
  c['T'] = c['t'] = 3;  // build the encoding table
  ts = (uint8_t *)malloc(tl);
  qs = (uint8_t *)malloc(ql);
  for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]];  // encode to 0/1/2/3
  for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
  ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
  /*
   * Just copy the codes above when using this.
   */

  char *buf = (char *)calloc(2048, sizeof(char));
  memset(buf, 0, sizeof(buf));
  for (i = 0; i < ez.n_cigar; ++i) {  // print CIGAR
    char tmpStr[10];
    // There is no itoa() under Linux. Use sprintf instead.
    sprintf(tmpStr, "%d", ez.cigar[i] >> 4);
    strcat(buf, tmpStr);
    sprintf(tmpStr, "%c", "MID"[ez.cigar[i] & 0xf]);
    strcat(buf, tmpStr);
  }
  // putchar('\n');
  // printf("score: %d\n", ez.score);

  // printf("cigar length: %ld\n", strlen(buf));
  char *cigar = (char *)calloc(strlen(buf), sizeof(char));
  strcpy(cigar, buf);

  ar->cigar = cigar;
  ar->score = ez.score;

  free(buf);
  free(ez.cigar);
  free(ts);
  free(qs);
}

AlignResult *align(const char *tseq, const char *qseq) {
  AlignResult *ar = init_AlignResult();
  _align_ksw2(tseq, qseq, 1, -2, 2, 1, ar);
  return ar;
}

static int _test_ksw2Alignment() {
  AlignResult *ar = align(
      "CGAAACTGGGCTACTCCATGACCAGGGGCAAAATAGGCTTTTAGCCGCTGCGTTCTGGGAGCTCCTCCCCCT"
      "TCTGGGAGCTCCTCCCCCTCCCCAGAAGGCCAAGGGATGTGGGGGCTGGGGGACTGGGAGGCCTGGCAGTCT"
      "T",
      "CGAAACTGGGCTACTCCATGACCAGGGGCAAAATAGGCTTTTAGCCGCTGCGTTCTGGGAGCTCCTCCCCCT"
      "CCCCAGAAGGCCAAGGGATGTTGGGG");
  destroy_AlignResult(ar);
  return 1;
}

void _testSet_alignment() { assert(_test_ksw2Alignment()); }