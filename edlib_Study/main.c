#include <stdio.h>
#include <string.h>
#include <cstdlib>  // this doesn't affect the compiler

#include "edlib.h"

void test1() {
  char* query = "ACCTCTG";
  char* target = "ACTCTGAAA";
  EdlibAlignResult result =
      edlibAlign(query, 7, target, 9, edlibDefaultAlignConfig());
  if (result.status == EDLIB_STATUS_OK) {
    printf("%d", result.editDistance);
  }
  edlibFreeAlignResult(result);
}

void test2() {
  EdlibAlignResult result =
      edlibAlign("hello", 5, "world!", 6, edlibDefaultAlignConfig());
  if (result.status == EDLIB_STATUS_OK) {
    printf("edit_distance('hello', 'world!') = %d\n", result.editDistance);
  }
  edlibFreeAlignResult(result);
}

void test3() {
  char* seq1 = "AGCTACCGCAT";
  char* seq2 = "ATAGCTAGCTAGCAT";
  EdlibAlignResult result = edlibAlign(
      seq1, strlen(seq1), seq2, strlen(seq2),
      edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
  if (result.status == EDLIB_STATUS_OK) {
    printf("%d\n", result.editDistance);
    printf("%d\n", result.alignmentLength);
    printf("%d\n", result.endLocations[0]);
  }
  char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength,
                                      EDLIB_CIGAR_STANDARD);
  printf("%s\n", cigar);
  free(cigar);
  edlibFreeAlignResult(result);
}

int main() {
  test1();
  test2();
  test3();
}