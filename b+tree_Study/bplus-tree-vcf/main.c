#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "genomeVcf_bPlus.h"

int main(int argc, char *argv[]) {
  clock_t time_start = clock();
  GenomeVcf_bplus *gv = genomeVcf_bplus_loadFile("merged.sorted.vcf", 5, 9);
  clock_t time_end = clock();

  printf("performance statistics - record cnt: %d, time: %fs\n", gv_cnt_rec(gv),
         ((float)(time_end - time_start)) / CLOCKS_PER_SEC);

  // genomeVcf_bplus_traverse(gv);

  // 47801942
  // 47800000
  int64_t query_pos = 47801944;
  const char *query_chromName = "21";
  RecVcf_bplus *rv =
      genomeVcf_bplus_getRecAfterPos(gv, query_chromName, query_pos);
  if (rv == NULL) {
    printf("... did not find requested record after pos %" PRId64 "\n",
           query_pos);
  } else {
    printf("... find records after %" PRId64 ":\n", query_pos);
  }
  while (rv != NULL) {
    genomeVcf_bplus_printRec(gv, rv);
    rv = next_RecVcf_bplus(rv);
  }

  destroy_GenomeVcf_bplus(gv);

  printf("press enter to quit the program\n");

  getchar();

  printf("(finished)hello,world\n");

  return 0;
}