#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "thread_study.h"

int main(int argc, char *argv[]) {
  int cnt_threads = 10;
  int64_t cnt_iteration = 2010010010;

  clock_t time_start = clock();
  int64_t tmp_iteration = 0;
  for (int64_t i = 0; i < cnt_iteration; i++) {
    tmp_iteration++;
  }
  clock_t time_end = clock();

  float time_taken = (float)(time_end - time_start) / CLOCKS_PER_SEC;
  // printf("time start: %lu, time end: %lu\n", time_start, time_end);
  printf("(directly) time taken: %fs\n", time_taken);

  threadTest(cnt_threads, cnt_iteration);

  printf("... exiting the program\n");

  return 0;
}