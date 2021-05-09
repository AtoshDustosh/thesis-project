#include <stdio.h>
#include <stdlib.h>

#include "binarySearch.h"

int main(int argc, char *argv[]) {
  char array_data[] = {'a', 'c', 'g', 't', 'n'};
  int array_key[] = {15, 25, 35, 45, 100};
  int cnt_key = sizeof(array_key) / sizeof(int);

  int searched_keys[] = {2, 15, 16, 21, 25, 35, 36, 40, 45, 100, 124};
  int cnt_searched = sizeof(searched_keys) / sizeof(int);

  int ret_idx = -1;

  for (int i = 0; i < cnt_searched; i++) {
    binarySearch_bplus_inner(array_key, cnt_key, searched_keys[i], &ret_idx);
    printf("query key: %d, ret_idx: %d, query data: %c\n\n", searched_keys[i],
           ret_idx,
           ret_idx != -1 && ret_idx < cnt_key ? array_data[ret_idx] : '*');
    ret_idx = -1;
  }

  return 0;
}