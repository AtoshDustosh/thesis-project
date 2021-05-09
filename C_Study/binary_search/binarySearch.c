#include "binarySearch.h"

void binarySearch_idx(int array_key[], int cnt_key, int query_key,
                      int *ret_idx) {
  int l_idx = 0;
  int r_idx = cnt_key - 1;

  int mid_idx = 0;
  while (l_idx <= r_idx) {
    mid_idx = (l_idx + r_idx) / 2;
    int mid_key = array_key[mid_idx];
    if (mid_key == query_key) {
      *ret_idx = mid_idx;
      return;
    } else if (mid_key > query_key) {
      r_idx = mid_idx - 1;
    } else {
      l_idx = mid_idx + 1;
    }
  }

  // printf("query key: %d, mid idx: %d\n", query_key, mid_idx);

  *ret_idx = -1;
}

void binarySearch_bplus_inner(int array_key[], int cnt_key, int query_key,
                              int *ret_idx) {
  int l_idx = 0;
  int r_idx = cnt_key - 1;

  while (l_idx <= r_idx) {
    int mid_idx = (l_idx + r_idx) / 2;
    int mid_key = array_key[mid_idx];
    if (mid_key == query_key) {
      *ret_idx = mid_idx;
      *ret_idx = *ret_idx + 1;
      return;
    } else if (mid_key > query_key) {
      r_idx = mid_idx - 1;
    } else {
      l_idx = mid_idx + 1;
    }
  }
  *ret_idx = l_idx - 1;
  *ret_idx = *ret_idx + 1;
}

void binarySearch_bplus_leaf(int array_key[], int cnt_key, int query_key,
                             int *ret_idx) {
  int l_idx = 0;
  int r_idx = cnt_key - 1;

  while (l_idx <= r_idx) {
    int mid_idx = (l_idx + r_idx) / 2;
    int mid_key = array_key[mid_idx];
    if (mid_key == query_key) {
      *ret_idx = mid_idx;
      // *ret_idx = *ret_idx + 1;
      return;
    } else if (mid_key > query_key) {
      r_idx = mid_idx - 1;
    } else {
      l_idx = mid_idx + 1;
    }
  }
  *ret_idx = l_idx - 1;
  *ret_idx = *ret_idx + 1;
}