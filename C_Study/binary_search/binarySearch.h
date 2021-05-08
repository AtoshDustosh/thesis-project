#ifndef BINARYSEARCH_H_INCLUDED
#define BINARYSEARCH_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>

void binarySearch_idx(int array_key[], int cnt_key, int query_key,
                      int *ret_idx);

void binarySearch_bplus_inner(int array_key[], int cnt_key, int query_key,
                              int *ret_idx);

void binarySearch_bplus_leaf(int array_key[], int cnt_key, int query_key,
                             int *ret_idx);

#endif