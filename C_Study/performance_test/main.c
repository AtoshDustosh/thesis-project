#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "performanceTest.h"

void glance_sizeof() {
  int size_uint8 = sizeof(uint8_t);
  int size_int8 = sizeof(int8_t);
  int size_uint32 = sizeof(uint32_t);
  int size_int = sizeof(int);
  int size_long = sizeof(long);
  int size_uint64 = sizeof(uint64_t);
  int size_char = sizeof(char);
  int size_float = sizeof(float);
  int size_double = sizeof(double);
  int size_pointer = sizeof(void *);
  printf("size of uint8_t: %d\n", size_uint8);
  printf("size of int8_t: %d\n", size_int8);
  printf("size of uint32_t: %d\n", size_uint32);
  printf("size of int: %d\n", size_int);
  printf("size of long: %d\n", size_long);
  printf("size of uint64_t: %d\n", size_uint64);
  printf("size of char: %d\n", size_char);
  printf("size of float: %d\n", size_float);
  printf("size of double: %d\n", size_double);
  printf("size of pointer: %d\n", size_pointer);
}

int main(int argc, char *argv[]) {
  test_linkedList_array(50000);
  
  return 0;
}