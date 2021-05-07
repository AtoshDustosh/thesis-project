#include "structA.h"

struct structA {
  bool ifModified;
  int data;
};

structA *init_structA(int data) {
  structA *strA = (structA *)calloc(1, sizeof(structA));
  strA->data = data;
  strA->ifModified = false;
}

void modify_structA(structA *strA, int new_data) {
  strA->data = new_data;
  strA->ifModified = true;
}

void print_structA(structA *strA) {
  printf("data: %8d\t(modified?:%s)\n", strA->data,
         strA->ifModified == true ? "true" : "false");
}

void destroy_structA(structA *strA) { free(strA); }

inline int data_structA(structA *strA) { return strA->data; }