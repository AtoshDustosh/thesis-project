#include "structA.h"

void new_structA_generated(){
  cnt_strA++;
}

struct structA {
  bool ifModified;
  int data;
  structA *next;
};

structA *init_structA(int data) {
  structA *strA = (structA *)calloc(1, sizeof(structA));
  strA->data = data;
  strA->ifModified = false;
  strA->next = NULL;
  new_structA_generated();
  return strA;
}

void modify_structA(structA *strA, int new_data) {
  static int lastUsed_data;
  printf("last used data: %d, cnt strA: %d\n", lastUsed_data, cnt_strA);
  if (strA->data < 10) {
    lastUsed_data = strA->data;
  }
  strA->data = new_data;
  strA->ifModified = true;
}

void append_structA(structA *strA, int new_data) {
  structA *next = strA->next;
  if (next == NULL) {
    strA->next = init_structA(new_data);
  } else {
    strA->next = init_structA(new_data);
    strA->next->next = next;
  }
}

void print_structA(structA *strA) {
  printf("data: %8d\t(modified?:%s)\n", strA->data,
         strA->ifModified == true ? "true" : "false");
}

void traverse_structA(structA *strA) {
  while (strA->next != NULL) {
    printf("(%d:%s) -> ", strA->data,
           strA->ifModified == true ? "true" : "false");
    strA = strA->next;
  }
  printf("\n");
}

void destroy_structA(structA *strA) { free(strA); }

inline int data_structA(structA *strA) { return strA->data; }