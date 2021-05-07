#include "pointer_conversion.h"

struct structA {
  int data_1;
  int data_2;
  int data_3;
};

struct structB {
  char data_1;
  char data_2;
};

struct structC {
  void **pointers;
  int cnt_pointer;
};

structA *init_structA(int data1, int data2, int data3) {
  structA *strA = (structA *)malloc(sizeof(structA));
  strA->data_1 = data1;
  strA->data_2 = data2;
  strA->data_3 = data3;
  return strA;
}

structB *init_structB(char data1, char data2) {
  structB *strB = (structB *)malloc(sizeof(structB));
  strB->data_1 = data1;
  strB->data_2 = data2;
  return strB;
}

structC *init_structC(void **pointers, int cnt_pointer) {
  structC *strC = (structC *)malloc(sizeof(structC));
  strC->pointers = (void **)calloc(cnt_pointer, sizeof(void *));
  for (int i = 0; i < cnt_pointer; i++) {
    strC->pointers[i] = pointers[i];
  }
  strC->cnt_pointer = cnt_pointer;
  return strC;
}

void print_structA(structA *str) {
  printf("strA: %d, %d, %d\n", str->data_1, str->data_2, str->data_3);
}

void print_structB(structB *str) {
  printf("strB: %c, %c\n", str->data_1, str->data_2);
}

void print_structC(structC *str, char *strTypes) {
  for (int i = 0; i < str->cnt_pointer; i++) {
    switch (strTypes[i]) {
      case 'A': {
        print_structA((structA *)str->pointers[i]);
        break;
      }
      case 'B': {
        print_structB((structB *)str->pointers[i]);
        break;
      }
      default: {
        fprintf(stderr, "Error: unexpected structure type. \n");
        exit(EXIT_FAILURE);
      }
    }
  }
}

void destroy_structA(structA *strA){
  free(strA);
}

void destroy_structB(structB *strB){
  free(strB);
}

void destroy_structC(structC *strC){
  for(int i = 0; i < strC->cnt_pointer; i++){
    free(strC->pointers[i]);
  }
}

void destroy_structC_specifiedType(structC *str, char *strTypes){
  for (int i = 0; i < str->cnt_pointer; i++) {
    switch (strTypes[i]) {
      case 'A': {
        free((structA *)str->pointers[i]);
        break;
      }
      case 'B': {
        free((structB *)str->pointers[i]);
        break;
      }
      default: {
        fprintf(stderr, "Error: unexpected structure type. \n");
        exit(EXIT_FAILURE);
      }
    }
  }
}