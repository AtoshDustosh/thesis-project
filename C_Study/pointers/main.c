#include <stdio.h>
#include <stdlib.h>

#include "pointer_conversion.h"

int main(int argc, char *argv[]) {
  structA *str1 = init_structA(11, 12, 13);
  structA *str2 = init_structA(21, 22, 23);
  structB *str3 = init_structB('1', '2');
  int cnt_pointer = 3;
  void **pointers = (void **)calloc(cnt_pointer, sizeof(void *));

  pointers[0] = (void *)str1;
  pointers[1] = (void *)str2;
  pointers[2] = (void *)str3;

  structC *strC = init_structC(pointers, cnt_pointer);

  char *strTypes = (char *)calloc(cnt_pointer, sizeof(char));
  strTypes[0] = 'A';
  strTypes[1] = 'A';
  strTypes[2] = 'B';

  print_structC(strC, strTypes);

  // This also works. But I'm afraid of memory leakage. So let's just use the later method from now on.
  // destroy_structC(strC, strTypes);

  destroy_structC_specifiedType(strC, strTypes);

  print_structC(strC, strTypes);
  
  // Double free
  // destroy_structC_specifiedType(strC, strTypes);

  free(pointers);
  free(strTypes);

  return 0;
}