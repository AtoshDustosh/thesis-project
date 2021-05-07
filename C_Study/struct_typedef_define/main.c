#include <stdio.h>

#include "structA.h"

int main(int argc, char* argv[]) {
  structA* strA = init_structA(0);

  char *testString = NULL;

  printf("string: %s\n", testString);

  print_structA(strA);

  for (int i = 14; i > 0; i--) {
    modify_structA(strA, i);
    append_structA(strA, i);
    traverse_structA(strA);
  }

  printf("inline method print - data: %d\n", data_structA(strA));

  destroy_structA(strA);

  return 0;
}