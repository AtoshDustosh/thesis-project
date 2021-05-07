#include <stdio.h>

#include "structA.h"

int main(int argc, char* argv[]){
  structA *strA = init_structA(10);

  print_structA(strA);

  modify_structA(strA, 1);

  print_structA(strA);

  printf("inline method print - data: %d\n", data_structA(strA));

  destroy_structA(strA);

  return 0;
}