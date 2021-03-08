#include "genomeMacros.h"

char charOfBase(Base bp) {
  switch (bp) {
    case BASE_A:
      return 'A';
    case BASE_C:
      return 'C';
    case BASE_G:
      return 'G';
    case BASE_T:
      return 'T';
    case BASE_N:
      return 'N';
    case BASE_INVALID:
      return '*';
    default:
      return '*';
  }
}

Base baseOfChar(char bp) {
  switch (bp) {
    case 'A':
      return BASE_A;
    case 'C':
      return BASE_C;
    case 'G':
      return BASE_G;
    case 'T':
      return BASE_T;
    case 'N':
      return BASE_N;
    case '*':
      return BASE_INVALID;
    default:
      return BASE_INVALID;
  }
}