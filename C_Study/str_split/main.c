#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mySplit.h"

static void _test_strsep() {
  /*
   * strsep() will take in the string and seprate it with all given chars from
   * the delim. But it will not treat delim as a whole. Instead, it uses every
   * unique type of character in it to split the string.
   *
   */
  char *str = ">This is a test for splitting strings. \n";
  char *delim = " ->";
  char *token;

  printf("str: %s\n", str);

  char *copiedStr = strdup(str);
  while ((token = strsep(&copiedStr, delim)) != NULL) {
    printf("token: \"%s\"\n", token);
  }
}

void _testSet() {
  _test_strsep();
  _testSet_mySplit();
}

int main() {
  // _testSet();

  char buf[10];
  memset(buf, 0, 10);

  strcat(buf, "Abc");
  printf("buf: %s\n", buf);

  // strcat(buf, NULL); // doesn't work
  printf("buf: %s\n", buf);

  return 0;
}
