#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {
  FILE* fp = fopen("input.txt", "r");

  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file. \n");
    exit(EXIT_FAILURE);
  }

  int chrom_id = 0;
  int pos = 0;
  char buf[25];
  char inputChar;
  char outputChar;

  memset(buf, 0, 25);

  int ret = 0;
  while ((ret = fscanf(fp, "[%d %d %s %c %c]\n", &chrom_id, &pos, buf,
                       &inputChar, &outputChar)) == 5) {
    fprintf(stdout, "(ret:%d) - [%d,%d,%s,%c,%c]\n", ret, chrom_id, pos, buf,
            inputChar, outputChar);

    // memset(buf, 0, 25);
  }

  fclose(fp);

  return 0;
}
