/*
 * @LastEditors: Atosh Dustosh
 */

#include "dataReader.h"



static void glimpseFirstLines(char *filePath, int n_firstlines)
{
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_SUCCESS);
  }

  char tempCh;
  while ((tempCh = fgetc(fp)) != EOF && n_firstlines >= 0)
  {
    if (tempCh == '\n')
      n_firstlines--;
    fprintf(stdout, "%c", tempCh);
  }
  free(fp);
}

static void glimpseCountRecords(char *filePath){
  getFileType(filePath);
}

void glimpseFile(const Options *opts)
{
  switch (opts->glimpse.type)
  {
  case GLIMPSE_NOT_DESIGNATED:
    printf("Warning: operation type of glimpse not designated. \n");
    return;
  case GLIMPSE_FIRST_LINES:
    glimpseFirstLines(opts->inputFilePath, opts->glimpse.n_firstLines);
    break;
  case GLIMPSE_COUNT_RECORDS:
    glimpseCountRecords(opts->inputFilePath);
    break;
  default:;
  }
}

void readBcfFile(const Options *opts, BcfData *bD)
{
  // TODO implement it.
}