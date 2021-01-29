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

static void countSamRecords(char *filePath)
{
  printf("counting sam records ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_SUCCESS);
  }
  int recCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex = 0;
  char tempCh;
  while ((tempCh = fgetc(fp)) != EOF)
  {
    if (tempCh != '\n')
    {
      buffer[bufIndex++] = tempCh;
      continue;
    }
    else
    {
      buffer[bufIndex] = '\0';
      bufIndex = 0;
      if (buffer[0] != '@')
      {
        recCnt++;
      }
    }
  }
  free(fp);

  fprintf(stdout, "Sam records count: %d\n", recCnt);
}

static void countVcfRecords(char *filePath)
{
  printf("counting vcf records ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_SUCCESS);
  }
  int recCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex = 0;
  char tempCh;
  while ((tempCh = fgetc(fp)) != EOF)
  {
    if (tempCh != '\n')
    {
      buffer[bufIndex++] = tempCh;
      continue;
    }
    else
    {
      buffer[bufIndex] = '\0';
      bufIndex = 0;
      if (buffer[0] != '#')
      {
        recCnt++;
      }
    }
  }
  free(fp);

  fprintf(stdout, "Vcf records count: %d\n", recCnt);
}

static void countRefLength(char *filePath)
{
  printf("counting ref length ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_SUCCESS);
  }
  int bpCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex = 0;
  char tempCh;
  while ((tempCh = fgetc(fp)) != EOF)
  {
    if (tempCh != '\n')
    {
      buffer[bufIndex++] = tempCh;
      bpCnt++;
      continue;
    }
    else
    {
      buffer[bufIndex] = '\0';
      bufIndex = 0;
      if (buffer[0] == '>')
      {
        fprintf(stdout, "bps count: %d\n", bpCnt);
        fprintf(stdout, "%s\n", buffer);
        bpCnt = 0;
      }
    }
  }
  free(fp);

  fprintf(stdout, "Ref bps count: %d\n", bpCnt);
}

static void countFastqRecords(char *filePath)
{
  printf("counting fastq records ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_SUCCESS);
  }
  int recCnt = 0;
  int unitLineCnt = 0;
  char tempCh;
  while ((tempCh = fgetc(fp)) != EOF)
  {
    if (tempCh == '\n')
    {
      unitLineCnt++;
    }
    if (unitLineCnt == 4)
    {
      unitLineCnt = 0;
      recCnt++;
    }
  }
  free(fp);

  fprintf(stdout, "Fastq records count: %d\n", recCnt);
}

static void glimpseCountRecords(char *filePath)
{
  switch (getFileType(filePath))
  {
  case FILE_REFERENCE:
    countRefLength(filePath);
    break;
  case FILE_READS:
    countFastqRecords(filePath);
    break;
  case FILE_ALIGNMENTS:
    countSamRecords(filePath);
    break;
  case FILE_VARIANTS:
    countVcfRecords(filePath);
    break;
  default:
    fprintf(stderr, "Warning: unknown file name suffix.\n");
  }
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
  default:
    fprintf(stderr, "Warning: unknown glimpse operation type.\n");
  }
}

void readBcfFile(const Options *opts, BcfData *bD)
{
  // TODO implement it.
}