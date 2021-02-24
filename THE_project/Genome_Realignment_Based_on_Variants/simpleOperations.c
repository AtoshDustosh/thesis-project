#include "simpleOperations.h"

#define MAX_LINE_BUFFER 2048

static void countRecFa(char *filePath)
{
}

static void countRecFastq(char *filePath)
{
}

static void countRecSam(char *filePath)
{
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }

  int recCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex;
  for (char tempCh = fgetc(fp); tempCh != EOF; tempCh = fgetc(fp))
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

  fclose(fp);

  fprintf(stdout, "Sam records count: %d\n", recCnt);
}

static void countRecVcf(char *filePath)
{
}

void countRec(Options *opts)
{
  if (opts->faFile != NULL)
    countRecFa(opts->faFile);
  if (opts->fastqFile != NULL)
    countRecFastq(opts->fastqFile);
  if (opts->samFile != NULL)
    countRecSam(opts->samFile);
  if (opts->vcfFile != NULL)
    countRecVcf(opts->vcfFile);
}

void firstLines(Options *opts)
{
  FileList *fileList = designatedFiles(opts);
  for (int i = 0; i < fileList->count; i++)
  {
    printf("First %d lines of file: %s\n", opts->firstLines, fileList->paths[i]);
    FILE *fp = fopen(fileList->paths[i], "r");
    if (fp == NULL)
    {
      fprintf(stderr, "Error: failed to open file %s\n", fileList->paths[i]);
      continue;
    }

    char tempCh;
    int n_firstLines = opts->firstLines;
    while ((tempCh = fgetc(fp)) != EOF && n_firstLines >= 0)
    {
      if (tempCh == '\n')
        n_firstLines--;
      fprintf(stdout, "%c", tempCh);
    }
    fclose(fp);
    printf("\n");
  }
  destroyFileList(fileList);
}