#include "simpleOperations.h"

#define MAX_LINE_BUFFER 2048

static void countRecFa(char *filePath)
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
  fclose(fp);
  fprintf(stdout, "Ref records count: %d\n", bpCnt);
}

static void countRecFastq(char *filePath)
{
  printf("Counting fastq records ...\n");
  FILE* fp  = fopen(filePath, "r");
  if(fp == NULL){
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }
  int recCnt = 0;
  for(char tempCh = fgetc(fp); tempCh != EOF; tempCh = fgetc(fp)){
    if(tempCh == '\n') recCnt++;
  }
  recCnt = recCnt / 4;
  fclose(fp);
  fprintf(stdout, "Fastq records count: %d\n", recCnt);
}

static void countRecSam(char *filePath)
{
  printf("Counting sam records ...\n");
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
  printf("Counting vcf records ...\n");
  FILE* fp = fopen(filePath, "r");
  if(fp == NULL){
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }
  int recCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex = 0;
  for(char tempCh = fgetc(fp); tempCh !=EOF; tempCh = fgetc(fp)){
    if(tempCh != '\n'){
      buffer[bufIndex++] = tempCh;
      continue;
    }else{
      buffer[bufIndex] = '\0';
      bufIndex = 0;
      if(buffer[0] != '#'){
        recCnt++;
      }
    }
  }
  fclose(fp);
  fprintf(stdout, "Vcf records count: %d\n", recCnt);
}

void countRec(Options *opts)
{
  if (getFaFile(opts) != NULL)
    countRecFa(getFaFile(opts));
  if (getFastqFile(opts) != NULL)
    countRecFastq(getFastqFile(opts));
  if (getSamFile(opts) != NULL)
    countRecSam(getSamFile(opts));
  if (getvcfFile(opts) != NULL)
    countRecVcf(getvcfFile(opts));
}

void firstLines(Options *opts)
{
  FileList *fileList = designatedFiles(opts);
  for (int i = 0; i < fileList->count; i++)
  {
    printf(" -- first %d lines of file: %s\n", opts->firstLines, fileList->paths[i]);
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