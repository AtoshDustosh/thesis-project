#include "simpleOperations.h"

#define MAX_LINE_BUFFER 4096

static void countRecFa(const char *filePath) {
  printf("counting ref length ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_SUCCESS);
  }
  int bpCnt = 0;
  char buf[MAX_LINE_BUFFER];
  int bufIndex = 0;
  char tmpCh;
  while ((tmpCh = fgetc(fp)) != EOF) {
    if (tmpCh != '\n') {
      buf[bufIndex++] = tmpCh;
      bpCnt++;
      continue;
    } else {
      buf[bufIndex] = '\0';
      bufIndex = 0;
      if (buf[0] == '>') {
        fprintf(stdout, "bps count: %d\n", bpCnt);
        fprintf(stdout, "%s\n", buf);
        bpCnt = 0;
      }
    }
  }
  fclose(fp);
  fprintf(stdout, "Ref records count: %d\n", bpCnt);
}

static void countRecFastq(const char *filePath) {
  printf("Counting fastq records ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }
  int recCnt = 0;
  for (char tempCh = fgetc(fp); tempCh != EOF; tempCh = fgetc(fp)) {
    if (tempCh == '\n') recCnt++;
  }
  recCnt = recCnt / 4;
  fclose(fp);
  fprintf(stdout, "Fastq records count: %d\n", recCnt);
}

static void countRecSam(const char *filePath) {
  printf("Counting sam records ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }
  int recCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex;
  for (char tempCh = fgetc(fp); tempCh != EOF; tempCh = fgetc(fp)) {
    if (tempCh != '\n') {
      buffer[bufIndex++] = tempCh;
      continue;
    } else {
      buffer[bufIndex] = '\0';
      bufIndex = 0;
      if (buffer[0] != '@') {
        recCnt++;
      }
    }
  }
  fclose(fp);
  fprintf(stdout, "Sam records count: %d\n", recCnt);
}

static void countRecVcf(const char *filePath) {
  printf("Counting vcf records ...\n");
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }
  int recCnt = 0;
  char buffer[MAX_LINE_BUFFER];
  int bufIndex = 0;
  for (char tempCh = fgetc(fp); tempCh != EOF; tempCh = fgetc(fp)) {
    if (tempCh != '\n') {
      buffer[bufIndex++] = tempCh;
      continue;
    } else {
      buffer[bufIndex] = '\0';
      bufIndex = 0;
      if (buffer[0] != '#') {
        recCnt++;
      }
    }
  }
  fclose(fp);
  fprintf(stdout, "Vcf records count: %d\n", recCnt);
}

void countRec(Options *opts) {
  if (getFaFile(opts) != NULL) countRecFa(getFaFile(opts));
  if (getFastqFile(opts) != NULL) countRecFastq(getFastqFile(opts));
  if (getSamFile(opts) != NULL) countRecSam(getSamFile(opts));
  if (getVcfFile(opts) != NULL) countRecVcf(getVcfFile(opts));
}

void firstLines(Options *opts) {
  FileList *fileList = designatedFiles(opts);
  for (int i = 0; i < fileList->count; i++) {
    printf(" -- first %d lines of file: %s\n", opts->firstLines,
           fileList->paths[i]);
    FILE *fp = fopen(fileList->paths[i], "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: failed to open file %s\n", fileList->paths[i]);
      continue;
    }

    char tempCh;
    int n_firstLines = opts->firstLines;
    while ((tempCh = fgetc(fp)) != EOF && n_firstLines >= 0) {
      if (tempCh == '\n') n_firstLines--;
      fprintf(stdout, "%c", tempCh);
    }
    fclose(fp);
    printf("\n");
  }
  destroyFileList(fileList);
}

void extractChrom(Options *opts) {
  // Check arguments
  if (opts->extractChrom <= 0) {
    fprintf(stderr, "Error: invalid index %d for a chromosome.\n",
            opts->extractChrom);
    exit(EXIT_FAILURE);
  }
  if (getFaFile(opts) == NULL) {
    fprintf(stderr, "Error: reference genome file (*.fa) not designated. \n");
    exit(EXIT_FAILURE);
  }
  if (getOutputFile(opts) == NULL) {
    fprintf(stderr, "Error: output file not designated. \n");
    exit(EXIT_FAILURE);
  }

  // Extract the chromosome's bases
  FILE *fp_fa = fopen(getFaFile(opts), "r");
  FILE *fp_op = fopen(getOutputFile(opts), "w");

  int foundChrom = 0;
  int chromIdx = 0;

  char buf[MAX_LINE_BUFFER];  // line buffer
  int bufIdx = 0;
  char tmpCh;
  // Locate the chromosome and fprint its info line
  while ((tmpCh = fgetc(fp_fa)) != EOF) {
    if (tmpCh != '\n') {
      buf[bufIdx++] = tmpCh;
      continue;
    } else {
      buf[bufIdx] = '\0';
      bufIdx = 0;
      if (buf[0] == '>') {
        chromIdx++;
        if (chromIdx == opts->extractChrom) {
          fprintf(fp_op, "%s\n", buf);
          foundChrom = 1;
          break;
        }
      }
    }
  }
  // Fprintf the chromosome's bases
  while ((tmpCh = fgetc(fp_fa)) != EOF) {
    if (tmpCh != '\n') {
      buf[bufIdx++] = tmpCh;
      continue;
    } else {
      buf[bufIdx] = '\0';
      bufIdx = 0;
      if (buf[0] != '>') {
        fprintf(fp_op, "%s\n", buf);
        continue;
      } else {
        break;
      }
    }
  }

  if (foundChrom == 0) {
    fprintf(stderr,
            "Warning: chrom #%d not found. Probably out of index boundary.\n",
            opts->extractChrom);
  }

  fclose(fp_fa);
  fclose(fp_op);
}