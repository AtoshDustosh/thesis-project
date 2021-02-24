/*
 * @Date: 2021-02-11 11:31:34
 * @LastEditors: AtoshDustosh
 * @LastEditTime: 2021-02-22 13:09:01
 * @FilePath: /Genome_Realignment_Based_on_Variants/main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <htslib/vcf.h>

#include "grbvOptions.h"
#include "simpleOperations.h"

const char *optStr = "";
int loptArg = 0;
static struct option optInitArray[] = {

    {"outputFile", required_argument, NULL, OPT_SET_OUTPUTFILE},
    {"faFile", required_argument, NULL, OPT_SET_FAFILE},
    {"fastqFile", required_argument, NULL, OPT_SET_FASTQFILE},
    {"samFile", required_argument, NULL, OPT_SET_SAMFILE},
    {"vcfFile", required_argument, NULL, OPT_SET_VCFFILE},

    {"countRec", no_argument, NULL, OPT_COUNTREC},
    {"firstLines", required_argument, NULL, OPT_FIRSTLINES},

    {"selectBadReads", required_argument, NULL, OPT_SELECTBADREADS},
    {0, 0, 0, 0},
};

static void Usage()
{
  printf("Usage: grbv [commands] [arguments]\n");
  printf("\n");

  printf("Commands:\n");
  printf(" -- Set files. Do this first!\n");
  printf("\toutputFile [filepath]\tset output file\n");
  printf("\tfaFile [filepath]\tset reference genome file\n");
  printf("\tfastqFile [filepath]\tset fastq file\n");
  printf("\tsamFile [filepath]\tset sam file\n");
  printf("\tvcfFile [filepath]\tset vcf file\n");
  printf("\n");

  printf(" -- Simple operations\n");
  printf("\tcountRec\tcount records in previously set filess\n");
  printf("\tfirstLines [number]\tprint the first [number] lines of all files to console\n");
  printf("\n");

  printf(" -- GRBV operations\n");
  printf("\tselectBadReads [MAPQ_threshold]\tselect reads with MAPQ lower than MAPQ_threshold from previously set sam file and then output them into the previously set output file\n");
}

static void printOptions(Options *opts)
{
  printf("faFile: %s\n", opts->faFile);
  printf("fastqFile: %s\n", opts->fastqFile);
  printf("samFile: %s\n", opts->samFile);
  printf("vcfFile: %s\n", opts->vcfFile);
  printf("outputFile: %s\n", opts->outputFile);
  printf("countRec: %d\n", opts->countRec);
  printf("firstLines: %d\n", opts->firstLines);
}

int main(int argc, char *argv[])
{
  Options options;

  options.faFile = NULL;
  options.fastqFile = NULL;
  options.samFile = NULL;
  options.vcfFile = NULL;
  options.outputFile = NULL;

  options.countRec = 0;
  options.firstLines = 0;

  options.selectBadReads = 0;

  int optRet = getopt_long(argc, argv, optStr, optInitArray, NULL);
  while (1)
  {
    printf("optRet: %d\n", optRet);
    switch (optRet)
    {
    case OPT_SET_OUTPUTFILE:
    {
      printf("Output file: %s\n", optarg);
      options.outputFile = optarg;
      break;
    }
    case OPT_SET_FAFILE:
    {
      printf("Fa/Fna (Reference Genome) file: %s\n", optarg);
      options.faFile = optarg;
      break;
    }
    case OPT_SET_FASTQFILE:
    {
      printf("Fastq (Runs) file: %s\n", optarg);
      options.fastqFile = optarg;
      break;
    }
    case OPT_SET_SAMFILE:
    {
      printf("Sam (alignment) file: %s\n", optarg);
      options.samFile = optarg;
      break;
    }
    case OPT_SET_VCFFILE:
    {
      printf("Vcf (variants) file: %s\n", optarg);
      options.vcfFile = optarg;
      break;
    }
    case OPT_COUNTREC:
    {
      printf("Count records of files.\n");
      options.countRec = 1;
      countRec(&options);
      break;
    }
    case OPT_FIRSTLINES:
    {
      printf("Print first %s lines of files.\n", optarg);
      options.firstLines = atoi(optarg);
      firstLines(&options);
      break;
    }
    default:
      Usage();
      break;
    }
    optRet = getopt_long(argc, argv, optStr, optInitArray, NULL);
    if(optRet == -1) break;
  }
  return 0;
}