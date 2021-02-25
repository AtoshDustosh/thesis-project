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
#include "grbvOperations.h"

const char *optStr = "";
int loptArg = 0;
static struct option optInitArray[] = {
    {"verbose", no_argument, NULL, OPT_VERBOSE},
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

  printf(" -- Program infos\n");
  printf("\tverbose\tverbose mode\n");
  printf("\n");

  printf(" -- Simple operations\n");
  printf("\tcountRec\tcount records in previously set files. Only works on the right file format\n");
  printf("\tfirstLines [number]\tprint the first [number] lines of all files to console\n");
  printf("\n");

  printf(" -- GRBV operations\n");
  printf("\tselectBadReads [MAPQ_threshold]\tselect reads with MAPQ lower than MAPQ_threshold from previously set sam file and then output them into the previously set output file\n");
}

int main(int argc, char *argv[])
{
  Options options;

  options.verbose = 0;

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
    switch (optRet)
    {
    case OPT_VERBOSE:
    {
      options.verbose = 1;
      break;
    }
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
    case OPT_SELECTBADREADS:
    {
      printf("Select bad reads with MAPQ lower than %s\n", optarg);
      options.selectBadReads = atoi(optarg);
      if (options.selectBadReads < 0)
      {
        printf("Arg invalid: [MAPQ] lower than 0\n");
        exit(EXIT_FAILURE);
      }
      else if (options.selectBadReads > 255)
      {
        printf("Arg warning: [MAPQ] higher than max [255]\n");
      }
      else
      {
        selectBadReads(&options);
      }
      break;
    }
    default:
      Usage();
      break;
    }
    optRet = getopt_long(argc, argv, optStr, optInitArray, NULL);
    if (optRet == -1)
      break;
  }
  return 0;
}