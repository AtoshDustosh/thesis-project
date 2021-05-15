/*
 * @Date: 2021-02-11 11:31:34
 * @LastEditors: AtoshDustosh
 * @LastEditTime: 2021-02-22 13:09:01
 * @FilePath: /Genome_Realignment_Based_on_Variants/main.c
 */

#include <getopt.h>
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "alignment.h"
#include "auxiliaryMethods.h"
#include "genomeFa.h"
#include "genomeSam.h"
#include "genomeVcf.h"
#include "grbvOperations.h"
#include "grbvOptions.h"
#include "simpleOperations.h"

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
    {"extractChrom", required_argument, NULL, OPT_EXTRACTCHROM},

    {"selectBadReads", required_argument, NULL, OPT_SELECTBADREADS},
    {"integrateVcfToSam", no_argument, NULL, OPT_INTEGRATEVCFTOSAM},
    {0, 0, 0, 0},
};

static void Usage() {
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
  printf(
      "\tcountRec\tcount records for all input files. Execute successfully "
      "only when the files' formats are correct\n");
  printf(
      "\tfirstLines [number]\tprint the first [number] lines for all files to "
      "console\n");
  printf(
      "\textractChrom [chrom_idx 1-based]\textract bases of the selected "
      "chromosome "
      "and write into designated output file together with the chromosome's "
      "info field\n");
  printf("\n");

  printf(" -- GRBV operations\n");
  printf(
      "\tselectBadReads [MAPQ_threshold]\tselect reads with MAPQ lower than "
      "MAPQ_threshold from previously set sam file and then output them into "
      "the previously set output file\n");
  printf(
      "\tintegrateVcfToSam\tintegrate variants from *.vcf file with *.sam "
      "file. This will perform realignment for all reads in the *.sam file "
      "with new created reference genome. It's actually one of the main "
      "purposes of the project. \n");
}

static int _testSet_full() {
  // TODO debug session ...
  _testSet_auxiliaryMethods();
  printf("... auxiliary methods test passed. \n");
  _testSet_genomeFa();
  printf("... genomeFa test passed. \n");
  _testSet_genomeSam();
  printf("... genomeSam test passed. \n");
  _testSet_genomeVcf();
  printf("... genomeVcf test passed. \n");
  _testSet_alignment();
  printf("... alignment test passed. \n");
  _testSet_grbvOperations();
  printf("... grbvOperation test passed. \n");
  printf("... all test passed :)\n");
  printf("\n");
  // printf("press \"Enter\" to continue. \n");
  // getchar();
  // ... debug session
  return 1;
}

int main(int argc, char *argv[]) {
  Options options;

  options.ifOptConflict = 0;
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
  while (1) {
    switch (optRet) {
      case OPT_VERBOSE: {
        options.verbose = 1;
        assert(_testSet_full());
        break;
      }
      case OPT_SET_OUTPUTFILE: {
        printf("Output file: %s\n", optarg);
        options.outputFile = optarg;
        break;
      }
      case OPT_SET_FAFILE: {
        printf("Fa/Fna (Reference Genome) file: %s\n", optarg);
        options.faFile = optarg;
        // Should not create gf here.
        // loadGenomeFaFromFile(gf, optarg);
        // printf("... genome data (%s) loaded successfully. \n", optarg);
        // printGenomeFa_brief(gf);
        break;
      }
      case OPT_SET_FASTQFILE: {
        printf("Fastq (Runs) file: %s\n", optarg);
        options.fastqFile = optarg;
        break;
      }
      case OPT_SET_SAMFILE: {
        printf("Sam (alignment) file: %s\n", optarg);
        options.samFile = optarg;
        // Should not create gs here.
        // GenomeSam *gs = init_GenomeSam();
        // loadGenomeSamFromFile(gs, optarg);
        // printf("... genome data (%s) loaded successfully. \n", optarg);
        // printGenomeSam(gs);
        // destroy_GenomeSam(gs);
        break;
      }
      case OPT_SET_VCFFILE: {
        printf("Vcf (variants) file: %s\n", optarg);
        options.vcfFile = optarg;
        // Should not create gv here.
        // loadGenomeVcfFromFile(gv, optarg);
        // printf("... genome data (%s) loaded successfully. \n", optarg);
        // printGenomeVcf(gv);
        break;
      }
      case OPT_COUNTREC: {
        optCheck_conflict(&options);
        printf("Count records of files.\n");
        options.countRec = 1;
        countRec(&options);
        break;
      }
      case OPT_FIRSTLINES: {
        optCheck_conflict(&options);
        printf("Print first %s lines of files.\n", optarg);
        options.firstLines = atoi(optarg);
        firstLines(&options);
        break;
      }
      case OPT_EXTRACTCHROM: {
        optCheck_conflict(&options);
        printf("Extract bases of the #%s chromosome. \n", optarg);
        options.extractChrom = atoi(optarg);
        extractChrom(&options);
        break;
      }
      case OPT_SELECTBADREADS: {
        optCheck_conflict(&options);
        printf("Select bad reads with MAPQ lower than %s\n", optarg);
        options.selectBadReads = atoi(optarg);
        if (options.selectBadReads < 0) {
          printf("Arg invalid: [MAPQ] lower than 0\n");
          exit(EXIT_FAILURE);
        } else if (options.selectBadReads > 255) {
          printf("Arg warning: [MAPQ] higher than max [255]\n");
        } else {
          selectBadReads(&options);
        }
        break;
      }
      case OPT_INTEGRATEVCFTOSAM: {
        optCheck_conflict(&options);
        integrateVcfToSam_refactored(&options);
        break;
      }
      default:
        Usage();
        break;
    }
    optRet = getopt_long(argc, argv, optStr, optInitArray, NULL);
    if (optRet == -1) break;
  }

  return 0;
}