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
#include "integrateVcfToSam.h"
#include "kmerGeneration.h"
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
    {"auxFile", required_argument, NULL, OPT_SET_AUXFILE},

    {"sv_min_len", required_argument, NULL, OPT_SET_SV_MIN_LEN},
    {"sv_max_len", required_argument, NULL, OPT_SET_SV_MAX_LEN},
    {"match", required_argument, NULL, OPT_SET_MATCH},
    {"mismatch", required_argument, NULL, OPT_SET_MISMATCH},
    {"gapOpen", required_argument, NULL, OPT_SET_GAPOPEN},
    {"gapExtension", required_argument, NULL, OPT_SET_GAPEXTENSION},

    {"countRec", no_argument, NULL, OPT_COUNTREC},
    {"firstLines", required_argument, NULL, OPT_FIRSTLINES},
    {"extractChrom", required_argument, NULL, OPT_EXTRACTCHROM},

    {"selectBadReads", required_argument, NULL, OPT_SELECTBADREADS},
    {"integrateVcfToSam", required_argument, NULL, OPT_INTEGRATEVCFTOSAM},
    {"threads", required_argument, NULL, OPT_THREADS},

    {"kmerGeneration", required_argument, NULL, OPT_KMERGENERATION},
    {0, 0, 0, 0},
};

static void Usage() {
  printf("Usage: grbv [commands] [arguments]\n");
  printf("Run one task at a time.\n");
  printf("\n");

  printf("Commands:\n");
  printf(" -- Set files. Do this first!\n");
  printf("\toutputFile [filepath]\tset output file. Default: %s\n",
         default_outputFile);
  printf("\tfaFile [filepath]\tset reference genome file\n");
  printf("\tfastqFile [filepath]\tset fastq file\n");
  printf("\tsamFile [filepath]\tset sam file\n");
  printf("\tvcfFile [filepath]\tset vcf file\n");
  printf(
      "\tauxFile [filepath]\tset auxiliary data file. Some options  require "
      "additional input files: kmerGeneration\n");
  printf(
      "\tsv_min_len [length]\tset minimal length for a SV. Designed for "
      "integration. Default: %d\n",
      default_sv_min_len);
  printf(
      "\tsv_max_len [length]\tset maximal length for a SV. Designed for "
      "integration. Do not set this parameter too big. That may cause the "
      "program running for decades! (combinations of too many variants "
      "generated) Default: %d\n",
      default_sv_max_len);
  printf("\tmatch [score]\tset score for match\n");
  printf("\tmismatch [score]\tset score for mismatch\n");
  printf("\tgapOpen [score]\tset score for gapOpen\n");
  printf("\tgapExtension [score]\tset score for gapExtension\n");
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
      "\tthreads [NUM_threads]\tuse multi-threads methods to run the program. "
      "This only works for integrateVcfToSam.\n");
  printf(
      "\tselectBadReads [MAPQ_threshold]\tselect reads with MAPQ lower than "
      "MAPQ_threshold from previously set sam file and then output them into "
      "the previously set output file\n");
  printf(
      "\tintegrateVcfToSam [integration_strategy]\tintegrate variants from "
      "*.vcf file with *.sam file. This will perform realignment for all reads "
      "in the *.sam file with new created reference genome. It's actually one "
      "of the main purposes of the project. \n");
  printf(
      "\t\t\t[integration_strategy]: [%d] SNP only; [%d] SV only; [%d] SNP and "
      "SV\n",
      _OPT_INTEGRATION_SNPONLY, _OPT_INTEGRATION_SVONLY, _OPT_INTEGRATION_ALL);
  printf(
      "\tkmerGeneration [length_kmer]\tRequested function: extract kmers from "
      "specified intervals on reference genome. And integrate variants during "
      "generation of kmers. Result will be output into specified file.\n");
  printf("\n");
}

static int _testSet_full() {
  // ... debug sector
  _testSet_auxiliaryMethods();
  printf("... auxiliary methods test passed. \n");
  // _testSet_hashTable();
  // printf("... hash table test passed. \n");
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
  _testSet_generateKmers();
  printf("... generateKmers test passed. \n");
  printf("... all test passed :)\n");
  printf("\n");
  // printf("press \"Enter\" to continue. \n");
  // getchar();
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
  options.auxFile = NULL;

  options.sv_min_len = default_sv_min_len;
  options.sv_max_len = default_sv_max_len;
  options.match = SCORE_DEFAULT_MATCH;
  options.mismatch = SCORE_DEFAULT_MISMATCH;
  options.gapOpen = SCORE_DEFAULT_GAPOPEN;
  options.gapExtension = SCORE_DEFAULT_GAPEXTENSION;

  options.countRec = 0;
  options.firstLines = 0;

  options.selectBadReads = 0;
  options.threads = 1;

  options.kmerGeneration = 0;

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
        // clock_t time_start = clock();
        // loadGenomeSamFromFile(gs, optarg);
        // clock_t time_end = clock();
        // printf("... genome data (%s) loaded successfully. Time: %fs\n", optarg,
        //        time_convert_clock2second(time_start, time_end));
        // printGenomeSam_brief(gs);
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
      case OPT_SET_AUXFILE: {
        printf("Auxiliary data file: %s\n", optarg);
        options.auxFile = optarg;
        break;
      }
      case OPT_SET_SV_MIN_LEN: {
        printf("Minimal SV length set as: %s\n", optarg);
        options.sv_min_len = atoi(optarg);
        break;
      }
      case OPT_SET_SV_MAX_LEN: {
        printf("Maximal SV length set as: %s\n", optarg);
        options.sv_max_len = atoi(optarg);
        break;
      }
      case OPT_SET_MATCH: {
        printf("score for match: %s\n", optarg);
        options.match = atoi(optarg);
        break;
      }
      case OPT_SET_MISMATCH: {
        printf("score for mismatch: %s\n", optarg);
        options.mismatch = atoi(optarg);
        break;
      }
      case OPT_SET_GAPOPEN: {
        printf("score for gapOpen: %s\n", optarg);
        options.gapOpen = atoi(optarg);
        break;
      }
      case OPT_SET_GAPEXTENSION: {
        printf("score for gapExtension: %s\n", optarg);
        options.gapExtension = atoi(optarg);
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
      case OPT_THREADS: {
        options.threads = atoi(optarg);
        break;
      }
      case OPT_INTEGRATEVCFTOSAM: {
        optCheck_conflict(&options);
        printf("Selected strategy for integration: ");
        options.integration = atoi(optarg);
        switch (options.integration) {
          case _OPT_INTEGRATION_SNPONLY: {
            printf("SNP only\n");
            break;
          }
          case _OPT_INTEGRATION_SVONLY: {
            printf("SV only\n");
            break;
          }
          case _OPT_INTEGRATION_ALL: {
            printf("all\n");
            break;
          }
          default: {
            fprintf(stderr, "Error: no such strategy for integration.\n");
            exit(EXIT_FAILURE);
          }
        }
        // integrateVcfToSam_refactored(&options);
        integration(&options);
        break;
      }
      case OPT_KMERGENERATION: {
        optCheck_conflict(&options);
        printf("Specified length for generated kmer: %s\n", optarg);
        generateKmers(&options);
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