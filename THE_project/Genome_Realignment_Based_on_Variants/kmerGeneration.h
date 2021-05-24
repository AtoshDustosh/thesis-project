#ifndef KMERGENERATION_H_INCLUDED
#define KMERGENERATION_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <inttypes.h>

#include "alignment.h"
#include "alleleCombinations.h"
#include "auxiliaryMethods.h"
#include "debug.h"
#include "genomeFa.h"
#include "genomeVcf_bPlus.h"
#include "grbvOptions.h"
#include "kmerHashTable.h"

#define MAX_LENGTH_AUXLINE 1024

static const char ignored_bases[] = {'N', 'M', 'R'};
static const int cnt_ignored_bases = sizeof(ignored_bases) / sizeof(char);

static inline void usage_auxFile() {
  printf("// ---------------------- format of auxFile ---------------------\n");
  printf("Format of interval: \"[#(id_chrom),#(pos_start),(#pos_end)]\"\n");
  printf("#(id_chrom): 1-based\n");
  printf("#(pos_start) and #(pos_end): 1-based.\n");
  printf("Each line stores a record for an interval. For example:\n");
  printf("[1,123211242,123214242]\n");
  printf("[1,123224242,123227242]\n");
  printf("[1,123247242,123249242]\n");
  printf("[5,444249242,444249242]\n");
  printf("\n");
  printf(
      "It is recommended to sort the intervals by their id_chrom and "
      "pos_start\n");
  printf("Output format of kmer: \"[#(pos_start),string_kmer]\"\n");
  printf("\n");
}

/**
 * @brief  Generate kmers.
 */
void generateKmers(Options *opts);

void _testSet_generateKmers();

#endif