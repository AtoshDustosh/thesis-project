#ifndef COMBINATIONSOFVARS_H_INCLUDED
#define COMBINATIONSOFVARS_H_INCLUDED

#pragma once

#include <inttypes.h>

#include "debug.h"
#include "genomeVcf.h"

typedef struct _define_ElementRecVcf {
  RecVcf *rv;
  /*
   * "alleleIdx" and "alleleCnt" are used to mark which alleles should be
   * integrated and how many of them should be integrated on this vcf record.
   */
  int *alleleIdx;
  int alleleCnt;
} ElementRecVcf;

/**
 * @brief  A structure for basic calculation of combinations. Note that memory
 * for arrays within this structure is assigned dynamically outside this
 * structure. If you malloc an array and used it to initialize this structure,
 * don't free the array using free(). Instead, use this structure's destructuion
 * methods.
 */
typedef struct _define_Combination {
  int **combis;  // combinations
  int combiSize;
  int combiCnt;
} Combinations;

/**
 * @brief  A structure for calculation of the combinations of alleles in RecVcf.
 * Note that memory for arrays within this structure is assigned dynamically
 * outside this structure. If you malloc an array and used it to initialize this
 * structure, don't free the array using free(). Instead, use this structure's
 * destructuion methods.
 */
typedef struct _define_AlleleCombinations {
  // combinations of selected vcf records' indexes
  int *rvCombi;
  // combinations of selected alleles' indexes in vcf records
  int **alleleCombis;
  // number of alleles' in a combination
  int combiSize;
  // count of combinations of selected alleles'
  int combiCnt;
} AlleleCombinations;

Combinations *combinations_init(int **combis, int combiSize, int combiCnt);

void combinations_print(Combinations *cbs);

void combinations_destroy(Combinations *cbs);

AlleleCombinations *init_AlleleCombinations(int *rvCombi, int **alleleCombis,
                                            int combiSize, int combiCnt);

void print_AlleleCombinations(AlleleCombinations *acbs);

void destroy_AlleleCombinations(AlleleCombinations *acbs);

/**
 * @brief  Get all combinations of combiSize from array. This method does not
 * handle duplicated elements in the input array.
 */
Combinations *combinations(int array[], int arraySize, int combiSize);

/**
 * * @brief  Permutate selected vcf records and output all combinations of their
 * alleles into the data structure  "AlleleCombinations". For vcf records with
 * multiple alleles, only 1 allele within the same vcf record can be selected.
 * Some vcf records may have long alleles and these alleles may cover the
 * following alleles. In such cases, if a long allele is selected, the alleles
 * covered by that long allele should not be selected. We call these alleles as
 * incompatible.
 * @note
 * @param  *ervArray[]: array of ElementRecVcf object
 * @param  ervCnt: length of ervArray / count of vcf records
 * @param  *rvCombi: combination of selected vcf records
 * @param  combiSize: size of combination / number of selected vcf records
 * @retval
 */
AlleleCombinations *alleleCombinations(ElementRecVcf *ervArray[], int ervCnt,
                                       int *rvCombi, int combiSize);

void _testSet_combinationsOfVars();

#endif