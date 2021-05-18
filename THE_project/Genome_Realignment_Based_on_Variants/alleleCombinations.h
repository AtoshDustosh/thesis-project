#ifndef ALLELECOMBINATIONS_H_INCLUDED
#define ALLELECOMBINATIONS_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "genomeVcf_bPlus.h"

typedef struct _define_Element_RecVcf {
  RecVcf_bplus *rv;
  /*
   * * "alleleIdx" and "alleleCnt" are used to mark which alleles should be
   * integrated and how many of them should be integrated on this vcf record.
   */
  int *alleleIdx;
  int alleleCnt;
} Element_RecVcf;

/**
 * @brief  A structure for basic calculation of combinations.
 * Note that memory for arrays within this structure is assigned dynamically
 * outside the initialization method. When you need to destroy this object, use
 * destroy_combinations() to free it. Note that this will also free the arrays
 * used to initialize this object.
 */
typedef struct _define_Combinations {
  int **combis;  // combinations
  int length;    // length of a combination
  int cnt;       // count of combinations kept in "combis"
} Combinations;

/**
 * @brief  Similar as Combinations. But this structure is used for calculating
 * combinatins of alleles.
 */
typedef struct _define_Combinations_alleles {
  // Indexes for selected vcf records in a combination
  int *combi_rv;
  // Indexes for selected alleles in selected vcf records
  int **combis_allele;
  // Number of alleles in a combination
  int length;
  // Count of combinations for alleles
  int cnt;
} Combinations_alleles;

Combinations *init_combinations(int **combis, int length, int cnt);

void print_combinations(Combinations *cbs);

void destroy_combinations(Combinations *cbs);

Combinations_alleles *init_combination_alleles(int *combi_rv,
                                               int **combis_allele, int length,
                                               int cnt);

void print_combinations_alleles(Combinations_alleles *acbs);

void destroy_combinations_alleles(Combinations_alleles *acbs);

/**
 * @brief  Get all combinations of combiSize from array. This method does not
 * handle duplicated elements in the input array.
 */
Combinations *calculate_combinations(int array[], int length_array,
                                     int length_combi);

/**
 * @brief  Permutate selected vcf records and output all combinations of their
 * alleles into the data structure  "AlleleCombinations". For vcf records with
 * multiple alleles, only 1 allele within the same vcf record can be selected.
 * Some vcf records may have long alleles and these alleles may cover the
 * following alleles. In such cases, if a long allele is selected, the alleles
 * covered by that long allele should not be selected. We call these alleles as
 * incompatible.
 * @param  *ervArray[]: array of ElementRecVcf object
 * @param  length_array: length of ervArray / count of vcf records
 * @param  *combi_rv: combination of selected vcf records
 * @param  length_combi: length of combination / number of selected vcf records
 */
Combinations_alleles *calculate_combinations_alleles(Element_RecVcf *ervArray[],
                                                     int length_array,
                                                     int *combi_rv,
                                                     int length_combi);

void _testSet_alleleCombinations();

#endif