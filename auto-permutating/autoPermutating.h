#ifndef AUTOPERMUTATING_H_INCLUDED
#define AUTOPERMUTATING_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct _define_CombinationResults {
  int **combis;  // compositions
  int combiSize;
  int combiCnt;
} Combinations;

typedef struct _define_Element {
  int recCnt;
  int recPos;
  int *recTypes;
  int *recLengths;
} Element;

typedef struct _define_EleCombinations {
  int *eleCombi;
  int **alleleCombis;
  int combiSize;
  int combiCnt;
} EleCombinations;

Element *element_init(int recPos, int *recTypes, int *recLengths, int recCnt);

void element_print(Element *e);

void element_destroy(Element *e);

Combinations *combinations_init(int **combis, int combiSize, int combiCnt);

void combinations_print(Combinations *cr);

void combinations_destroy(Combinations *cr);

EleCombinations *eleCombinations_init(int *eleCombi, int **alleleCombis,
                                      int combiSize, int combiCnt);

void eleCombinations_print(EleCombinations *ecr);

void eleCombinations_destroy(EleCombinations *ecr);

/**
 * @brief  Get all compositions of combiSize from array. This method does not
 * handle duplicated elements in the input array.
 */
Combinations *combinations(int array[], int arraySize, int combiSize);

/**
 * @brief  Permutate all elements and output all compositions into the data
 * structure  "EleCombinations". For elements with multiple records, only 1
 * record within the same element can be selected.
 * Some elements may have long records and these records may cover the following
 * records. In such cases, if a long record is selected, the record covered by
 * the long record should not be selected.
 * @param *combi combination of elements' indexes in the array of elements
 * @param combiSize size of combinations / number of selected elements in a
 * combination
 * @param *eles[] array of elements
 */
EleCombinations *eleCombinations(int *eleCombi, int combiSize, Element *eles[]);

void _testSet_autoPermutating();

#endif